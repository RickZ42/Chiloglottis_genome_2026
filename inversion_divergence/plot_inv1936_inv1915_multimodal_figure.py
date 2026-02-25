#!/usr/bin/env python3
import argparse
import csv
import math
import os
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch
from PIL import Image, ImageChops


DEFAULT_INV_IDS = ["INV1936", "INV1915"]


def parse_args():
    p = argparse.ArgumentParser(description="Build publication-style composite inversion figure for INV1936 and INV1915")
    p.add_argument("--inv-tsv", default="Output/20260127Genome/H1/inversion_rna_minimap2/inv_rna_expr_minimap2_all.selected_inversions.for_hic.tsv")
    p.add_argument("--inv-ids", default=",".join(DEFAULT_INV_IDS))
    p.add_argument("--hic-summary", default="Output/20260127Genome/hic_inv_context_joint_INV1936_INV1915/hic_inv_context.summary.tsv")
    p.add_argument("--hic-png-dir", default="Output/20260127Genome/hic_inv_context_joint_INV1936_INV1915/per_inversion_png")
    p.add_argument("--rna-genes", default="Output/20260127Genome/H1/inversion_rna_minimap2/inv_rna_expr_minimap2_all.genes_expression.tsv")
    p.add_argument("--anchor-pairs", default="Output/20260127Genome/H1H2_inversion_divergence/selected_inversions_quick/INV1936_INV1915/selected_inversions.anchor_divergence_pairs.tsv")
    p.add_argument("--gene-counts", default="Output/20260127Genome/H1H2_inversion_divergence/selected_inversions_quick/INV1936_INV1915/selected_inversions.gene_counts.tsv")
    p.add_argument("--h1-gtf", default="compare_H1_vs_Ophrys/New_H1_H2_ophrys_Arobx/H1t1.gtf")
    p.add_argument("--h2-gtf", default="compare_H1_vs_Ophrys/New_H1_H2_ophrys_Arobx/H2t1.gtf")
    p.add_argument("--h1-func", default="Output/20260127Genome/functional_annotation/H1.Functional_Annotation_results.txt")
    p.add_argument("--h2-func", default="Output/20260127Genome/functional_annotation/H2.Functional_Annotation_results.txt")
    p.add_argument("--outdir", default="Output/20260127Genome/inversion_pubfig_INV1936_INV1915")
    p.add_argument("--prefix", default="INV1936_INV1915_multimodal_pubfig")
    p.add_argument("--gene-flank-max", type=int, default=500000)
    p.add_argument("--gene-flank-min", type=int, default=250000)
    return p.parse_args()


def clean_str(x):
    if x is None:
        return ""
    s = str(x)
    if s == "nan":
        return ""
    s = s.strip()
    if len(s) >= 2 and s[0] == "'" and s[-1] == "'":
        s = s[1:-1]
    if len(s) >= 2 and s[0] == '"' and s[-1] == '"':
        s = s[1:-1]
    return s


def load_inversions(path):
    df = pd.read_csv(path, sep='\t', dtype={'RefChr': str, 'QryChr': str})
    if 'ID' not in df.columns and 'INV_ID' in df.columns:
        df = df.rename(columns={'INV_ID': 'ID'})
    need = {'ID', 'RefChr', 'RefStart', 'RefEnd', 'QryChr', 'QryStart', 'QryEnd'}
    miss = need - set(df.columns)
    if miss:
        raise ValueError(f"Missing columns in inversion TSV: {sorted(miss)}")
    for c in ['RefStart', 'RefEnd', 'QryStart', 'QryEnd']:
        df[c] = pd.to_numeric(df[c])
    return df


def load_hic_summary(path):
    df = pd.read_csv(path, sep='\t')
    return {r['ID']: r for _, r in df.iterrows()}


def parse_gtf_transcripts(gtf_path, chroms):
    rows = []
    n_tx_total = 0
    with open(gtf_path) as fh:
        for ln in fh:
            if not ln or ln.startswith('#'):
                continue
            parts = ln.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, phase, attrs = parts
            if chrom not in chroms:
                continue
            if feature != 'transcript':
                continue
            n_tx_total += 1
            tid = None
            gid = None
            a = attrs.strip()
            m = re.search(r'transcript_id\s+"([^"]+)"', a)
            if m:
                tid = m.group(1)
            if tid is None:
                m = re.search(r'ID=([^;\s]+)', a)
                if m:
                    tid = m.group(1)
            if tid is None:
                token = a.split(';', 1)[0].strip()
                if token:
                    tid = token.split()[0]
            m = re.search(r'gene_id\s+"([^"]+)"', a)
            if m:
                gid = m.group(1)
            if gid is None:
                m = re.search(r'geneID=([^;\s]+)', a)
                if m:
                    gid = m.group(1)
            if tid is None:
                continue
            if gid is None:
                gid = re.sub(r'\.t\d+$', '', tid)
            rows.append({
                'Chr': chrom,
                'Start': int(start),
                'End': int(end),
                'Strand': strand,
                'TranscriptID': tid,
                'GeneID': gid,
                'Is_t1': bool(re.search(r'\.t1$', tid)),
                'Is_tx': True,
            })
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    return df


def subset_window(df, chrom, start, end):
    if df.empty:
        return df.copy()
    m = (df['Chr'] == chrom) & (df['End'] >= start) & (df['Start'] <= end)
    return df.loc[m].copy().sort_values(['Start', 'End', 'TranscriptID'])


def pack_lanes(starts, ends, pad_bp=5000):
    lane_ends = []
    lanes = []
    for s, e in zip(starts, ends):
        assigned = None
        for i, le in enumerate(lane_ends):
            if s > le + pad_bp:
                assigned = i
                lane_ends[i] = e
                break
        if assigned is None:
            assigned = len(lane_ends)
            lane_ends.append(e)
        lanes.append(assigned)
    return lanes


def load_expression(path, inv_ids):
    df = pd.read_csv(path, sep='\t')
    df = df[df['INV_ID'].isin(inv_ids)].copy()
    sample_cols = [c for c in df.columns if c.endswith('.bam')]
    df[sample_cols] = df[sample_cols].apply(pd.to_numeric, errors='coerce').fillna(0.0)
    # collapse potential duplicates by GeneID within inversion (should generally not happen)
    agg_first = ['RefChr', 'RefStart', 'RefEnd', 'QryChr', 'QryStart', 'QryEnd', 'GeneChr', 'GeneStart', 'GeneEnd']
    agg = {c: 'first' for c in agg_first}
    agg.update({c: 'mean' for c in sample_cols})
    df = df.groupby(['INV_ID', 'GeneID'], as_index=False).agg(agg)
    df['MeanExpr'] = df[sample_cols].mean(axis=1)
    df['Log2MeanExpr'] = np.log2(df['MeanExpr'] + 1.0)
    df['Gene_t1'] = df['GeneID'].astype(str) + '.t1'
    sample_names = [Path(c).name.replace('.sorted.bam', '') for c in sample_cols]
    rename_map = {c: s for c, s in zip(sample_cols, sample_names)}
    df = df.rename(columns=rename_map)
    sample_cols = sample_names
    out = {}
    for inv in inv_ids:
        sub = df[df['INV_ID'] == inv].copy().sort_values(['GeneStart', 'GeneEnd', 'GeneID'])
        out[inv] = sub.reset_index(drop=True)
    return out, sample_cols


def load_anchor_pairs(path, inv_ids):
    df = pd.read_csv(path, sep='\t')
    df = df[df['INV_ID'].isin(inv_ids)].copy()
    # sanitize values
    for c in ['aa_identity', 'dN_dS_NG86']:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors='coerce')
    return df


def load_gene_counts(path):
    if not os.path.exists(path):
        return {}
    df = pd.read_csv(path, sep='\t')
    return {r['INV_ID']: r for _, r in df.iterrows()}


def load_function_map(path):
    candidates = defaultdict(list)
    with open(path, newline='') as fh:
        rd = csv.DictReader(fh, delimiter='\t')
        for r in rd:
            ident = clean_str(r.get('IDENTIFIER', ''))
            if not ident or not ident.endswith('.t1'):
                continue
            root = ident[:-3]
            name = clean_str(r.get('NAME', ''))
            desc = clean_str(r.get('DESCRIPTION', ''))
            typ = clean_str(r.get('TYPE', ''))
            bincode = clean_str(r.get('BINCODE', ''))
            if not name and not desc:
                continue
            candidates[root].append({
                'tid': ident,
                'name': name,
                'desc': desc,
                'type': typ,
                'bincode': bincode,
            })

    def choose_best(rows):
        def concise_name(name):
            if not name:
                return ''
            parts = [p for p in name.split('.') if p]
            if not parts:
                return name
            if len(parts) >= 2:
                tail = '.'.join(parts[-2:])
            else:
                tail = parts[-1]
            tail = tail.replace('*( ', '*(')
            return tail

        scored = []
        for r in rows:
            name = r['name']
            desc = r['desc']
            low = (name + ' ' + desc).lower()
            unchar = ('uncharacterised' in low) or ('no annotation' in low and not name)
            pref_t = 0 if r['type'] == 'T' else 1
            cname = concise_name(name) if name else ''
            if not cname and desc:
                cname = desc.split('&')[0].strip()
            cname = cname.replace('mercator4v8.0:', '').strip()
            cname = re.sub(r'\s+', ' ', cname)
            if not cname:
                cname = 'Unannotated'
            scored.append((pref_t, unchar, len(cname), cname, r))
        scored.sort(key=lambda x: (x[0], x[1], x[2], x[3]))
        _, _, _, cname, best = scored[0]
        return {
            'label': cname,
            'bincode': best['bincode'],
            'name': best['name'],
            'desc': best['desc'],
            'type': best['type'],
        }

    return {gene: choose_best(rows) for gene, rows in candidates.items()}


def trim_white_margin(img_path, pad=12):
    im = Image.open(img_path).convert('RGB')
    bg = Image.new('RGB', im.size, (255, 255, 255))
    diff = ImageChops.difference(im, bg)
    bbox = diff.getbbox()
    if bbox is None:
        return np.asarray(im)
    l, t, r, b = bbox
    l = max(0, l - pad)
    t = max(0, t - pad)
    r = min(im.size[0], r + pad)
    b = min(im.size[1], b + pad)
    return np.asarray(im.crop((l, t, r, b)))


def norm_pos(pos, start, end):
    if end <= start:
        return 0.5
    return (pos - start) / float(end - start)


def short_num_bp(n):
    n = int(n)
    if n >= 1_000_000:
        return f"{n/1e6:.2f} Mb"
    if n >= 1_000:
        return f"{n/1e3:.1f} kb"
    return f"{n} bp"


def truncate_text(s, max_len=52):
    s = s or 'Unannotated'
    s = re.sub(r'\s+', ' ', str(s)).strip()
    return s if len(s) <= max_len else s[:max_len - 1] + '…'


def build_gene_window(inv_row, flank_min, flank_max):
    ref_len = abs(int(inv_row['RefEnd']) - int(inv_row['RefStart'])) + 1
    qry_len = abs(int(inv_row['QryEnd']) - int(inv_row['QryStart'])) + 1
    min_len = min(ref_len, qry_len)
    flank = int(min(flank_max, max(flank_min, round(min_len * 0.2))))
    return {
        'H1_start': max(1, int(inv_row['RefStart']) - flank),
        'H1_end': int(inv_row['RefEnd']) + flank,
        'H2_start': max(1, int(inv_row['QryStart']) - flank),
        'H2_end': int(inv_row['QryEnd']) + flank,
        'flank': flank,
    }


def add_panel_title(ax, text):
    ax.text(0, 1.02, text, transform=ax.transAxes, ha='left', va='bottom', fontsize=9, fontweight='bold')


def _draw_axis_ticks_text(ax, window_start, window_end, y, where='below', n=5, fs=6.2):
    ticks = np.linspace(0, 1, n)
    labs = [f"{(window_start + t*(window_end-window_start))/1e6:.2f}" for t in ticks]
    if where == 'below':
        va = 'top'
        ytxt = y
    else:
        va = 'bottom'
        ytxt = y
    for t, lab in zip(ticks, labs):
        ax.text(t, ytxt, lab, ha='center', va=va, fontsize=fs, color='#444')
    ax.text(0.995, ytxt, ' Mb', ha='left', va=va, fontsize=fs, color='#444')


def draw_overview_track(ax, transcripts, window_start, window_end, inv_start, inv_end,
                        hap_label, color_tick, color_density, y0, y1):
    ax.add_patch(patches.Rectangle((0, y0), 1, y1-y0, facecolor='white', edgecolor='#d9d9d9', lw=0.6, zorder=0))
    x_inv1 = norm_pos(inv_start, window_start, window_end)
    x_inv2 = norm_pos(inv_end, window_start, window_end)
    xlo, xhi = sorted([x_inv1, x_inv2])
    ax.add_patch(patches.Rectangle((xlo, y0), max(1e-6, xhi-xlo), y1-y0,
                                   facecolor='#fde8e7', edgecolor='none', alpha=0.95, zorder=0.2))
    for xv in (xlo, xhi):
        ax.plot([xv, xv], [y0, y1], color='#cf3f3f', lw=0.8, ls='--', zorder=1)

    ax.text(-0.015, (y0+y1)/2, hap_label, ha='right', va='center', fontsize=8, fontweight='bold', transform=ax.transAxes)
    if transcripts.empty:
        ax.text(0.5, (y0+y1)/2, 'No t1 transcripts', ha='center', va='center', fontsize=7, color='#777')
        return

    mids = ((transcripts['Start'].to_numpy(dtype=float) + transcripts['End'].to_numpy(dtype=float)) / 2.0)
    xs = np.clip((mids - window_start) / float(max(1, window_end - window_start)), 0, 1)

    # density curve (genes/bin)
    bins = 40
    hist, edges = np.histogram(xs, bins=bins, range=(0, 1))
    if hist.max() > 0:
        hnorm = hist / hist.max()
        xx = (edges[:-1] + edges[1:]) / 2.0
        base = y0 + 0.10*(y1-y0)
        amp = 0.38*(y1-y0)
        yy = base + hnorm * amp
        ax.fill_between(xx, base, yy, color=color_density, alpha=0.22, lw=0, zorder=1.1)
        ax.plot(xx, yy, color=color_density, lw=1.0, zorder=1.2)

    # gene ticks
    tick_y0 = y0 + 0.53*(y1-y0)
    tick_y1 = y0 + 0.90*(y1-y0)
    for x in xs:
        ax.plot([x, x], [tick_y0, tick_y1], color=color_tick, lw=0.55, alpha=0.65, zorder=2)

    ax.text(0.005, y0 + 0.05*(y1-y0), f"{len(transcripts)} t1 transcripts in local window", fontsize=6.0, color='#555', ha='left', va='bottom')


def draw_zoom_gene_track(ax, transcripts, window_start, window_end, inv_start, inv_end,
                         overlap_roots, anchor_roots, hap_label, color_main, color_anchor, color_bg,
                         y0, y1, label_roots=None, label_above=True):
    ax.add_patch(patches.Rectangle((0, y0), 1, y1-y0, facecolor='white', edgecolor='#d9d9d9', lw=0.6, zorder=0))
    x_inv1 = norm_pos(inv_start, window_start, window_end)
    x_inv2 = norm_pos(inv_end, window_start, window_end)
    xlo, xhi = sorted([x_inv1, x_inv2])
    ax.add_patch(patches.Rectangle((xlo, y0), max(1e-6, xhi-xlo), y1-y0, facecolor='#fde8e7', edgecolor='none', alpha=0.98, zorder=0.2))
    for xv in (xlo, xhi):
        ax.plot([xv, xv], [y0, y1], color='#cf3f3f', lw=0.8, ls='--', zorder=1.5)

    ax.text(-0.015, (y0+y1)/2, hap_label, ha='right', va='center', fontsize=8, fontweight='bold', transform=ax.transAxes)
    if transcripts.empty:
        ax.text(0.5, (y0+y1)/2, 'No t1 transcripts in zoom window', ha='center', va='center', fontsize=7, color='#777')
        return {}

    tx = transcripts.sort_values(['Start', 'End', 'TranscriptID']).copy()
    tx['lane'] = pack_lanes(tx['Start'].tolist(), tx['End'].tolist(), pad_bp=max(200, int((window_end-window_start)*0.006)))
    n_lanes = int(tx['lane'].max()) + 1 if len(tx) else 1
    pad_top = 0.03
    pad_bot = 0.03
    band_h = (y1 - y0) - (pad_top + pad_bot)
    lane_h = max(0.018, band_h / max(1, n_lanes))
    draw_h = min(0.028, lane_h * 0.65)

    centers = {}
    label_roots = set(label_roots or [])
    for _, r in tx.iterrows():
        x1 = norm_pos(r['Start'], window_start, window_end)
        x2 = norm_pos(r['End'], window_start, window_end)
        x1, x2 = sorted([x1, x2])
        lane = int(r['lane'])
        y = y0 + pad_bot + lane * lane_h + (lane_h - draw_h)/2
        gene = r['GeneID']
        is_overlap = gene in overlap_roots
        is_anchor = gene in anchor_roots
        fc = color_bg
        ec = '#b5b5b5'
        alpha = 0.55
        lw = 0.35
        z = 2.0
        if is_overlap:
            fc = color_main
            ec = color_main
            alpha = 0.9
            lw = 0.6
            z = 3.0
        if is_anchor:
            fc = color_anchor
            ec = color_anchor
            alpha = 0.98
            lw = 0.7
            z = 3.2
        width = max(0.0012, x2 - x1)
        if width > 0.009:
            head = min(width * 0.28, 0.010)
            arr_x = x1 if r['Strand'] == '+' else x2
            arr_dx = width if r['Strand'] == '+' else -width
            ax.add_patch(patches.FancyArrow(
                arr_x, y + draw_h/2, arr_dx, 0,
                width=draw_h*0.62, head_width=draw_h*0.95, head_length=head,
                length_includes_head=True, facecolor=fc, edgecolor=ec, linewidth=lw, alpha=alpha, zorder=z
            ))
        else:
            ax.add_patch(patches.Rectangle((x1, y), width, draw_h, facecolor=fc, edgecolor=ec, lw=lw, alpha=alpha, zorder=z))

        centers[gene] = ((x1 + x2) / 2.0, y + draw_h/2)

    # Label selected genes (paper-style: highlight key genes only)
    if label_roots:
        sel = [r for _, r in tx.iterrows() if r['GeneID'] in label_roots]
        sel = sorted(sel, key=lambda r: (r['Start'], r['End'], r['GeneID']))
        for idx, r in enumerate(sel):
            gene = r['GeneID']
            cx, cy = centers.get(gene, (None, None))
            if cx is None:
                continue
            y_text = (y1 + 0.014 + (idx % 2) * 0.016) if label_above else (y0 - 0.014 - (idx % 2) * 0.016)
            va = 'bottom' if label_above else 'top'
            ax.plot([cx, cx], [cy, y1 if label_above else y0], color='#777', lw=0.4, alpha=0.8, zorder=4)
            ax.text(cx, y_text, f"{gene}.t1", ha='center', va=va, fontsize=5.6,
                    color='#1f1f1f', rotation=45 if len(sel) > 8 else 0, zorder=4)

    return centers


def draw_anchor_ribbons(ax, pairs_df, h1_pos, h2_pos, h1w1, h1w2, h2w1, h2w2, y_top, y_bottom):
    ax.add_patch(patches.Rectangle((0, y_bottom), 1, y_top-y_bottom, facecolor='white', edgecolor='none', zorder=0))
    if pairs_df.empty:
        return 0
    n = 0
    for _, r in pairs_df.iterrows():
        h1_root = re.sub(r'\.t\d+$', '', str(r['H1_gene']))
        h2_root = re.sub(r'\.t\d+$', '', str(r['H2_gene']))
        if h1_root not in h1_pos or h2_root not in h2_pos:
            continue
        n += 1
        h1r = h1_pos[h1_root]
        h2r = h2_pos[h2_root]
        x1 = norm_pos((h1r['Start'] + h1r['End'])/2, h1w1, h1w2)
        x2 = norm_pos((h2r['Start'] + h2r['End'])/2, h2w1, h2w2)
        same_strand = (h1r['Strand'] == h2r['Strand'])
        aa_id = pd.to_numeric(r.get('aa_identity', np.nan), errors='coerce') if hasattr(pd, 'to_numeric') else np.nan
        lw = 0.9 if not np.isfinite(aa_id) else (0.8 + 1.6 * max(0.0, min(1.0, float(aa_id) - 0.97) / 0.03))
        color = '#2a9d8f' if same_strand else '#e76f51'
        verts = [
            (x1, y_top),
            (x1, (y_top+y_bottom)/2 + 0.04),
            (x2, (y_top+y_bottom)/2 - 0.04),
            (x2, y_bottom),
        ]
        codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4]
        ax.add_patch(PathPatch(MplPath(verts, codes), facecolor='none', edgecolor=color, lw=lw, alpha=0.7, zorder=1.5))
    return n


def draw_synteny_panel(ax, inv_id, inv_row, gene_windows, h1_tx_all, h2_tx_all, h1_tx_t1, h2_tx_t1, pairs_df):
    ref_chr = inv_row['RefChr']
    qry_chr = inv_row['QryChr']
    h1w1, h1w2 = gene_windows['H1_start'], gene_windows['H1_end']
    h2w1, h2w2 = gene_windows['H2_start'], gene_windows['H2_end']
    h1i1, h1i2 = int(inv_row['RefStart']), int(inv_row['RefEnd'])
    h2i1, h2i2 = int(inv_row['QryStart']), int(inv_row['QryEnd'])

    h1_win_all = subset_window(h1_tx_all, ref_chr, h1w1, h1w2)
    h2_win_all = subset_window(h2_tx_all, qry_chr, h2w1, h2w2)
    h1_win = subset_window(h1_tx_t1, ref_chr, h1w1, h1w2)
    h2_win = subset_window(h2_tx_t1, qry_chr, h2w1, h2w2)

    h1_overlap = set(h1_win[(h1_win['End'] >= min(h1i1, h1i2)) & (h1_win['Start'] <= max(h1i1, h1i2))]['GeneID'])
    h2_overlap = set(h2_win[(h2_win['End'] >= min(h2i1, h2i2)) & (h2_win['Start'] <= max(h2i1, h2i2))]['GeneID'])

    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Prepare anchor sets
    anchor_roots_h1 = set()
    anchor_roots_h2 = set()
    for _, r in pairs_df.iterrows():
        anchor_roots_h1.add(re.sub(r'\.t\d+$', '', str(r['H1_gene'])))
        anchor_roots_h2.add(re.sub(r'\.t\d+$', '', str(r['H2_gene'])))

    # Zoom windows (paper-style: detailed view centered on inversion with smaller flank)
    zoom_flank_h1 = int(min(200_000, max(50_000, 0.10 * abs(h1i2 - h1i1 + 1))))
    zoom_flank_h2 = int(min(200_000, max(50_000, 0.10 * abs(h2i2 - h2i1 + 1))))
    h1z1 = max(1, min(h1i1, h1i2) - zoom_flank_h1)
    h1z2 = max(h1i1, h1i2) + zoom_flank_h1
    h2z1 = max(1, min(h2i1, h2i2) - zoom_flank_h2)
    h2z2 = max(h2i1, h2i2) + zoom_flank_h2

    h1_zoom = subset_window(h1_tx_t1, ref_chr, h1z1, h1z2)
    h2_zoom = subset_window(h2_tx_t1, qry_chr, h2z1, h2z2)
    h1_zoom_pos = {r['GeneID']: r for _, r in h1_zoom.iterrows()}
    h2_zoom_pos = {r['GeneID']: r for _, r in h2_zoom.iterrows()}

    # Layout blocks: overview + zoom (common in high-quality SV/genome figures)
    y_ov_h1 = (0.80, 0.90)
    y_ov_h2 = (0.66, 0.76)
    y_zoom_h1 = (0.42, 0.58)
    y_rib = (0.26, 0.40)
    y_zoom_h2 = (0.08, 0.24)

    draw_overview_track(ax, h1_win, h1w1, h1w2, h1i1, h1i2, 'H1', '#4c78a8', '#4c78a8', y_ov_h1[0], y_ov_h1[1])
    draw_overview_track(ax, h2_win, h2w1, h2w2, h2i1, h2i2, 'H2', '#f28e2b', '#f28e2b', y_ov_h2[0], y_ov_h2[1])
    _draw_axis_ticks_text(ax, h1w1, h1w2, y_ov_h1[1] + 0.008, where='above', n=5, fs=5.8)
    _draw_axis_ticks_text(ax, h2w1, h2w2, y_ov_h2[0] - 0.008, where='below', n=5, fs=5.8)

    # Zoom labels: all inversion-overlap genes if manageable, otherwise anchor-only + some longest
    h1_labels = set(h1_overlap)
    h2_labels = set(h2_overlap)
    if len(h1_labels) > 16:
        keep = set(list(sorted(anchor_roots_h1 & h1_overlap))[:16])
        if len(keep) < 16:
            lengths = []
            for _, r in h1_zoom.iterrows():
                if r['GeneID'] in h1_overlap and r['GeneID'] not in keep:
                    lengths.append((r['End'] - r['Start'], r['GeneID']))
            lengths.sort(reverse=True)
            keep.update([g for _, g in lengths[:16-len(keep)]])
        h1_labels = keep
    if len(h2_labels) > 16:
        keep = set(list(sorted(anchor_roots_h2 & h2_overlap))[:16])
        if len(keep) < 16:
            lengths = []
            for _, r in h2_zoom.iterrows():
                if r['GeneID'] in h2_overlap and r['GeneID'] not in keep:
                    lengths.append((r['End'] - r['Start'], r['GeneID']))
            lengths.sort(reverse=True)
            keep.update([g for _, g in lengths[:16-len(keep)]])
        h2_labels = keep

    draw_zoom_gene_track(
        ax, h1_zoom, h1z1, h1z2, h1i1, h1i2,
        h1_overlap, anchor_roots_h1, 'H1', '#2f6fab', '#174f8a', '#cfd9e8',
        y0=y_zoom_h1[0], y1=y_zoom_h1[1], label_roots=h1_labels, label_above=True
    )
    draw_zoom_gene_track(
        ax, h2_zoom, h2z1, h2z2, h2i1, h2i2,
        h2_overlap, anchor_roots_h2, 'H2', '#e58a2c', '#bd6316', '#edd9bf',
        y0=y_zoom_h2[0], y1=y_zoom_h2[1], label_roots=h2_labels, label_above=False
    )
    _draw_axis_ticks_text(ax, h1z1, h1z2, y_zoom_h1[1] + 0.030, where='above', n=5, fs=6.0)
    _draw_axis_ticks_text(ax, h2z1, h2z2, y_zoom_h2[0] - 0.030, where='below', n=5, fs=6.0)

    n_ribbons = draw_anchor_ribbons(ax, pairs_df, h1_zoom_pos, h2_zoom_pos, h1z1, h1z2, h2z1, h2z2, y_top=y_rib[1], y_bottom=y_rib[0])

    # Titles and legends
    ax.text(0.01, 0.985, f"{inv_id} gene organization across inversion (overview + zoom, t1 only)",
            ha='left', va='top', fontsize=8.6, fontweight='bold')
    ax.text(0.01, 0.958,
            f"Overview windows: H1 {ref_chr}:{h1w1:,}-{h1w2:,} | H2 {qry_chr}:{h2w1:,}-{h2w2:,}",
            ha='left', va='top', fontsize=6.9, color='#333')
    ax.text(0.01, 0.942,
            f"Zoom windows: H1 {ref_chr}:{h1z1:,}-{h1z2:,} | H2 {qry_chr}:{h2z1:,}-{h2z2:,}",
            ha='left', va='top', fontsize=6.9, color='#333')
    ax.text(0.99, 0.958, f"local flank={short_num_bp(gene_windows['flank'])}; zoom flanks H1/H2={short_num_bp(zoom_flank_h1)}/{short_num_bp(zoom_flank_h2)}",
            ha='right', va='top', fontsize=6.8, color='#333')

    # Compact legend
    lx, ly = 0.01, 0.612
    ax.add_patch(patches.Rectangle((lx, ly), 0.016, 0.013, facecolor='#fde8e7', edgecolor='none', transform=ax.transAxes))
    ax.text(lx+0.020, ly+0.0065, 'Inversion span', va='center', fontsize=6.2, transform=ax.transAxes)
    ax.add_patch(patches.Rectangle((lx+0.17, ly), 0.016, 0.013, facecolor='#2f6fab', edgecolor='none', transform=ax.transAxes))
    ax.text(lx+0.190, ly+0.0065, 'overlap genes', va='center', fontsize=6.2, transform=ax.transAxes)
    ax.add_patch(patches.Rectangle((lx+0.33, ly), 0.016, 0.013, facecolor='#174f8a', edgecolor='none', transform=ax.transAxes))
    ax.text(lx+0.350, ly+0.0065, 'anchor-paired genes', va='center', fontsize=6.2, transform=ax.transAxes)
    ax.plot([lx+0.55, lx+0.58], [ly+0.0065, ly+0.0065], color='#e76f51', lw=1.2, transform=ax.transAxes)
    ax.text(lx+0.585, ly+0.0065, 'anchor pair opposite strand', va='center', fontsize=6.2, transform=ax.transAxes)
    ax.plot([lx+0.86, lx+0.89], [ly+0.0065, ly+0.0065], color='#2a9d8f', lw=1.2, transform=ax.transAxes)
    ax.text(lx+0.895, ly+0.0065, 'same strand', va='center', fontsize=6.2, transform=ax.transAxes)
    ax.text(0.01, 0.595, f"Zoom labels show inversion-overlap genes (all if <=16; otherwise anchor-priority subset). Ribbons shown when both anchor genes fall in zoom windows (n={n_ribbons}).",
            ha='left', va='top', fontsize=5.8, color='#555')

    stats = {
        'H1_window_tx_all': int(len(h1_win_all)),
        'H1_window_tx_t1': int(len(h1_win)),
        'H2_window_tx_all': int(len(h2_win_all)),
        'H2_window_tx_t1': int(len(h2_win)),
        'H1_overlap_t1': int(len(h1_overlap)),
        'H2_overlap_t1': int(len(h2_overlap)),
        'anchor_pairs_plotted': int(n_ribbons),
    }
    return stats, h1_overlap, h2_overlap


def draw_expression_panel(ax_bar, ax_heat, ax_func, expr_df, sample_cols, func_map_h1, anchor_pairs_df):
    ax_bar.set_facecolor('white')
    ax_heat.set_facecolor('white')
    ax_func.set_facecolor('white')

    if expr_df.empty:
        for a in (ax_bar, ax_heat, ax_func):
            a.axis('off')
            a.text(0.5, 0.5, 'No RNA expression rows', ha='center', va='center', fontsize=8)
        return

    df = expr_df.copy().sort_values(['GeneStart', 'GeneEnd', 'GeneID']).reset_index(drop=True)
    gene_roots = df['GeneID'].tolist()
    mat = df[sample_cols].to_numpy(dtype=float)
    logm = np.log2(mat + 1.0)
    row_mean = logm.mean(axis=1, keepdims=True)
    row_std = logm.std(axis=1, ddof=0, keepdims=True)
    row_std[row_std == 0] = 1.0
    z = (logm - row_mean) / row_std
    z = np.clip(z, -2.5, 2.5)

    # anchor flags and function strings
    anchor_roots = set(re.sub(r'\.t\d+$', '', g) for g in anchor_pairs_df['H1_gene'].tolist()) if not anchor_pairs_df.empty else set()
    func_labels = []
    for g in gene_roots:
        f = func_map_h1.get(g, {'label': 'Unannotated'})
        func_labels.append(truncate_text(f.get('label', 'Unannotated'), 58))

    # Mean expression bar uses log2(count+1) to compress dynamic range and keep zeros.
    mean_vals = np.log2(df['MeanExpr'].to_numpy(dtype=float) + 1.0)
    y = np.arange(len(df))
    bar_colors = ['#1f5b99' if g in anchor_roots else '#8fb1d8' for g in gene_roots]
    ax_bar.barh(y, mean_vals, color=bar_colors, edgecolor='none', height=0.82)
    ax_bar.set_ylim(len(df)-0.5, -0.5)
    ax_bar.set_xlabel('Mean RNA-seq count (log2 scale; +1 pseudocount)', fontsize=6.6)
    ax_bar.tick_params(axis='x', labelsize=6)
    ax_bar.tick_params(axis='y', left=False, labelleft=False)
    for sp in ['top', 'right', 'left']:
        ax_bar.spines[sp].set_visible(False)
    ax_bar.grid(axis='x', color='#dddddd', lw=0.5, alpha=0.8)

    im = ax_heat.imshow(z, aspect='auto', interpolation='nearest', cmap='RdBu_r', vmin=-2.5, vmax=2.5)
    ax_heat.set_yticks(np.arange(len(df)))
    ylabels = [f"{g}.t1" if g not in anchor_roots else f"{g}.t1*" for g in gene_roots]
    ax_heat.set_yticklabels(ylabels, fontsize=6.2)
    ax_heat.set_xticks(np.arange(len(sample_cols)))
    ax_heat.set_xticklabels(sample_cols, rotation=90, fontsize=5.6)
    ax_heat.tick_params(axis='both', length=0)
    for sp in ax_heat.spines.values():
        sp.set_visible(False)
    ax_heat.set_xlabel('RNA-seq samples (heatmap = row z-score after log2(count+1))', fontsize=6.6)

    ax_func.set_xlim(0, 1)
    ax_func.set_ylim(len(df)-0.5, -0.5)
    ax_func.axis('off')
    ax_func.text(0.0, -0.9, 'Function (H1 t1; MapMan/Mercator representative label)', fontsize=7, fontweight='bold')
    for i, (g, label) in enumerate(zip(gene_roots, func_labels)):
        is_anchor = g in anchor_roots
        dot_color = '#1f5b99' if is_anchor else '#aaaaaa'
        ax_func.add_patch(patches.Circle((0.015, i), radius=0.12, transform=ax_func.transData,
                                         facecolor=dot_color, edgecolor='none'))
        ax_func.text(0.04, i, label, va='center', fontsize=6.0, color='#222')

    return im


def draw_hic_panel(ax, img_path, hic_row, inv_id):
    ax.set_facecolor('white')
    if not os.path.exists(img_path):
        ax.axis('off')
        ax.text(0.5, 0.5, f'Missing Hi-C panel\n{inv_id}', ha='center', va='center', fontsize=9)
        return
    arr = trim_white_margin(img_path)
    ax.imshow(arr)
    ax.axis('off')
    if hic_row is not None:
        title = f"Hi-C inter-haplotype context ({hic_row.get('Norm','NA')}, {int(hic_row.get('bin_size',0))/1000:.0f} kb bins)"
        subtitle = f"Contacts={int(hic_row.get('Inter_contacts', 0)):,}; flank={short_num_bp(hic_row.get('flank_bp', 0))}"
    else:
        title = 'Hi-C inter-haplotype context'
        subtitle = 'No summary row found'
    ax.text(0.0, 1.01, title, ha='left', va='bottom', transform=ax.transAxes, fontsize=8.5, fontweight='bold')
    ax.text(0.0, 0.985, subtitle, ha='left', va='top', transform=ax.transAxes, fontsize=7.0, color='#333')


def main():
    args = parse_args()
    inv_ids = [x.strip() for x in args.inv_ids.split(',') if x.strip()]
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    inv_df = load_inversions(args.inv_tsv)
    inv_df = inv_df[inv_df['ID'].isin(inv_ids)].copy()
    if inv_df.empty:
        raise SystemExit(f"No requested inversions found in {args.inv_tsv}: {inv_ids}")
    inv_df['ID'] = pd.Categorical(inv_df['ID'], categories=inv_ids, ordered=True)
    inv_df = inv_df.sort_values('ID').reset_index(drop=True)

    expr_by_inv, sample_cols = load_expression(args.rna_genes, inv_ids)
    anchor_df = load_anchor_pairs(args.anchor_pairs, inv_ids)
    gene_counts = load_gene_counts(args.gene_counts)
    hic_summary = load_hic_summary(args.hic_summary)
    func_h1 = load_function_map(args.h1_func)
    func_h2 = load_function_map(args.h2_func)  # loaded for completeness; future use

    target_chroms_h1 = set(inv_df['RefChr'])
    target_chroms_h2 = set(inv_df['QryChr'])
    h1_tx_all = parse_gtf_transcripts(args.h1_gtf, target_chroms_h1)
    h2_tx_all = parse_gtf_transcripts(args.h2_gtf, target_chroms_h2)
    h1_tx_t1 = h1_tx_all[h1_tx_all['Is_t1']].copy()
    h2_tx_t1 = h2_tx_all[h2_tx_all['Is_t1']].copy()

    # Figure layout: two rows (one per inversion), 4 main panels per row
    nrows = len(inv_df)
    fig = plt.figure(figsize=(22, 8.4 * nrows), dpi=300, facecolor='white')
    outer = fig.add_gridspec(nrows=nrows, ncols=4, width_ratios=[8.8, 1.7, 6.8, 5.8], hspace=0.30, wspace=0.18)

    # header
    fig.text(0.012, 0.995, 'Comparative inversion panels for INV1936 and INV1915 (H1 vs H2)',
             ha='left', va='top', fontsize=15, fontweight='bold')
    fig.text(0.012, 0.980,
             'Left panel uses overview+zoom gene organization (t1-only). Middle panels summarize H1 inversion-overlap RNA counts (bar = mean log2(count+1); heatmap = within-gene z-score after log2 transform). Right panel shows inter-haplotype Hi-C context.',
             ha='left', va='top', fontsize=8.5, color='#333')

    summary_rows = []
    last_heat_im = None

    for i, (_, inv) in enumerate(inv_df.iterrows()):
        inv_id = str(inv['ID'])
        gene_windows = build_gene_window(inv, args.gene_flank_min, args.gene_flank_max)
        pairs_sub = anchor_df[anchor_df['INV_ID'] == inv_id].copy()

        ax_syn = fig.add_subplot(outer[i, 0])
        stats, h1_overlap_t1, h2_overlap_t1 = draw_synteny_panel(
            ax_syn, inv_id, inv, gene_windows, h1_tx_all, h2_tx_all, h1_tx_t1, h2_tx_t1, pairs_sub
        )

        ax_bar = fig.add_subplot(outer[i, 1])
        mid = outer[i, 2].subgridspec(1, 2, width_ratios=[4.6, 5.4], wspace=0.03)
        ax_heat = fig.add_subplot(mid[0, 0])
        ax_func = fig.add_subplot(mid[0, 1], sharey=ax_heat)
        expr_sub = expr_by_inv.get(inv_id, pd.DataFrame())
        last_heat_im = draw_expression_panel(ax_bar, ax_heat, ax_func, expr_sub, sample_cols, func_h1, pairs_sub)

        hic_ax = fig.add_subplot(outer[i, 3])
        hic_png = os.path.join(args.hic_png_dir, f"{inv_id}.hic_context.juicer.png")
        draw_hic_panel(hic_ax, hic_png, hic_summary.get(inv_id), inv_id)

        # row summary header above bar/heat panels
        row_top = ax_syn.get_position().y1
        row_left = ax_syn.get_position().x0
        row_right = hic_ax.get_position().x1
        ref_len = abs(int(inv['RefEnd']) - int(inv['RefStart'])) + 1
        qry_len = abs(int(inv['QryEnd']) - int(inv['QryStart'])) + 1
        gc_row = gene_counts.get(inv_id)
        extra = ''
        if gc_row is not None:
            extra = (f"; overlap genes={int(gc_row['overlap_genes_total'])}, anchor-paired={int(gc_row['anchor_paired_genes_used'])}, "
                     f"missing anchor mates={int(gc_row['missing_from_anchor_pair_set'])}")
        fig.text(row_left, row_top + 0.010,
                 f"{inv_id}: H1 {inv['RefChr']} {int(inv['RefStart']):,}-{int(inv['RefEnd']):,} ({short_num_bp(ref_len)})  |  "
                 f"H2 {inv['QryChr']} {int(inv['QryStart']):,}-{int(inv['QryEnd']):,} ({short_num_bp(qry_len)})"
                 f"{extra}",
                 ha='left', va='bottom', fontsize=9.5, fontweight='bold')

        # add a small note on isoform filtering in row
        fig.text(row_right, row_top + 0.010,
                 f"t1 filter in local windows: H1 {stats['H1_window_tx_t1']}/{stats['H1_window_tx_all']} tx, H2 {stats['H2_window_tx_t1']}/{stats['H2_window_tx_all']} tx",
                 ha='right', va='bottom', fontsize=7.5, color='#444')

        summary_rows.append({
            'INV_ID': inv_id,
            'RefChr': inv['RefChr'], 'RefStart': int(inv['RefStart']), 'RefEnd': int(inv['RefEnd']),
            'QryChr': inv['QryChr'], 'QryStart': int(inv['QryStart']), 'QryEnd': int(inv['QryEnd']),
            'RefLen_bp': ref_len, 'QryLen_bp': qry_len,
            'GeneWindowFlank_bp': gene_windows['flank'],
            'H1_window_transcripts_all': stats['H1_window_tx_all'],
            'H1_window_transcripts_t1': stats['H1_window_tx_t1'],
            'H2_window_transcripts_all': stats['H2_window_tx_all'],
            'H2_window_transcripts_t1': stats['H2_window_tx_t1'],
            'H1_overlap_genes_t1_from_gtf': stats['H1_overlap_t1'],
            'H2_overlap_genes_t1_from_gtf': stats['H2_overlap_t1'],
            'H1_overlap_genes_in_rna_table': int(len(expr_sub)),
            'AnchorPairs_plotted_t1': stats['anchor_pairs_plotted'],
            'HiC_panel_png': hic_png,
            'HiC_inter_contacts': int(hic_summary.get(inv_id, {}).get('Inter_contacts', 0)) if inv_id in hic_summary else np.nan,
        })

    # shared heatmap colorbar
    if last_heat_im is not None:
        cax = fig.add_axes([0.915, 0.965 - 0.015*nrows, 0.07, 0.012])
        cb = fig.colorbar(last_heat_im, cax=cax, orientation='horizontal')
        cb.set_label('Heatmap color: row z-score after log2(count+1)', fontsize=7)
        cb.ax.tick_params(labelsize=6, length=2)

    # footer notes
    fig.text(0.012, 0.012,
             '* Anchor-paired genes are from selected_inversions.anchor_divergence_pairs.tsv. H1/H2 gene tracks were explicitly filtered to .t1 transcripts at plot time. Middle bar plot shows mean raw count on log2 scale with +1 pseudocount; heatmap colors show relative sample-to-sample changes within each gene (row z-score), not absolute expression magnitude.',
             ha='left', va='bottom', fontsize=7.3, color='#333')

    pdf_path = outdir / f"{args.prefix}.pdf"
    png_path = outdir / f"{args.prefix}.png"
    fig.savefig(pdf_path, bbox_inches='tight', facecolor='white')
    fig.savefig(png_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)

    summary_tsv = outdir / f"{args.prefix}.summary.tsv"
    pd.DataFrame(summary_rows).to_csv(summary_tsv, sep='\t', index=False)

    # Export per-inversion function table for traceability
    func_export_rows = []
    for inv_id in inv_ids:
        expr_sub = expr_by_inv.get(inv_id, pd.DataFrame())
        if expr_sub.empty:
            continue
        pairs_sub = anchor_df[anchor_df['INV_ID'] == inv_id].copy()
        anchor_roots = set(re.sub(r'\.t\d+$', '', g) for g in pairs_sub['H1_gene'].tolist())
        for _, r in expr_sub.sort_values(['GeneStart', 'GeneEnd', 'GeneID']).iterrows():
            g = str(r['GeneID'])
            f = func_h1.get(g, {'label': 'Unannotated'})
            func_export_rows.append({
                'INV_ID': inv_id,
                'GeneID': g,
                'GeneID_t1': g + '.t1',
                'GeneChr': r['GeneChr'],
                'GeneStart': int(r['GeneStart']),
                'GeneEnd': int(r['GeneEnd']),
                'MeanExpr': float(r['MeanExpr']),
                'Log2MeanExpr': float(r['Log2MeanExpr']),
                'AnchorPaired': g in anchor_roots,
                'FunctionLabel_H1_t1': f.get('label', 'Unannotated'),
            })
    func_tsv = outdir / f"{args.prefix}.gene_function_expression_table.tsv"
    pd.DataFrame(func_export_rows).to_csv(func_tsv, sep='\t', index=False)

    print(f"[OK] Figure PDF: {pdf_path}")
    print(f"[OK] Figure PNG: {png_path}")
    print(f"[OK] Summary TSV: {summary_tsv}")
    print(f"[OK] Gene function table: {func_tsv}")


if __name__ == '__main__':
    main()
