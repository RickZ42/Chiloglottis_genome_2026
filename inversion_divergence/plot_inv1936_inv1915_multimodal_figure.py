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


def draw_gene_track_block(ax, transcripts, window_start, window_end, inv_start, inv_end,
                          highlight_roots, anchor_roots, hap_label, color_main, color_anchor,
                          color_bg, y0, y1, tick_side='top'):
    # background band and inversion shading
    ax.add_patch(patches.Rectangle((0, y0), 1, y1-y0, facecolor='#f7f7f7', edgecolor='none', zorder=0))
    x_inv1 = norm_pos(inv_start, window_start, window_end)
    x_inv2 = norm_pos(inv_end, window_start, window_end)
    xlo, xhi = sorted([x_inv1, x_inv2])
    ax.add_patch(patches.Rectangle((xlo, y0), max(1e-6, xhi-xlo), y1-y0, facecolor='#ffe4e1', edgecolor='none', alpha=0.9, zorder=0.2))
    for xv in [xlo, xhi]:
        ax.plot([xv, xv], [y0, y1], color='#c33', lw=0.9, ls='--', zorder=2)

    ax.text(-0.015, (y0+y1)/2, hap_label, ha='right', va='center', fontsize=8, fontweight='bold', transform=ax.transAxes)
    ax.plot([0, 1], [y0+0.02, y0+0.02], color='#9a9a9a', lw=0.5, zorder=1)

    if transcripts.empty:
        ax.text(0.5, (y0+y1)/2, 'No t1 transcripts in window', ha='center', va='center', fontsize=7, color='#777')
        return

    tx = transcripts.sort_values(['Start', 'End', 'TranscriptID']).copy()
    tx['lane'] = pack_lanes(tx['Start'].tolist(), tx['End'].tolist(), pad_bp=max(1000, int((window_end-window_start)*0.003)))
    n_lanes = int(tx['lane'].max()) + 1 if len(tx) else 1
    lane_pad = 0.01
    band_h = (y1 - y0) - 0.06
    lane_h = max(0.025, band_h / max(1, n_lanes))
    draw_h = min(0.05, lane_h * 0.7)

    for _, r in tx.iterrows():
        x1 = norm_pos(r['Start'], window_start, window_end)
        x2 = norm_pos(r['End'], window_start, window_end)
        x1, x2 = sorted([x1, x2])
        lane = int(r['lane'])
        y = y0 + 0.04 + lane * lane_h + lane_pad
        gene = r['GeneID']
        is_highlight = gene in highlight_roots
        is_anchor = gene in anchor_roots
        fc = color_bg
        ec = '#999999'
        alpha = 0.75
        lw = 0.4
        if is_highlight:
            fc = color_main
            ec = color_main
            alpha = 0.88
            lw = 0.6
        if is_anchor:
            fc = color_anchor
            ec = color_anchor
            alpha = 0.95
            lw = 0.7
        width = max(0.0015, x2 - x1)
        if width > 0.01:
            head = min(width * 0.25, 0.012)
            if r['Strand'] == '+':
                arr_x = x1
                arr_dx = width
            else:
                arr_x = x2
                arr_dx = -width
            ax.add_patch(patches.FancyArrow(
                arr_x, y + draw_h/2, arr_dx, 0,
                width=draw_h*0.65, head_width=draw_h*0.95, head_length=head,
                length_includes_head=True, facecolor=fc, edgecolor=ec, linewidth=lw, alpha=alpha, zorder=3
            ))
        else:
            ax.add_patch(patches.Rectangle((x1, y), width, draw_h, facecolor=fc, edgecolor=ec, lw=lw, alpha=alpha, zorder=3))

    # Mb labels for this haplotype window
    ticks = np.linspace(0, 1, 5)
    tick_labels = [f"{(window_start + t*(window_end-window_start))/1e6:.2f}" for t in ticks]
    ytxt = y1 + 0.005 if tick_side == 'top' else y0 - 0.018
    va = 'bottom' if tick_side == 'top' else 'top'
    for t, lab in zip(ticks, tick_labels):
        ax.text(t, ytxt, lab, ha='center', va=va, fontsize=6.5, color='#444')
    unit_y = y1 + 0.028 if tick_side == 'top' else y0 - 0.04
    ax.text(0.995, unit_y, 'Mb', ha='right', va='bottom' if tick_side=='top' else 'top', fontsize=6.5, color='#444')


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

    h1_pos = {r['GeneID']: r for _, r in h1_win.iterrows()}
    h2_pos = {r['GeneID']: r for _, r in h2_win.iterrows()}
    anchor_roots_h1 = set()
    anchor_roots_h2 = set()

    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Draw ribbons first (under genes)
    y_h1 = (0.57, 0.94)
    y_h2 = (0.06, 0.43)
    for _, r in pairs_df.iterrows():
        h1_root = re.sub(r'\.t\d+$', '', str(r['H1_gene']))
        h2_root = re.sub(r'\.t\d+$', '', str(r['H2_gene']))
        if h1_root not in h1_pos or h2_root not in h2_pos:
            continue
        anchor_roots_h1.add(h1_root)
        anchor_roots_h2.add(h2_root)
        h1r = h1_pos[h1_root]
        h2r = h2_pos[h2_root]
        x1 = norm_pos((h1r['Start'] + h1r['End'])/2, h1w1, h1w2)
        x2 = norm_pos((h2r['Start'] + h2r['End'])/2, h2w1, h2w2)
        y1 = y_h1[0] + 0.02
        y2 = y_h2[1] - 0.02
        same_strand = (h1r['Strand'] == h2r['Strand'])
        color = '#2a9d8f' if same_strand else '#e76f51'
        alpha = 0.65
        verts = [
            (x1, y1),
            (x1, 0.52),
            (x2, 0.48),
            (x2, y2),
        ]
        codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4]
        patch = PathPatch(MplPath(verts, codes), facecolor='none', edgecolor=color, lw=1.2, alpha=alpha, zorder=1)
        ax.add_patch(patch)

    draw_gene_track_block(
        ax, h1_win, h1w1, h1w2, h1i1, h1i2,
        h1_overlap, anchor_roots_h1,
        'H1', color_main='#4c78a8', color_anchor='#1f5b99', color_bg='#d0d7df',
        y0=y_h1[0], y1=y_h1[1], tick_side='top'
    )
    draw_gene_track_block(
        ax, h2_win, h2w1, h2w2, h2i1, h2i2,
        h2_overlap, anchor_roots_h2,
        'H2', color_main='#f28e2b', color_anchor='#c76a12', color_bg='#ead8c2',
        y0=y_h2[0], y1=y_h2[1], tick_side='bottom'
    )

    ax.text(0.01, 0.985, f"{inv_id} synteny + t1 gene layout", ha='left', va='top', fontsize=8.5, fontweight='bold')
    ax.text(0.01, 0.955, f"{ref_chr}:{h1i1:,}-{h1i2:,} (H1) vs {qry_chr}:{h2i1:,}-{h2i2:,} (H2)",
            ha='left', va='top', fontsize=7.2, color='#333')
    ax.text(0.99, 0.955, f"Local gene window flank = {short_num_bp(gene_windows['flank'])}",
            ha='right', va='top', fontsize=7.0, color='#333')

    # legend
    lx, ly = 0.01, 0.50
    ax.add_patch(patches.Rectangle((lx, ly), 0.018, 0.018, facecolor='#ffe4e1', edgecolor='none', transform=ax.transAxes))
    ax.text(lx+0.024, ly+0.009, 'Inversion interval', va='center', fontsize=6.5, transform=ax.transAxes)
    ax.plot([lx+0.24, lx+0.27], [ly+0.009, ly+0.009], color='#e76f51', lw=1.2, transform=ax.transAxes)
    ax.text(lx+0.275, ly+0.009, 'Anchor pair (opposite strand)', va='center', fontsize=6.5, transform=ax.transAxes)
    ax.plot([lx+0.66, lx+0.69], [ly+0.009, ly+0.009], color='#2a9d8f', lw=1.2, transform=ax.transAxes)
    ax.text(lx+0.695, ly+0.009, 'Anchor pair (same strand)', va='center', fontsize=6.5, transform=ax.transAxes)

    stats = {
        'H1_window_tx_all': int(len(h1_win_all)),
        'H1_window_tx_t1': int(len(h1_win)),
        'H2_window_tx_all': int(len(h2_win_all)),
        'H2_window_tx_t1': int(len(h2_win)),
        'H1_overlap_t1': int(len(h1_overlap)),
        'H2_overlap_t1': int(len(h2_overlap)),
        'anchor_pairs_plotted': int(len(anchor_roots_h1)),
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

    # mean expression bar (raw mean counts, log2-scaled axis label)
    mean_vals = np.log2(df['MeanExpr'].to_numpy(dtype=float) + 1.0)
    y = np.arange(len(df))
    bar_colors = ['#1f5b99' if g in anchor_roots else '#8fb1d8' for g in gene_roots]
    ax_bar.barh(y, mean_vals, color=bar_colors, edgecolor='none', height=0.82)
    ax_bar.set_ylim(len(df)-0.5, -0.5)
    ax_bar.set_xlabel('log2(mean count+1)', fontsize=7)
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
    ax_heat.set_xlabel('RNA-seq samples (raw counts; row z-score heatmap)', fontsize=7)

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
    fig = plt.figure(figsize=(20, 8.2 * nrows), dpi=300)
    outer = fig.add_gridspec(nrows=nrows, ncols=4, width_ratios=[7.0, 1.7, 7.0, 6.2], hspace=0.30, wspace=0.18)

    # header
    fig.text(0.012, 0.995, 'Comparative inversion panels for INV1936 and INV1915 (H1 vs H2)',
             ha='left', va='top', fontsize=15, fontweight='bold')
    fig.text(0.012, 0.980,
             'Tracks show local synteny and t1-only transcript distribution; RNA heatmaps use H1 inversion-overlapping genes (row z-score on raw counts); Hi-C panels show inter-haplotype context.',
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
        cb.set_label('Row z-score of log2(count+1)', fontsize=7)
        cb.ax.tick_params(labelsize=6, length=2)

    # footer notes
    fig.text(0.012, 0.012,
             '* Anchor-paired genes (from selected_inversions.anchor_divergence_pairs.tsv). Functional labels are representative H1 t1 Mercator/MapMan annotations; H1/H2 gene-track isoforms were explicitly filtered to .t1 at plot time.',
             ha='left', va='bottom', fontsize=7.3, color='#333')

    pdf_path = outdir / f"{args.prefix}.pdf"
    png_path = outdir / f"{args.prefix}.png"
    fig.savefig(pdf_path, bbox_inches='tight')
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
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
