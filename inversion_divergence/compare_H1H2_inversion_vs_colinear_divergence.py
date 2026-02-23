#!/usr/bin/env python3
import math
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
from scipy.stats import mannwhitneyu

# Paths for current project
BASE = Path('/g/data/xf3/zz3507')
ANCHORS = BASE / 'compare_H1_vs_Ophrys/New_H1_H2_ophrys_Arobx/H1t1.H2t1.anchors'
H1_CDS = BASE / 'compare_H1_vs_Ophrys/New_H1_H2_ophrys_Arobx/H1t1.cds'
H2_CDS = BASE / 'compare_H1_vs_Ophrys/New_H1_H2_ophrys_Arobx/H2t1.cds'
INV_H1_GENES = BASE / 'Output/20260127Genome/H1/breaker/H1_INV_overlapped_genes.unique.txt'
OUTDIR = BASE / 'Output/20260127Genome/H1H2_inversion_divergence'

VALID_NT = set('ACGTN')
STOP_CODONS = {'TAA', 'TAG', 'TGA'}
AA_GAP = '-'
AA_ALIGN_PARAMS = dict(match=2, mismatch=-1, open=-10, extend=-0.5)


def load_fasta(path: Path):
    seqs = {}
    cur = None
    parts = []
    with path.open() as f:
        for line in f:
            if line.startswith('>'):
                if cur is not None:
                    seqs[cur] = ''.join(parts)
                cur = line[1:].strip().split()[0]
                parts = []
            else:
                parts.append(line.strip())
        if cur is not None:
            seqs[cur] = ''.join(parts)
    return seqs


def parse_anchors(path: Path):
    pairs = []
    with path.open() as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.split()
            if len(fields) < 3:
                continue
            h1, h2, score = fields[:3]
            try:
                score = float(score)
            except ValueError:
                score = math.nan
            pairs.append((h1, h2, score))
    return pairs


def clean_cds(seq: str):
    s = re.sub(r'\s+', '', seq).upper()
    illegal = sum(1 for c in s if c not in VALID_NT)
    if illegal:
        s = ''.join(c if c in VALID_NT else 'N' for c in s)
    trailing_stop_removed = 0
    while len(s) >= 3 and len(s) % 3 == 0 and s[-3:] in STOP_CODONS:
        s = s[:-3]
        trailing_stop_removed += 1
    return s, illegal, trailing_stop_removed


def translate_cds(cds: str):
    # Biopython translates N-containing codons to X
    return str(Seq(cds).translate(table=1, to_stop=False))


def aa_global_align(a: str, b: str):
    if len(a) == len(b):
        return a, b, False
    aln = pairwise2.align.globalms(
        a,
        b,
        AA_ALIGN_PARAMS['match'],
        AA_ALIGN_PARAMS['mismatch'],
        AA_ALIGN_PARAMS['open'],
        AA_ALIGN_PARAMS['extend'],
        one_alignment_only=True,
    )
    if not aln:
        raise RuntimeError('protein alignment failed')
    return aln[0].seqA, aln[0].seqB, True


def backtranslate_codon_alignment(cds1: str, cds2: str, aa1_aln: str, aa2_aln: str):
    i1 = 0
    i2 = 0
    cod1_parts = []
    cod2_parts = []
    for x, y in zip(aa1_aln, aa2_aln):
        if x == AA_GAP:
            cod1_parts.append('---')
        else:
            cod1_parts.append(cds1[i1:i1 + 3])
            i1 += 3
        if y == AA_GAP:
            cod2_parts.append('---')
        else:
            cod2_parts.append(cds2[i2:i2 + 3])
            i2 += 3
    if i1 != len(cds1) or i2 != len(cds2):
        raise RuntimeError(f'backtranslate index mismatch: {i1}/{len(cds1)} {i2}/{len(cds2)}')
    return ''.join(cod1_parts), ''.join(cod2_parts)


def pair_metrics(h1_id, h2_id, h1_cds_raw, h2_cds_raw, score, inv_roots):
    h1_root = h1_id.split('.')[0]
    h2_root = h2_id.split('.')[0]
    rec = {
        'H1_gene': h1_id,
        'H2_gene': h2_id,
        'H1_root': h1_root,
        'H2_root': h2_root,
        'anchor_score': score,
        'is_inversion_gene': h1_root in inv_roots,
    }

    if h1_cds_raw is None or h2_cds_raw is None:
        rec['status'] = 'missing_cds'
        return rec

    h1_cds, h1_illegal_nt, h1_trail_stop = clean_cds(h1_cds_raw)
    h2_cds, h2_illegal_nt, h2_trail_stop = clean_cds(h2_cds_raw)
    rec.update({
        'H1_cds_len_bp': len(h1_cds),
        'H2_cds_len_bp': len(h2_cds),
        'H1_illegal_nt_replaced': h1_illegal_nt,
        'H2_illegal_nt_replaced': h2_illegal_nt,
        'H1_trailing_stop_codons_removed': h1_trail_stop,
        'H2_trailing_stop_codons_removed': h2_trail_stop,
    })

    if len(h1_cds) == 0 or len(h2_cds) == 0:
        rec['status'] = 'empty_cds_after_cleaning'
        return rec
    if len(h1_cds) % 3 != 0 or len(h2_cds) % 3 != 0:
        rec['status'] = 'cds_not_mod3'
        return rec

    try:
        h1_aa = translate_cds(h1_cds)
        h2_aa = translate_cds(h2_cds)
    except Exception:
        rec['status'] = 'translation_error'
        return rec

    rec['H1_aa_len'] = len(h1_aa)
    rec['H2_aa_len'] = len(h2_aa)

    # Internal stop codons can break codon-based dN/dS assumptions; keep record and skip dN/dS.
    h1_internal_stop = '*' in h1_aa
    h2_internal_stop = '*' in h2_aa
    rec['H1_internal_stop'] = h1_internal_stop
    rec['H2_internal_stop'] = h2_internal_stop

    # Replace stops with X for protein alignment only (if any).
    h1_aa_aln_in = h1_aa.replace('*', 'X')
    h2_aa_aln_in = h2_aa.replace('*', 'X')

    try:
        aa1_aln, aa2_aln, used_alignment = aa_global_align(h1_aa_aln_in, h2_aa_aln_in)
    except Exception:
        rec['status'] = 'protein_alignment_error'
        return rec
    rec['used_pairwise_aa_alignment'] = used_alignment

    try:
        cod1_aln, cod2_aln = backtranslate_codon_alignment(h1_cds, h2_cds, aa1_aln, aa2_aln)
    except Exception:
        rec['status'] = 'backtranslate_error'
        return rec

    # Metrics excluding codons with gaps in either seq
    codon_pairs = []
    for i in range(0, len(cod1_aln), 3):
        c1 = cod1_aln[i:i+3]
        c2 = cod2_aln[i:i+3]
        if '-' in c1 or '-' in c2:
            continue
        codon_pairs.append((c1, c2))

    if not codon_pairs:
        rec['status'] = 'no_ungapped_codons'
        return rec

    codon_sites = len(codon_pairs)
    nt_sites = codon_sites * 3
    codon_matches = sum(1 for c1, c2 in codon_pairs if c1 == c2)
    nt_matches = sum(sum(1 for a, b in zip(c1, c2) if a == b) for c1, c2 in codon_pairs)

    aa_pairs = [(a, b) for a, b in zip(aa1_aln, aa2_aln) if a != '-' and b != '-']
    aa_sites = len(aa_pairs)
    aa_matches = sum(1 for a, b in aa_pairs if a == b)

    rec.update({
        'aligned_codon_sites': codon_sites,
        'aligned_nt_sites': nt_sites,
        'aligned_aa_sites': aa_sites,
        'codon_identity': codon_matches / codon_sites if codon_sites else np.nan,
        'cds_nt_identity': nt_matches / nt_sites if nt_sites else np.nan,
        'cds_nt_pdistance': 1 - (nt_matches / nt_sites) if nt_sites else np.nan,
        'aa_identity': aa_matches / aa_sites if aa_sites else np.nan,
        'aa_pdistance': 1 - (aa_matches / aa_sites) if aa_sites else np.nan,
        'codon_differences': codon_sites - codon_matches,
        'nt_differences': nt_sites - nt_matches,
        'aa_differences': aa_sites - aa_matches,
    })

    # dN/dS calculation (NG86)
    dN = np.nan
    dS = np.nan
    dNdS = np.nan
    ds_status = 'not_run'
    if h1_internal_stop or h2_internal_stop:
        ds_status = 'internal_stop_skip'
    else:
        try:
            _dN, _dS = cal_dn_ds(CodonSeq(cod1_aln), CodonSeq(cod2_aln), method='NG86')
            # Biopython can return sentinel negatives (e.g., -1) for undefined values.
            if _dN is not None and np.isfinite(_dN) and _dN >= 0:
                dN = float(_dN)
            if _dS is not None and np.isfinite(_dS) and _dS >= 0:
                dS = float(_dS)
            if np.isfinite(dN) and np.isfinite(dS):
                dNdS = np.inf if dS == 0 and dN > 0 else (dN / dS if dS > 0 else np.nan)
            ds_status = 'ok' if np.isfinite(dS) else 'undefined_ds'
        except Exception:
            ds_status = 'dn_ds_error'

    rec.update({'dN_NG86': dN, 'dS_NG86': dS, 'dN_dS_NG86': dNdS, 'dnds_status': ds_status})
    rec['status'] = 'ok'
    return rec


def safe_mwu(x, y):
    x = np.asarray([v for v in x if pd.notna(v)], dtype=float)
    y = np.asarray([v for v in y if pd.notna(v)], dtype=float)
    if len(x) == 0 or len(y) == 0:
        return np.nan, np.nan
    try:
        res = mannwhitneyu(x, y, alternative='two-sided')
        return float(res.statistic), float(res.pvalue)
    except Exception:
        return np.nan, np.nan


def fmt(v, digits=6):
    if v is None or (isinstance(v, float) and (math.isnan(v) or math.isinf(v))):
        return 'NA'
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    return f'{float(v):.{digits}f}'


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    inv_roots = {line.strip() for line in INV_H1_GENES.open() if line.strip()}
    h1_cds = load_fasta(H1_CDS)
    h2_cds = load_fasta(H2_CDS)
    anchors = parse_anchors(ANCHORS)

    h1_counts = Counter(h1 for h1, _, _ in anchors)
    h2_counts = Counter(h2 for _, h2, _ in anchors)
    unique_pairs = [(h1, h2, s) for h1, h2, s in anchors if h1_counts[h1] == 1 and h2_counts[h2] == 1]

    rows = []
    for i, (h1, h2, score) in enumerate(unique_pairs, 1):
        rows.append(pair_metrics(h1, h2, h1_cds.get(h1), h2_cds.get(h2), score, inv_roots))
        if i % 500 == 0:
            print(f'[INFO] processed {i}/{len(unique_pairs)} unique anchor pairs', file=sys.stderr)

    df = pd.DataFrame(rows)
    # enforce columns ordering (core first)
    core_cols = [
        'H1_gene', 'H2_gene', 'H1_root', 'H2_root', 'anchor_score', 'is_inversion_gene', 'status', 'dnds_status',
        'H1_cds_len_bp', 'H2_cds_len_bp', 'H1_aa_len', 'H2_aa_len',
        'aligned_codon_sites', 'aligned_nt_sites', 'aligned_aa_sites',
        'cds_nt_identity', 'cds_nt_pdistance', 'aa_identity', 'aa_pdistance', 'codon_identity',
        'nt_differences', 'aa_differences', 'codon_differences',
        'dN_NG86', 'dS_NG86', 'dN_dS_NG86',
        'H1_internal_stop', 'H2_internal_stop', 'used_pairwise_aa_alignment',
        'H1_illegal_nt_replaced', 'H2_illegal_nt_replaced',
        'H1_trailing_stop_codons_removed', 'H2_trailing_stop_codons_removed',
    ]
    cols = [c for c in core_cols if c in df.columns] + [c for c in df.columns if c not in core_cols]
    df = df[cols]

    out_tsv = OUTDIR / 'H1H2_anchor_unique_pairs.inversion_vs_colinear.divergence.tsv'
    df.to_csv(out_tsv, sep='\t', index=False)

    # Summaries
    ok = df[df['status'] == 'ok'].copy()
    inv = ok[ok['is_inversion_gene']].copy()
    col = ok[~ok['is_inversion_gene']].copy()
    ds_ok = ok[pd.notna(ok['dS_NG86'])].copy()
    inv_ds = ds_ok[ds_ok['is_inversion_gene']]
    col_ds = ds_ok[~ds_ok['is_inversion_gene']]

    inv_genes_anchored = set(inv['H1_root'])
    inv_roots_not_anchored = sorted(inv_roots - inv_genes_anchored)

    metrics = ['cds_nt_pdistance', 'cds_nt_identity', 'aa_pdistance', 'aa_identity', 'dN_NG86', 'dS_NG86']
    summary_rows = []
    for m in metrics:
        xi = inv[m].dropna().astype(float)
        xc = col[m].dropna().astype(float)
        stat, p = safe_mwu(xi, xc)
        summary_rows.append({
            'metric': m,
            'inv_n': len(xi),
            'inv_median': float(np.median(xi)) if len(xi) else np.nan,
            'inv_mean': float(np.mean(xi)) if len(xi) else np.nan,
            'col_n': len(xc),
            'col_median': float(np.median(xc)) if len(xc) else np.nan,
            'col_mean': float(np.mean(xc)) if len(xc) else np.nan,
            'mannwhitney_u': stat,
            'pvalue_two_sided': p,
        })
    summary_df = pd.DataFrame(summary_rows)
    summary_tsv = OUTDIR / 'H1H2_inversion_vs_colinear.metric_summary.tsv'
    summary_df.to_csv(summary_tsv, sep='\t', index=False)

    # Top divergent inversion genes (anchor-paired)
    inv_rank_nt = inv.sort_values(['cds_nt_pdistance', 'aa_pdistance'], ascending=False)
    inv_rank_ds = inv_ds.sort_values('dS_NG86', ascending=False)
    top_nt_tsv = OUTDIR / 'H1H2_inversion_pairs.top20_by_cds_nt_pdistance.tsv'
    top_ds_tsv = OUTDIR / 'H1H2_inversion_pairs.top20_by_dS_NG86.tsv'
    inv_rank_nt.head(20).to_csv(top_nt_tsv, sep='\t', index=False)
    inv_rank_ds.head(20).to_csv(top_ds_tsv, sep='\t', index=False)

    # Text summary for easy reporting
    txt = OUTDIR / 'H1H2_inversion_vs_colinear.divergence_summary.txt'
    with txt.open('w') as w:
        w.write('H1-H2 inversion vs colinear (anchor-defined syntenic) gene divergence summary\n')
        w.write('====================================================================\n')
        w.write(f'Anchors file: {ANCHORS}\n')
        w.write(f'H1 inversion-overlap genes: {INV_H1_GENES}\n')
        w.write(f'H1 inversion-overlap gene roots (input): {len(inv_roots)}\n')
        w.write(f'All anchor pairs: {len(anchors)}\n')
        w.write(f'Unique one-to-one anchor pairs used (H1 and H2 unique in anchors): {len(unique_pairs)}\n')
        w.write(f'Pairs with status=ok: {len(ok)}\n')
        w.write(f'Inversion-overlap anchor pairs (status=ok): {len(inv)}\n')
        w.write(f'Non-inversion anchor pairs (status=ok): {len(col)}\n')
        w.write(f'Unique inversion H1 genes represented in anchor-paired set: {len(inv_genes_anchored)} / {len(inv_roots)}\n')
        w.write(f'Inversion H1 genes not found in unique anchor-paired set: {len(inv_roots_not_anchored)}\n')
        if inv_roots_not_anchored:
            w.write('First 30 inversion genes not anchor-paired (unique anchors): ' + ','.join(inv_roots_not_anchored[:30]) + '\n')
        w.write('\nStatus counts:\n')
        for k, v in df['status'].value_counts(dropna=False).items():
            w.write(f'  {k}: {v}\n')
        w.write('\ndN/dS status counts (within status=ok):\n')
        for k, v in ok['dnds_status'].value_counts(dropna=False).items():
            w.write(f'  {k}: {v}\n')
        w.write('\nMetric comparison (Mann-Whitney U, two-sided):\n')
        for _, r in summary_df.iterrows():
            w.write(
                f"  {r['metric']}: inv n={int(r['inv_n'])}, median={fmt(r['inv_median'])}, mean={fmt(r['inv_mean'])}; "
                f"col n={int(r['col_n'])}, median={fmt(r['col_median'])}, mean={fmt(r['col_mean'])}; "
                f"p={fmt(r['pvalue_two_sided'])}\n"
            )
        w.write('\nTop 10 inversion pairs by CDS nucleotide p-distance:\n')
        cols_show = ['H1_gene', 'H2_gene', 'cds_nt_pdistance', 'aa_pdistance', 'dS_NG86', 'dN_NG86', 'aligned_codon_sites']
        top = inv_rank_nt.head(10)[[c for c in cols_show if c in inv_rank_nt.columns]]
        w.write(top.to_string(index=False))
        w.write('\n')

    print(f'[INFO] wrote pair table: {out_tsv}')
    print(f'[INFO] wrote metric summary: {summary_tsv}')
    print(f'[INFO] wrote text summary: {txt}')
    print(f'[INFO] wrote top tables: {top_nt_tsv}, {top_ds_tsv}')


if __name__ == '__main__':
    main()
