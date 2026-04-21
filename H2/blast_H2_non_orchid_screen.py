#!/g/data/xf3/miniconda/bin/python
import argparse
import concurrent.futures
import os
import sqlite3
import subprocess
from collections import defaultdict
from pathlib import Path


ORCHIDACEAE_TAXID = 4747
VIRIDIPLANTAE_TAXID = 33090


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Screen an H2 genome against nt, summarize taxonomy for each contig, "
            "and write a filtered FASTA copy without candidate non-plant contigs."
        )
    )
    parser.add_argument("--genome", required=True, help="Input genome FASTA.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--blastn-bin", required=True, help="Path to blastn binary.")
    parser.add_argument("--blast-db", required=True, help="BLAST nt database prefix.")
    parser.add_argument(
        "--taxonomy-sqlite",
        required=True,
        help="SQLite DB with TaxidInfo(taxid,parent).",
    )
    parser.add_argument("--threads", type=int, default=48, help="Total CPU budget.")
    parser.add_argument(
        "--parallel-jobs",
        type=int,
        default=8,
        help="Number of blastn processes to run in parallel.",
    )
    parser.add_argument(
        "--min-len",
        type=int,
        default=5000,
        help="Skip contigs shorter than this length.",
    )
    parser.add_argument(
        "--query-limit",
        type=int,
        default=100000,
        help="For contigs >= this size, blast only the first N bp.",
    )
    parser.add_argument(
        "--max-target-seqs",
        type=int,
        default=20,
        help="Maximum target sequences to keep per query.",
    )
    parser.add_argument(
        "--remove-scope",
        choices=["none", "non_plant", "non_orchid"],
        default="non_plant",
        help=(
            "Automatic removal scope for the filtered FASTA copy. "
            "'non_plant' is the safe default."
        ),
    )
    return parser.parse_args()


def iter_fasta(path):
    header = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def write_fasta(records, path, width=60):
    with open(path, "w") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for i in range(0, len(seq), width):
                handle.write(seq[i : i + width] + "\n")


def split_queries(genome, split_dir, min_len, query_limit):
    split_dir.mkdir(parents=True, exist_ok=True)
    contigs = []
    for contig_id, seq in iter_fasta(genome):
        seq_len = len(seq)
        entry = {
            "contig": contig_id,
            "length": seq_len,
            "query_id": contig_id,
            "query_len": seq_len,
            "query_path": None,
            "skipped": False,
        }
        if seq_len < min_len:
            entry["skipped"] = True
            contigs.append(entry)
            continue

        if seq_len >= query_limit:
            query_id = f"{contig_id}_first100k"
            query_seq = seq[:query_limit]
        else:
            query_id = contig_id
            query_seq = seq

        query_path = split_dir / f"{query_id}.fasta"
        if not query_path.exists() or query_path.stat().st_size == 0:
            write_fasta([(query_id, query_seq)], query_path)

        entry["query_id"] = query_id
        entry["query_len"] = len(query_seq)
        entry["query_path"] = str(query_path)
        contigs.append(entry)
    return contigs


def run_one_blast(task):
    query_path = task["query_path"]
    output_path = task["output_path"]
    blastn_bin = task["blastn_bin"]
    blast_db = task["blast_db"]
    blast_threads = task["blast_threads"]
    max_target_seqs = task["max_target_seqs"]
    outfmt = (
        "6 qseqid sseqid stitle pident length mismatch gapopen "
        "qstart qend sstart send sscinames staxids evalue bitscore"
    )

    if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
        return query_path, "skipped_existing"

    env = os.environ.copy()
    blastdb_dir = os.path.dirname(blast_db)
    current = env.get("BLASTDB", "")
    env["BLASTDB"] = f"{current}:{blastdb_dir}" if current else blastdb_dir

    cmd = [
        blastn_bin,
        "-db",
        blast_db,
        "-query",
        query_path,
        "-out",
        output_path,
        "-outfmt",
        outfmt,
        "-perc_identity",
        "75",
        "-evalue",
        "1e-5",
        "-max_target_seqs",
        str(max_target_seqs),
        "-num_threads",
        str(blast_threads),
    ]
    subprocess.run(cmd, check=True, env=env)
    return query_path, "done"


class TaxonomyLookup:
    def __init__(self, sqlite_path):
        self.conn = sqlite3.connect(sqlite_path)
        self.parent_cache = {}
        self.descendant_cache = {}

    def parent(self, taxid):
        if taxid in self.parent_cache:
            return self.parent_cache[taxid]
        row = self.conn.execute(
            "SELECT parent FROM TaxidInfo WHERE taxid = ?", (taxid,)
        ).fetchone()
        parent = row[0] if row else None
        self.parent_cache[taxid] = parent
        return parent

    def is_descendant(self, taxid, ancestor_taxid):
        key = (taxid, ancestor_taxid)
        if key in self.descendant_cache:
            return self.descendant_cache[key]

        original = taxid
        visited = set()
        result = False
        while taxid is not None and taxid not in visited:
            if taxid == ancestor_taxid:
                result = True
                break
            visited.add(taxid)
            parent = self.parent(taxid)
            if parent is None or parent == taxid:
                break
            taxid = parent

        self.descendant_cache[key] = result
        return result


def parse_first_taxid(raw):
    if not raw:
        return None
    for token in str(raw).replace(",", ";").split(";"):
        token = token.strip()
        if token.isdigit():
            return int(token)
    return None


def parse_blast_file(path, taxdb):
    hits = []
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return hits

    with open(path) as handle:
        for rank, line in enumerate(handle, start=1):
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 15:
                continue
            taxid = parse_first_taxid(cols[12])
            hit = {
                "rank": rank,
                "qseqid": cols[0],
                "sseqid": cols[1],
                "stitle": cols[2],
                "pident": float(cols[3]),
                "length": int(cols[4]),
                "mismatch": int(cols[5]),
                "gapopen": int(cols[6]),
                "qstart": int(cols[7]),
                "qend": int(cols[8]),
                "sstart": int(cols[9]),
                "send": int(cols[10]),
                "sscinames": cols[11],
                "staxids": cols[12],
                "taxid": taxid,
                "evalue": cols[13],
                "bitscore": float(cols[14]),
                "is_orchidaceae": False,
                "is_viridiplantae": False,
            }
            if taxid is not None:
                hit["is_orchidaceae"] = taxdb.is_descendant(taxid, ORCHIDACEAE_TAXID)
                hit["is_viridiplantae"] = taxdb.is_descendant(taxid, VIRIDIPLANTAE_TAXID)
            hits.append(hit)
    return hits


def classify_contig(entry, hits):
    top1 = hits[0] if hits else None
    has_orchidaceae = any(h["is_orchidaceae"] for h in hits)
    has_viridiplantae = any(h["is_viridiplantae"] for h in hits)

    if entry["skipped"]:
        classification = "skipped_short"
    elif not hits:
        classification = "no_hit_manual_review"
    elif has_orchidaceae:
        classification = "keep_orchidaceae"
    elif has_viridiplantae:
        classification = "manual_review_non_orchid_plant"
    else:
        classification = "candidate_remove_non_plant"

    summary = {
        "contig": entry["contig"],
        "contig_len": entry["length"],
        "query_id": entry["query_id"],
        "query_len": entry["query_len"],
        "classification": classification,
        "has_hits": bool(hits),
        "has_orchidaceae_hit": has_orchidaceae,
        "has_viridiplantae_hit": has_viridiplantae,
        "top1_sscinames": top1["sscinames"] if top1 else "",
        "top1_taxid": top1["taxid"] if top1 and top1["taxid"] is not None else "",
        "top1_pident": top1["pident"] if top1 else "",
        "top1_bitscore": top1["bitscore"] if top1 else "",
        "top1_stitle": top1["stitle"] if top1 else "",
    }
    return summary


def write_tsv(path, header, rows):
    with open(path, "w") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(col, "")) for col in header) + "\n")


def decide_auto_remove(contig_rows, scope):
    to_remove = []
    for row in contig_rows:
        cls = row["classification"]
        if scope == "none":
            continue
        if scope == "non_plant" and cls == "candidate_remove_non_plant":
            to_remove.append(row["contig"])
        elif scope == "non_orchid" and cls in {
            "candidate_remove_non_plant",
            "manual_review_non_orchid_plant",
        }:
            to_remove.append(row["contig"])
    return sorted(set(to_remove))


def filter_genome(genome, remove_ids, kept_out, removed_out):
    remove_set = set(remove_ids)
    kept_records = []
    removed_records = []
    for contig_id, seq in iter_fasta(genome):
        if contig_id in remove_set:
            removed_records.append((contig_id, seq))
        else:
            kept_records.append((contig_id, seq))
    write_fasta(kept_records, kept_out)
    write_fasta(removed_records, removed_out)


def main():
    args = parse_args()

    genome = Path(args.genome)
    outdir = Path(args.outdir)
    split_dir = outdir / "split_for_blast"
    blast_dir = outdir / "Blast"
    summary_dir = outdir / "summary"

    outdir.mkdir(parents=True, exist_ok=True)
    split_dir.mkdir(parents=True, exist_ok=True)
    blast_dir.mkdir(parents=True, exist_ok=True)
    summary_dir.mkdir(parents=True, exist_ok=True)

    contigs = split_queries(genome, split_dir, args.min_len, args.query_limit)
    active = [c for c in contigs if not c["skipped"]]

    parallel_jobs = max(1, min(args.parallel_jobs, len(active) or 1, args.threads))
    blast_threads = max(1, args.threads // parallel_jobs)

    print(f"[INFO] Genome: {genome}")
    print(f"[INFO] Output dir: {outdir}")
    print(f"[INFO] Active contigs to blast: {len(active)}")
    print(f"[INFO] Parallel blast jobs: {parallel_jobs}")
    print(f"[INFO] Threads per blast job: {blast_threads}")
    print(f"[INFO] Remove scope: {args.remove_scope}")

    tasks = []
    for entry in active:
        tasks.append(
            {
                "query_path": entry["query_path"],
                "output_path": str(blast_dir / f"{entry['query_id']}.blast"),
                "blastn_bin": args.blastn_bin,
                "blast_db": args.blast_db,
                "blast_threads": blast_threads,
                "max_target_seqs": args.max_target_seqs,
            }
        )

    with concurrent.futures.ThreadPoolExecutor(max_workers=parallel_jobs) as pool:
        futures = [pool.submit(run_one_blast, task) for task in tasks]
        for future in concurrent.futures.as_completed(futures):
            query_path, status = future.result()
            print(f"[BLAST] {status}: {query_path}")

    taxdb = TaxonomyLookup(args.taxonomy_sqlite)
    all_hits = []
    contig_rows = []

    for entry in contigs:
        blast_path = blast_dir / f"{entry['query_id']}.blast"
        hits = []
        if not entry["skipped"]:
            hits = parse_blast_file(blast_path, taxdb)
        contig_rows.append(classify_contig(entry, hits))
        for hit in hits:
            all_hits.append(
                {
                    "contig": entry["contig"],
                    "contig_len": entry["length"],
                    "query_id": entry["query_id"],
                    "rank": hit["rank"],
                    "sscinames": hit["sscinames"],
                    "taxid": hit["taxid"] if hit["taxid"] is not None else "",
                    "is_viridiplantae": int(hit["is_viridiplantae"]),
                    "is_orchidaceae": int(hit["is_orchidaceae"]),
                    "pident": hit["pident"],
                    "length": hit["length"],
                    "bitscore": hit["bitscore"],
                    "evalue": hit["evalue"],
                    "stitle": hit["stitle"],
                    "sseqid": hit["sseqid"],
                }
            )

    hit_header = [
        "contig",
        "contig_len",
        "query_id",
        "rank",
        "sscinames",
        "taxid",
        "is_viridiplantae",
        "is_orchidaceae",
        "pident",
        "length",
        "bitscore",
        "evalue",
        "stitle",
        "sseqid",
    ]
    contig_header = [
        "contig",
        "contig_len",
        "query_id",
        "query_len",
        "classification",
        "has_hits",
        "has_orchidaceae_hit",
        "has_viridiplantae_hit",
        "top1_sscinames",
        "top1_taxid",
        "top1_pident",
        "top1_bitscore",
        "top1_stitle",
    ]
    write_tsv(summary_dir / "blast_hits_annotated.tsv", hit_header, all_hits)
    write_tsv(summary_dir / "contig_classification.tsv", contig_header, contig_rows)

    candidate_remove_non_plant = sorted(
        row["contig"]
        for row in contig_rows
        if row["classification"] == "candidate_remove_non_plant"
    )
    manual_review_non_orchid = sorted(
        row["contig"]
        for row in contig_rows
        if row["classification"] == "manual_review_non_orchid_plant"
    )
    no_hit_review = sorted(
        row["contig"]
        for row in contig_rows
        if row["classification"] == "no_hit_manual_review"
    )

    (summary_dir / "candidate_remove_non_plant.txt").write_text(
        "\n".join(candidate_remove_non_plant) + ("\n" if candidate_remove_non_plant else "")
    )
    (summary_dir / "manual_review_non_orchid_plant.txt").write_text(
        "\n".join(manual_review_non_orchid) + ("\n" if manual_review_non_orchid else "")
    )
    (summary_dir / "no_hit_manual_review.txt").write_text(
        "\n".join(no_hit_review) + ("\n" if no_hit_review else "")
    )

    auto_remove = decide_auto_remove(contig_rows, args.remove_scope)
    (summary_dir / f"auto_remove_{args.remove_scope}.txt").write_text(
        "\n".join(auto_remove) + ("\n" if auto_remove else "")
    )

    base_name = genome.name.rsplit(".fa", 1)[0]
    kept_out = outdir / f"{base_name}.filtered_{args.remove_scope}.fa"
    removed_out = outdir / f"{base_name}.removed_{args.remove_scope}.fa"
    filter_genome(genome, auto_remove, kept_out, removed_out)

    print(f"[INFO] Annotated hits: {summary_dir / 'blast_hits_annotated.tsv'}")
    print(f"[INFO] Contig classification: {summary_dir / 'contig_classification.tsv'}")
    print(f"[INFO] Candidate non-plant removals: {summary_dir / 'candidate_remove_non_plant.txt'}")
    print(f"[INFO] Manual review non-orchid plants: {summary_dir / 'manual_review_non_orchid_plant.txt'}")
    print(f"[INFO] No-hit contigs for manual review: {summary_dir / 'no_hit_manual_review.txt'}")
    print(f"[INFO] Auto-remove list ({args.remove_scope}): {summary_dir / f'auto_remove_{args.remove_scope}.txt'}")
    print(f"[INFO] Filtered FASTA: {kept_out}")
    print(f"[INFO] Removed FASTA: {removed_out}")


if __name__ == "__main__":
    main()
