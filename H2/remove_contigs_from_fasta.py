#!/g/data/xf3/miniconda/bin/python
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="Remove contigs listed in a text file from a FASTA without touching the original."
    )
    parser.add_argument("--input", required=True, help="Input FASTA.")
    parser.add_argument("--remove-list", required=True, help="One contig ID per line.")
    parser.add_argument("--output", required=True, help="Filtered FASTA.")
    parser.add_argument("--removed-output", required=True, help="Removed contigs FASTA.")
    parser.add_argument("--summary", required=True, help="Summary TSV.")
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


def main():
    args = parse_args()

    remove_ids = {
        line.strip()
        for line in open(args.remove_list)
        if line.strip() and not line.startswith("#")
    }

    kept = []
    removed = []
    for contig_id, seq in iter_fasta(args.input):
        if contig_id in remove_ids:
            removed.append((contig_id, seq))
        else:
            kept.append((contig_id, seq))

    write_fasta(kept, args.output)
    write_fasta(removed, args.removed_output)

    removed_set = {name for name, _ in removed}
    missing = sorted(remove_ids - removed_set)
    with open(args.summary, "w") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"requested_remove_ids\t{len(remove_ids)}\n")
        handle.write(f"removed_contigs\t{len(removed)}\n")
        handle.write(f"kept_contigs\t{len(kept)}\n")
        handle.write(f"missing_remove_ids\t{len(missing)}\n")
        if missing:
            handle.write("missing_ids\t" + ",".join(missing) + "\n")


if __name__ == "__main__":
    main()
