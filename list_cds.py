#!/usr/bin/env python3
"""Run Prodigal on a FASTA file and summarize predicted CDS."""
import argparse
import subprocess
from typing import List, Dict


def run_prodigal(fasta: str, prefix: str, mode: str = "meta") -> tuple[str, str]:
    """Run prodigal and return paths to GFF and protein FASTA."""
    gff = f"{prefix}.gff"
    faa = f"{prefix}.faa"
    subprocess.run([
        "prodigal",
        "-i",
        fasta,
        "-p",
        mode,
        "-a",
        faa,
        "-f",
        "gff",
        "-o",
        gff,
        "-q",
    ], check=True)
    return gff, faa


def parse_gff(gff: str) -> List[Dict[str, str]]:
    """Parse Prodigal GFF and return CDS info."""
    cds = []
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            seqid, _source, type_, start, end, _score, strand, _phase, attrs = parts
            if type_ != "CDS":
                continue
            attr_dict: Dict[str, str] = {}
            for item in attrs.split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    attr_dict[key] = value
            cds.append({
                "seqid": seqid,
                "start": start,
                "end": end,
                "strand": strand,
                "id": attr_dict.get("ID", ""),
            })
    return cds


def summarize_cds(cds: List[Dict[str, str]]) -> None:
    """Print a short summary of CDS."""
    print(f"{len(cds)} CDS predicted")
    strand_counts: Dict[str, int] = {}
    for c in cds:
        strand_counts[c["strand"]] = strand_counts.get(c["strand"], 0) + 1
    for strand, count in strand_counts.items():
        print(f"{count} CDS on strand {strand}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Prodigal and list CDS")
    parser.add_argument("fasta", help="Input nucleotide FASTA")
    parser.add_argument(
        "-o", "--out-prefix", default="prodigal", help="Output prefix"
    )
    parser.add_argument(
        "-m",
        "--mode",
        choices=["meta", "single"],
        default="meta",
        help="Prodigal mode",
    )
    args = parser.parse_args()

    gff, faa = run_prodigal(args.fasta, args.out_prefix, args.mode)
    cds = parse_gff(gff)
    summarize_cds(cds)
    print(f"CDS annotations written to {gff}")
    print(f"Protein translations written to {faa}")


if __name__ == "__main__":
    main()
