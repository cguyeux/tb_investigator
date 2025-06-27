#!/usr/bin/env python3
"""Run Prodigal on a FASTA file and summarize predicted CDS."""
import argparse
from typing import List, Dict
from analyse_seq import run_prodigal, parse_gff, summarize_cds


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
