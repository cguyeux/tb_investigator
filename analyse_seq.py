#!/usr/bin/env python3
"""Analyse sequences and check their origin.

The script computes the GC content of an input FASTA file and can optionally
test whether the contigs match known plasmid, phage or transposon sequences.
Detection relies on external tools and databases such as BLAST, PLSDB,
PHASTER, ISEScan and TransposonPSI.
"""
import argparse
import os
import subprocess
from typing import List


def compute_gc(fasta: str) -> float:
    """Return the overall GC percentage of all sequences in *fasta*."""
    total_len = 0
    gc_count = 0
    seq = []
    with open(fasta) as fh:
        for line in fh:
            if line.startswith('>'):
                if seq:
                    sequence = ''.join(seq).upper()
                    gc_count += sequence.count('G') + sequence.count('C')
                    total_len += len(sequence)
                    seq = []
            else:
                seq.append(line.strip())
        if seq:
            sequence = ''.join(seq).upper()
            gc_count += sequence.count('G') + sequence.count('C')
            total_len += len(sequence)
    if total_len == 0:
        return 0.0
    return gc_count / total_len * 100


def blast_hits(fasta: str, db: str, evalue: float = 1e-5) -> List[str]:
    """Return query identifiers with hits against *db* using BLAST."""
    cmd = [
        "blastn",
        "-query",
        fasta,
        "-db",
        db,
        "-max_target_seqs",
        "1",
        "-evalue",
        str(evalue),
        "-outfmt",
        "6 qseqid sseqid pident length bitscore evalue",
    ]
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True
        )
    except FileNotFoundError as exc:
        raise RuntimeError("blastn not found. Install BLAST+ to use this option.") from exc

    hits = []
    for line in result.stdout.strip().splitlines():
        parts = line.split("\t")
        if parts:
            hits.append(parts[0])
    return hits


def run_isescan(fasta: str, isescan: str, outdir: str) -> bool:
    """Run ISEScan on *fasta*. Return True if predictions were produced."""
    os.makedirs(outdir, exist_ok=True)
    try:
        subprocess.run([
            isescan,
            fasta,
            outdir,
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except FileNotFoundError as exc:
        raise RuntimeError(f"{isescan} not found. Install ISEScan to use this option.") from exc
    gff = os.path.join(outdir, "prediction.gff3")
    return os.path.exists(gff)


def run_transposonpsi(fasta: str, script: str, outdir: str) -> bool:
    """Run TransposonPSI on *fasta*. Return True if prediction files were created."""
    os.makedirs(outdir, exist_ok=True)
    try:
        subprocess.run([
            script,
            fasta,
        ], cwd=outdir, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"{script} not found. Install TransposonPSI to use this option."
        ) from exc
    for fname in os.listdir(outdir):
        if fname.endswith("_transposonPSI.gff3"):
            return True
    return False


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Compute GC content and optionally search for plasmid, phage or "
            "transposon sequences"
        )
    )
    parser.add_argument("fasta", help="Input FASTA file")
    parser.add_argument("--plasmid-db", help="BLAST database built from PLSDB")
    parser.add_argument("--phage-db", help="BLAST database of phage sequences")
    parser.add_argument("--isescan", help="Path to isescan.py executable")
    parser.add_argument(
        "--transposonpsi", help="Path to TransposonPSI.pl script"
    )
    parser.add_argument(
        "--tmpdir",
        default="tmp",
        help="Directory for intermediate files from external tools",
    )
    parser.add_argument(
        "--evalue",
        type=float,
        default=1e-5,
        help="E-value threshold for BLAST searches",
    )
    args = parser.parse_args()

    gc = compute_gc(args.fasta)
    print(f"GC content: {gc:.2f}%")

    if args.plasmid_db:
        try:
            hits = blast_hits(args.fasta, args.plasmid_db, args.evalue)
        except RuntimeError as err:
            print(err)
        else:
            if hits:
                print("Possible plasmid origin for:", ", ".join(hits))
            else:
                print("No plasmid match found.")

    if args.phage_db:
        try:
            hits = blast_hits(args.fasta, args.phage_db, args.evalue)
        except RuntimeError as err:
            print(err)
        else:
            if hits:
                print("Possible phage origin for:", ", ".join(hits))
            else:
                print("No phage match found.")

    if args.isescan:
        try:
            found = run_isescan(args.fasta, args.isescan, os.path.join(args.tmpdir, "isescan"))
            msg = "ISEScan detected IS elements." if found else "ISEScan found no IS elements."
            print(msg)
        except RuntimeError as err:
            print(err)

    if args.transposonpsi:
        try:
            found = run_transposonpsi(
                args.fasta,
                args.transposonpsi,
                os.path.join(args.tmpdir, "transposonpsi"),
            )
            msg = (
                "TransposonPSI detected transposons." if found else "TransposonPSI found no transposons."
            )
            print(msg)
        except RuntimeError as err:
            print(err)


if __name__ == '__main__':
    main()
