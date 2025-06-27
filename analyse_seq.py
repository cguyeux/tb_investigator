#!/usr/bin/env python3
"""Compute GC content of a FASTA file."""
import argparse


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


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute GC content of a FASTA file"
    )
    parser.add_argument('fasta', help='Input FASTA file')
    args = parser.parse_args()
    gc = compute_gc(args.fasta)
    print(f"{gc:.2f}")


if __name__ == '__main__':
    main()
