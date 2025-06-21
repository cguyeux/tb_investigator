#!/usr/bin/env python3
"""Trim adapters and low-quality bases using cutadapt."""
import argparse
import subprocess


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Trim adapters and filter low-quality reads using cutadapt"
    )
    parser.add_argument(
        "-1",
        "--reads1",
        required=True,
        help="Input FASTQ file for read 1",
    )
    parser.add_argument(
        "-2",
        "--reads2",
        help="Input FASTQ file for read 2 (optional)",
    )
    parser.add_argument(
        "-o",
        "--output-prefix",
        required=True,
        help="Prefix for the trimmed output FASTQ files",
    )
    parser.add_argument(
        "--adapter-fwd",
        default="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        help="Adapter sequence for read 1",
    )
    parser.add_argument(
        "--adapter-rev",
        default="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        help="Adapter sequence for read 2",
    )
    parser.add_argument(
        "--quality",
        type=int,
        default=20,
        help="Trim low-quality bases below this Phred score",
    )
    parser.add_argument(
        "--length",
        type=int,
        default=50,
        help="Discard reads shorter than this length after trimming",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of threads to use",
    )

    args = parser.parse_args()

    cmd = [
        "cutadapt",
        "-j",
        str(args.threads),
        "-q",
        str(args.quality),
        "-m",
        str(args.length),
        "-a",
        args.adapter_fwd,
    ]

    if args.reads2:
        cmd += ["-A", args.adapter_rev]

    cmd += ["-o", f"{args.output_prefix}_1.fastq"]

    if args.reads2:
        cmd += ["-p", f"{args.output_prefix}_2.fastq", args.reads1, args.reads2]
    else:
        cmd += [args.reads1]

    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
