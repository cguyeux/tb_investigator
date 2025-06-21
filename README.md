# TB Investigator Utilities

This repository provides small scripts used in tuberculosis genome analyses.

## Preprocess reads

`preprocess_reads.py` trims adapter sequences and filters low-quality reads
using **cutadapt**. Install cutadapt and run the script before assembly to
improve the quality of your data.

### Installation

```bash
pip install cutadapt
```

### Usage

```bash
./preprocess_reads.py -1 sample_R1.fastq -2 sample_R2.fastq \
    -o trimmed/sample --quality 20 --length 50 --threads 4
```

This will produce `trimmed/sample_1.fastq` and `trimmed/sample_2.fastq` with
adapter sequences removed and low-quality bases trimmed.
