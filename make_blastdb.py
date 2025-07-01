import os
import pickle
import subprocess
import unicodedata
from typing import List, Optional, Tuple


def preprocess_nanopore(read_fastq: str, out_prefix: str, threads: int = 16) -> str:
    """Trim ONT adapters and filter low-quality reads.

    The cleaned FASTQ path is returned.
    """
    trimmed = f"{out_prefix}_trimmed.fastq"
    cleaned = f"{out_prefix}_1.fastq"
    subprocess.run(
        f"porechop -i {read_fastq} -o {trimmed} -t {threads}",
        shell=True,
        check=True,
    )
    subprocess.run(
        f"filtlong --min_length 1000 --min_mean_q 7 {trimmed} > {cleaned}",
        shell=True,
        check=True,
    )
    return cleaned


def detect_fastq_files(srr: str) -> Tuple[str, Optional[str]]:
    """Return paths to raw FASTQ files for ``srr``.

    The function deals with varying naming conventions observed after
    ``fasterq-dump --split-files``. It looks for paired ``_1``/``_2`` files,
    single ``.fastq`` files or cases where the second mate is labelled
    ``_4.fastq``.
    """

    r1 = os.path.join("fastq", f"{srr}_1.fastq")
    r2 = os.path.join("fastq", f"{srr}_2.fastq")

    if os.path.exists(r1) and os.path.exists(r2):
        return r1, r2

    alt2 = os.path.join("fastq", f"{srr}_4.fastq")
    if os.path.exists(r1) and os.path.exists(alt2):
        return r1, alt2

    single = os.path.join("fastq", f"{srr}.fastq")
    if os.path.exists(single):
        return single, None

    if os.path.exists(r1):
        return r1, None

    raise FileNotFoundError(f"No FASTQ files found for {srr}")

import pandas as pd
from Bio import SeqIO


def get_read_length(fastq_file: str) -> int:
    """Return the length of the first read found in ``fastq_file``."""
    with open(fastq_file) as fh:
        fh.readline()
        seq = fh.readline().strip()
    return len(seq)


def choose_k_values(fastq_file: str) -> str:
    """Select SPAdes -k values according to read length."""
    length = get_read_length(fastq_file)
    candidates = [21, 33, 55, 77, 99, 127]
    selected = [k for k in candidates if k < length]
    return ",".join(str(k) for k in selected)


def sanitize_header(header: str) -> str:
    """Return an ASCII-only header safe for makeblastdb."""
    header = unicodedata.normalize("NFKD", header)
    header = header.encode("ascii", "ignore").decode("ascii")
    header = header.replace(" ", "_").replace("/", "_")
    return header


def deduplicate_fasta(fasta_path: str) -> None:
    """Ensure sequence identifiers are unique in ``fasta_path``."""

    if not os.path.exists(fasta_path):
        return

    entries = []
    header = None
    seq_lines = []

    with open(fasta_path) as fh:
        last = None
        for line in fh:
            last = line
            if line.startswith(">"):
                if header is not None:
                    entries.append((header, "".join(seq_lines)))
                header = sanitize_header(line[1:].split()[0])
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            entries.append((header, "".join(seq_lines)))
        if last is not None and not last.endswith("\n"):
            entries[-1] = (entries[-1][0], entries[-1][1] + "\n")

    counts = {}
    with open(fasta_path, "w") as out:
        for head, seq in entries:
            idx = counts.get(head, 0)
            counts[head] = idx + 1
            clean = head if idx == 0 else f"{head}_{idx+1}"
            out.write(f">{clean}\n")
            out.write(seq)


def build_rd_fasta(rd_dir: str, rd_fasta: str) -> None:
    """Concatenate all RD FASTA files found in ``rd_dir``.
    Headers are sanitized to avoid non-ASCII characters."""

    counts = {}
    with open(rd_fasta, "w") as out:
        for fname in sorted(os.listdir(rd_dir)):
            if not fname.endswith(".fasta"):
                continue
            path = os.path.join(rd_dir, fname)
            last = None
            with open(path) as fh:
                for line in fh:
                    last = line
                    if line.startswith(">"):
                        base = sanitize_header(line[1:].strip())
                        idx = counts.get(base, 0)
                        counts[base] = idx + 1
                        clean = base if idx == 0 else f"{base}_{idx+1}"
                        out.write(f">{clean}\n")
                    else:
                        out.write(line)
            if last is not None and not last.endswith("\n"):
                out.write("\n")


def analyze_absent_regions(srr: str, bam: str, h37rv_fasta: str, rd_dir: str,
                           results_dir: str = "results") -> None:
    """Extract H37Rv regions absent from ``bam`` and BLAST them against RDs."""

    os.makedirs(results_dir, exist_ok=True)
    zero_cov = os.path.join(results_dir, f"{srr}_zero_coverage.bed")
    merged_bed = os.path.join(results_dir, f"{srr}_zero_coverage_merged.bed")
    absent_fasta = os.path.join(results_dir, f"{srr}_absent_regions.fasta")
    rd_fasta = os.path.join(results_dir, "RD_all.fasta")
    blast_db = os.path.join(results_dir, "RD_db")
    blast_out = os.path.join(results_dir, f"{srr}_blast.txt")
    csv_out = os.path.join(results_dir, f"{sra_list[srr].replace(' ', '')}_{srr}_rd.csv")

    # Coverage and regions absent from H37Rv
    # Extract zero coverage positions as a proper BED file
    subprocess.run(
        f"bedtools genomecov -d -ibam {bam} | "
        "awk '$3==0 {printf(\"%s\\t%d\\t%d\\n\", $1, $2-1, $2)}' "
        f"> {zero_cov}",
        shell=True,
        check=True,
    )
    subprocess.run(f"bedtools merge -i {zero_cov} > {merged_bed}", shell=True, check=True)
    subprocess.run(
        f"bedtools getfasta -fi {h37rv_fasta} -bed {merged_bed} -fo {absent_fasta}",
        shell=True,
        check=True,
    )

    # Prepare RD database
    build_rd_fasta(rd_dir, rd_fasta)
    subprocess.run(
        f"makeblastdb -in {rd_fasta} -dbtype nucl -out {blast_db} -parse_seqids",
        shell=True,
        check=True,
    )
    subprocess.run(
        f"blastn -query {absent_fasta} -db {blast_db} "
        "-outfmt '6 qseqid sseqid pident length bitscore evalue qstart qend sstart send' "
        f"-max_target_seqs 1 -evalue 1e-5 -out {blast_out}",
        shell=True,
        check=True,
    )

    cols = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "bitscore",
        "evalue",
        "qstart",
        "qend",
        "sstart",
        "send",
    ]
    blast_df = pd.read_csv(blast_out, sep="\t", names=cols)
    blast_df.rename(columns={"length": "aln_length"}, inplace=True)

    bed_df = pd.read_csv(merged_bed, sep="\t", names=["chrom", "start", "end"])
    bed_df["qseqid"] = [rec.id for rec in SeqIO.parse(absent_fasta, "fasta")]
    bed_df["region_length"] = bed_df["end"] - bed_df["start"]

    final_df = bed_df.merge(blast_df, on="qseqid", how="left")
    col_order = [
        "chrom",
        "start",
        "end",
        "region_length",
        "qseqid",
        "sseqid",
        "pident",
        "aln_length",
        "bitscore",
        "evalue",
        "qstart",
        "qend",
        "sstart",
        "send",
    ]
    final_df = final_df[col_order]
    final_df.to_csv(csv_out, index=False)

sra_list = {
    "SRR32024533": "M. riyadhense", 
    "SRR30899481": "M. shinjukuense", # ok
    "DRR161275": "M. shinjukuense", # ok
    "DRR161212": "M. shinjukuense", # ok
    "SRR30899492": "M. lacus", # ok
    "DRR161283": "M. lacus", # ok
    "DRR161220": "M. lacus", # ok
    "SRR14802015": "M. decipiens", # ok
    "SRR12650503": "Canettii", # ok
    "ERR5159441": "Canettii", # ok
    "SRR12680971": "Canettii", # ok
    "ERR7451300": "Canettii", # ok
    "ERR552306": "Canettii", # ok
    "SRR10828835": "L8", # ok
    "SRR1173284": "L8", # ok
    "ERR12115321": "L8", # ok
    "ERR2513437": "L1", # ok
    "ERR4822123": "L1", # ok
    "ERR2512595": "L1", # ok
    "ERR4822131": "L1", # ok
    "ERR4799389": "L1", # ok
    "ERR756346": "L7", # ok
    "ERR234681": "L7", # ok
    "ERR1971864": "L7", # ok
    "ERR040140": "L3", # ok
    "ERR2514441": "L3", # ok
    "ERR3077965": "L3", # ok
    "ERR2513449": "L3", # ok
    "ERR4816244": "L3", # ok
    "ERR4831956": "L2.1 proto", # ok
    "ERR4831637": "L2.1 proto", # ok
    "SRR5065500": "L2.2", # ok
    "ERR551809": "L2.2", # ok 
    "ERR4810586": "L2.2", # ok
    "ERR234202": "L5", # ok
    "ERR5336143": "L5",
    "ERR1082132": "L5", # ok
    "ERR6866750": "L4.1",
    "SRR26112465": "L4.1", # ok
    "ERR11466175": "L4.1", # ok
    "SRR16087965": "L4.1", # ok
    "SRR23861992": "L4.2", # ok
    "SRR3082077": "L4.2", # ok
    "ERR5865631": "L4.2", # ok
    "SRR28498548": "L4.2", # ok
    "ERR2791547": "L4.17", # ok
    "ERR1193816": "L4.17", # ok
    "ERR4821894": "L4.17", # ok
    "SRR21678367": "L4.18", # ok
    "ERR4822579": "L4.18", # ok
    "ERR2516366": "L4.18", # ok
    "SRR4033537": "L4.3", # ok
    "ERR2503415": "L4.3", # ok
    "ERR11488624": "L4.3", # ok
    "SRR18337407": "L4.3", # ok
    "ERR4827827": "L4.5 pre", # ok
    "ERR133829": "L4.5 pre", # ok
    "ERR5865010": "L4.5 pre", # ok
    "ERR4811090": "L4.5", # ok
    "ERR4831031": "L4.5", # ok
    "ERR4808923": "L4.5", # ok
    "SRR6824457": "L4.5", # ok
    "SRR15368984": "L4.14", # ok
    "SRR21864655": "L4.14", # ok
    "SRR18548992": "L4.14", # ok
    "ERR2517475": "L4.13", # ok
    "ERR2516176": "L4.13", # ok
    "ERR550924": "L4.13", # ok
    "ERR4813747": "L4.4", # ok
    "SRR14519870": "L4.4", # ok
    "ERR8975764": "L4.4", # ok
    "SRR4037728": "L4.4", # ok
    "ERR2707060": "L4.6.1", # ok
    "ERR1035303": "L4.6.1", # ok
    "ERR6358988": "L4.6.1", # ok
    "SRR26112020": "L4.6.2", # ok
    "SRR7517773": "L4.6.2", # ok
    "ERR551552": "L4.6.2", # ok
    "SRR19276639": "L4.11", # ok
    "ERR4553824": "L4.11", # ok
    "ERR4457827": "L4.11", # ok
    "SRR16370260": "L4.12", # ok
    "ERR1034925": "L4.12", # ok
    "ERR161019": "L4.12",
    "ERR4821775": "L4.16",  # ok
    "ERR10465960": "L4.16",
    "ERR4815931": "L4.16",
    "SRR6480475": "L4.15", # ok
    "ERR4811456": "L4.15",
    "SRR16370338": "L4.15",
    "ERR7363138": "L4.7", # ok
    "ERR1465939": "L4.7",
    "ERR2653121": "L4.7",
    "ERR4831403": "L4.7",
    "ERR4818890": "L4.7",
    "ERR4829326": "L4.7",
    "SRR3085327": "L4.8",
    "SRR18337381": "L4.8", # ok
    "ERR4822085": "L4.8",
    "ERR11466257": "L4.9",
    "ERR2513858": "L4.9", # ok
    "ERR4553799": "L4.9",
    "SRR13736067": "L4.9 H37Rv", # ok
    "ERR5979827": "L4.9 H37Rv",
    "SRR23080319": "L4.9 H37Rv",
    "ERR552768": "Pinipedii", # ok
    "SRR1239339": "Pinipedii",
    "SRR7693090": "Pinipedii",
    "ERR552037": "Microti", # ok
    "ERR4627394": "Microti",
    "ERR2659166": "Microti",
    "SRR16058599": "Orygis La3", # ok
    "ERR2517593": "Orygis La3",
    "SRR10251199": "Orygis La3",
    "SRR7617526": "Bovis La1", # ok
    "ERR6337391": "Bovis La1",
    "SRR1791869": "Bovis La1",
    "SRR14782713": "Bovis La1.2 BCG", # ok
    "SRR17111243": "Caprae La2", # ok
    "ERR551704": "Caprae La2",
    "SRR13888777": "Caprae La2",
    "SRR3135069": "La4", # ok
    "DRR120408": "La4",
    "SRR16503843": "La4",
    "ERR2707158": "L10", # ok 
    "ERR2516384": "L10", # ok 
    "ERR4815138": "L6.1", # ok
    "ERR9787035": "L6.1",
    "ERR1023243": "L6.1",
    "ERR10680128": "L6.2", # ok
    "ERR3170426": "L6.2",
    "ERR3800960": "L6.2",
    "ERR3806601": "L6.3", # ok
    "ERR439992": "L6.3",
    "ERR1023286": "L6.3",
    "ERR181315": "L9", # ok
    "ERR2514391": "L9",
    "SRR24827312": "L9",
    "SRR3500411": "Mungi", # ok
    "SRR3745458": "Dassie", # ok
    "SRR3745456": "Dassie",
    "SRR3745457": "Dassie",
    "ERR713575": "Chimpanze",
    "ERR970412": "Suricattae", # ok
    "ERR970409": "Suricattae",
    "ERR970410": "Suricattae",
}

# bwa index H37Rv.fasta
os.system("mkdir -p fastq")
os.system("mkdir -p fastq/cleaned")
os.system("mkdir -p alignments")
os.system("mkdir -p unmapped")
os.system("mkdir -p assemblies_unmapped")
os.system("mkdir -p contigs_unmapped")
os.system("mkdir -p data/RD")
os.system("mkdir -p mapped")
os.system("mkdir -p assemblies_mapped")
os.system("mkdir -p contigs_mapped")
os.system("mkdir -p bdd")

def concat(fasta_file, motif_source='', motif_cible=''):
    with open(fasta_file, 'r') as f:
        fasta = f.read()
    if len(motif_source)>0:
        fasta = fasta.replace(motif_source, motif_cible)
    with open('data/all_contigs.fasta', 'a') as f:
        f.write('\n')
        f.write(fasta)

def append_to_fasta(fasta_file: str, dest: str, motif_source: str = "", motif_cible: str = "") -> None:
    """Append ``fasta_file`` content to ``dest`` applying header replacement."""
    if not os.path.exists(fasta_file):
        return
    with open(fasta_file, "r") as fh:
        content = fh.read()
    if motif_source:
        content = content.replace(motif_source, motif_cible)
    with open(dest, "a") as out:
        out.write("\n")
        out.write(content)

def update_lineage_db(srr: str, lineage: str, mapped: str, unmapped: str) -> None:
    """Create or update the BLAST DB for ``lineage`` with ``srr`` contigs."""

    clean = lineage.replace(" ", "")
    fasta_path = os.path.join("bdd", f"{clean}.fasta")
    header_lineage = sanitize_header(lineage)

    existing = ""
    if os.path.exists(fasta_path):
        with open(fasta_path) as fh:
            existing = fh.read()

    if f">{srr}_" not in existing:
        append_to_fasta(
            unmapped,
            fasta_path,
            motif_source=">NODE_",
            motif_cible=f">{srr}_{header_lineage}_unmapped_NODE_",
        )
        append_to_fasta(
            mapped,
            fasta_path,
            motif_source=">NODE_",
            motif_cible=f">{srr}_{header_lineage}_mapped_NODE_",
        )

    # Ensure each sequence ID is unique before running makeblastdb
    deduplicate_fasta(fasta_path)

    subprocess.run(
        f"makeblastdb -in {fasta_path} -dbtype nucl -out bdd/{clean} -parse_seqids",
        shell=True,
        check=True,
    )

with open('done.pkl', 'rb') as f:
     done = pickle.load(f)


for SRR in [u for u in sra_list if u not in done]:
    print(f" - On téléchage {SRR} ({sra_list[SRR]})")
    os.system(f"fasterq-dump -e 16 --split-files --force --outdir fastq {SRR}")

    try:
        raw_r1, raw_r2 = detect_fastq_files(SRR)
    except FileNotFoundError as e:
        print(e)
        continue

    paired = raw_r2 is not None
    print(f" - Nettoyage {SRR} ({sra_list[SRR]})")

    if paired:
        cmd = [
            "python3",
            "preprocess_reads.py",
            "-1",
            raw_r1,
            "-o",
            f"fastq/cleaned/{SRR}",
            "--threads",
            "16",
        ]
        cmd.extend(["-2", raw_r2])
        subprocess.run(cmd, check=True)
    else:
        preprocess_nanopore(
            read_fastq=raw_r1,
            out_prefix=f"fastq/cleaned/{SRR}",
            threads=16,
        )

    print(f" - On aligne {SRR} ({sra_list[SRR]})")
    if paired:
        os.system(
            f"bwa mem -t 16 data/H37Rv.fasta fastq/cleaned/{SRR}_1.fastq fastq/cleaned/{SRR}_2.fastq > alignments/{SRR}.sam"
        )
    else:
        os.system(
            f"bwa mem -x ont2d -t 16 data/H37Rv.fasta fastq/cleaned/{SRR}_1.fastq > alignments/{SRR}.sam"
        )
    os.system(f"samtools view -bS alignments/{SRR}.sam | samtools sort -o alignments/{SRR}_sorted.bam")
    os.system(f"rm -f alignments/{SRR}.sam")
    os.system(f"samtools index alignments/{SRR}_sorted.bam")

    analyze_absent_regions(
        srr=SRR,
        bam=f"alignments/{SRR}_sorted.bam",
        h37rv_fasta="data/H37Rv.fasta",
        rd_dir="data/RD",
    )

    if paired:
        # Extraction en BAM des reads non mappés (les deux mates non mappés)
        os.system(f"samtools view -b -f 12 -F 256 alignments/{SRR}_sorted.bam > unmapped/{SRR}_unmapped.bam")
        # Extraction en fastq directement avec samtools
        os.system(
            f"samtools fastq -1 unmapped/{SRR}_unmapped_1.fastq -2 unmapped/{SRR}_unmapped_2.fastq -0 /dev/null -s /dev/null -n unmapped/{SRR}_unmapped.bam"
        )
        outdir = f"assemblies_unmapped/{SRR}"
        contigs = f"contigs_unmapped/{SRR}_unmapped_contigs.fasta"
        k_values = choose_k_values(f"unmapped/{SRR}_unmapped_1.fastq")
        os.system(
            f"spades.py --careful -t 16 --cov-cutoff auto -k {k_values} -1 unmapped/{SRR}_unmapped_1.fastq -2 unmapped/{SRR}_unmapped_2.fastq -o {outdir}"
        )
    else:
        os.system(f"samtools view -b -f 4 alignments/{SRR}_sorted.bam > unmapped/{SRR}_unmapped.bam")
        os.system(f"samtools fastq unmapped/{SRR}_unmapped.bam > unmapped/{SRR}_unmapped_1.fastq")
        outdir = f"assemblies_unmapped/{SRR}"
        contigs = f"contigs_unmapped/{SRR}_unmapped_contigs.fasta"
        os.system(
            f"spades.py --careful -t 16 --cov-cutoff auto --nanopore unmapped/{SRR}_unmapped_1.fastq -o {outdir}"
        )
    assembled = f"{outdir}/contigs.fasta"
    if os.path.exists(assembled):
        os.rename(assembled, contigs)
    os.system(f"rm -f fastq/{SRR}*fastq")
    os.system(f"rm -f fastq/cleaned/{SRR}*fastq")

    sorted_bam = f"alignments/{SRR}_sorted.bam"
    
    # 1) extraction des reads mappés
    mapped_bam = f"mapped/{SRR}_mapped.bam"
    if not os.path.exists(mapped_bam):
        os.system(f"samtools view -b -F 4 {sorted_bam} > {mapped_bam}")

    mapped_fastq1 = f"mapped/{SRR}_mapped_1.fastq"
    mapped_fastq2 = f"mapped/{SRR}_mapped_2.fastq"
    if paired:
        if not (os.path.exists(mapped_fastq1) and os.path.exists(mapped_fastq2)):
            os.system(
                f"samtools fastq -1 {mapped_fastq1} -2 {mapped_fastq2} -0 /dev/null -s /dev/null -n {mapped_bam}"
            )
    else:
        if not os.path.exists(mapped_fastq1):
            os.system(f"samtools fastq {mapped_bam} > {mapped_fastq1}")

    # 3) assemblage SPAdes
    outdir_m = f"assemblies_mapped/{SRR}"
    contigs_m = f"contigs_mapped/{SRR}_mapped_contigs.fasta"
    if not os.path.exists(contigs_m):
        if paired:
            k_values = choose_k_values(mapped_fastq1)
            os.system(
                f"spades.py -k {k_values} -1 {mapped_fastq1} -2 {mapped_fastq2} -o {outdir_m}"
            )
        else:
            os.system(
                f"spades.py --careful -t 16 --cov-cutoff auto --nanopore {mapped_fastq1} -o {outdir_m}"
            )
        assembled_m = f"{outdir_m}/contigs.fasta"
        if os.path.exists(assembled_m):
            os.rename(assembled_m, contigs_m)

    update_lineage_db(
        srr=SRR,
        lineage=sra_list[SRR],
        mapped=contigs_m,
        unmapped=contigs,
    )

    os.system('rm -f data/all_contigs.fasta')

    for rep in ['data/sequences', 'data/sequences/IS', 'data/sequences/CDS']:
        for seq in os.listdir(rep):
            if seq.endswith('fasta'):
                concat(f'{rep}/{seq}')

    for contigs in os.listdir("contigs_unmapped"):
        sra = contigs.split('_')[0]
        concat(
            f"contigs_unmapped/{contigs}", 
            motif_source='>NODE_', 
            motif_cible=f'>{sra}_{sra_list[sra]}_unmapped_NODE_'
        )

    for contigs in os.listdir("contigs_mapped"):
        sra = contigs.split('_')[0]
        concat(
            f"contigs_mapped/{contigs}",
            motif_source='>NODE_',
            motif_cible=f'>{sra}_{sra_list[sra]}_mapped_NODE_'
        )
        
    concat("data/H37Rv.fasta")

    H37Rv = open('data/H37Rv.fasta').read()
    H37Rv = ''.join(H37Rv.split('\n')[1:])

    with open('data/rd.bed') as f:
        RDs = f.read()
        for rd in RDs.split('\n')[:-1]:
            _, debut, fin, nom = rd.split('\t')
            debut = int(debut)
            fin = int(fin)
            nom = nom.rstrip(' ')
            nom = unicodedata.normalize('NFKD', nom).encode('ascii', 'ignore').decode('ascii')
            nom = nom.replace(' ', '_').replace('/', '-')
            assert debut < fin
            with open(f'data/RD/{nom}.fasta', 'w') as f:
                f.write(f'>{nom}\n')
                f.write(H37Rv[debut-1:fin-1])    
            concat(f"data/RD/{nom}.fasta")
        

    # Sanitize headers and ensure unique identifiers before creating the
    # consolidated BLAST database
    deduplicate_fasta("data/all_contigs.fasta")
    os.system(
        "makeblastdb -in data/all_contigs.fasta -dbtype nucl -out bdd/mydb -parse_seqids"
    )
    done.append(SRR)
    with open('done.pkl', 'wb') as f:
        pickle.dump(done, f)
        exit()    

# blastn -query data/sequences/TbD1.fasta -db bdd/mydb
# blastn -query data/sequences/TbD1.fasta -db bdd/mydb -outfmt '6 qseqid sseqid pident length bitscore evalue qstart qend sstart send' -max_target_seqs 1 -evalue 1e-5
