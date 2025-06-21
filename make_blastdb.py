import os
import subprocess
import unicodedata


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

sra_list = {
    "SRR32024533": "M. riyadhense",
    "SRR30899481": "M. shinjukuense",
    "SRR30899492": "M. lacus",
    "SRR14802015": "M. decipiens",
    "SRR12650503": "Canettii",
    "ERR5159441": "Canettii",
    "SRR12680971": "Canettii",
    "ERR7451300": "Canettii",
    "ERR552306": "Canettii",
    "SRR10828835": "L8",
    "SRR1173284": "L8",
    "ERR12115321": "L8",
    "ERR756346": "L7",
    "ERR234681": "L7",
    "ERR1971864": "L7",    
    "ERR2513437": "L1",
    "ERR4822123": "L1",
    "ERR2512595": "L1",
    "ERR4822131": "L1",
    "ERR4799389": "L1",
    "ERR040140": "L3",
    "ERR2514441": "L3",
    "ERR3077965": "L3",
    "ERR2513449": "L3",
    "ERR4816244": "L3",
    "ERR4831956": "L2.1 proto",
    "SRR1710059": "L2.1 proto",
    "ERR4831637": "L2.1 proto",
    "SRR5065500": "L2.2",
    "ERR551809": "L2.2",
    "ERR019862": "L2.2",
    "SRR12420419": "L2.2",
    "ERR4810586": "L2.2",
    "ERR6866750": "L4.1",
    "SRR26112465": "L4.1",
    "ERR11466175": "L4.1",
    "SRR16087965": "L4.1",
    "SRR23861992": "L4.2",
    "SRR3082077": "L4.2",
    "ERR5865631": "L4.2",
    "SRR28498548": "L4.2",
    "ERR1193816": "L4.17",
    "ERR2791547": "L4.17",
    "ERR4821894": "L4.17",
    "ERR4822579": "L4.18",
    "SRR21678367": "L4.18",
    "ERR2516366": "L4.18",
    "ERR2503415": "L4.3",
    "SRR4033537": "L4.3",
    "ERR11488624": "L4.3",
    "SRR18337407": "L4.3",
    "ERR4827827": "L4.5 pre",
    "ERR133829": "L4.5 pre",
    "ERR5865010": "L4.5 pre",
    "ERR4811090": "L4.5",
    "ERR4831031": "L4.5",
    "ERR4808923": "L4.5",
    "SRR6824457": "L4.5",
    "SRR21864655": "L4.14",
    "SRR15368984": "L4.14",
    "SRR18548992": "L4.14",
    "ERR2516176": "L4.13",
    "ERR2517475": "L4.13",
    "ERR550924": "L4.13",
    "SRR14519870": "L4.4",
    "ERR4813747": "L4.4",
    "ERR8975764": "L4.4",
    "SRR4037728": "L4.4",
    "ERR2707060": "L4.6.1",
    "ERR1035303": "L4.6.1",
    "ERR6358988": "L4.6.1",
    "SRR7517773": "L4.6.2",
    "SRR26112020": "L4.6.2",
    "ERR551552": "L4.6.2",
    "SRR19276639": "L4.11",
    "ERR4553824": "L4.11",
    "ERR4457827": "L4.11",
    "SRR16370260": "L4.12",
    "ERR1034925": "L4.12",
    "ERR161019": "L4.12",
    "ERR4821775": "L4.16",
    "ERR10465960": "L4.16",
    "ERR4815931": "L4.16",
    "SRR6480475": "L4.15",
    "ERR4811456": "L4.15",
    "SRR16370338": "L4.15",
    "ERR7363138": "L4.7",
    "ERR1465939": "L4.7",
    "ERR2653121": "L4.7",
    "ERR4831403": "L4.7",
    "ERR4818890": "L4.7",
    "ERR4829326": "L4.7",
    "SRR18337381": "L4.8",
    "SRR3085327": "L4.8",
    "ERR4822085": "L4.8",
    "ERR2513858": "L4.9",
    "ERR11466257": "L4.9",
    "ERR4553799": "L4.9",
    "SRR13736067": "L4.9 H37Rv",
    "ERR5979827": "L4.9 H37Rv",
    "SRR23080319": "L4.9 H37Rv",
    "ERR234202": "L5",
    "ERR017801": "L5",
    "ERR5336143": "L5",
    "ERR1082132": "L5",
    "ERR552768": "Pinipedii",
    "SRR1239339": "Pinipedii",
    "SRR7693090": "Pinipedii",
    "ERR552037": "Microti",
    "ERR4627394": "Microti",
    "ERR2659166": "Microti",
    "SRR16058599": "Orygis La3",
    "ERR2517593": "Orygis La3",
    "SRR10251199": "Orygis La3",
    "SRR7617526": "Bovis La1",
    "ERR6337391": "Bovis La1",
    "SRR1791869": "Bovis La1",
    "SRR14782713": "Bovis La1.2 BCG",
    "SRR17111243": "Caprae La2",
    "ERR551704": "Caprae La2",
    "SRR13888777": "Caprae La2",
    "SRR3135069": "La4",
    "DRR120408": "La4",
    "SRR16503843": "La4",
    "ERR2707158": "L10",
    "ERR2516384": "L10",
    "ERR4815138": "L6.1",
    "ERR9787035": "L6.1",
    "ERR1023243": "L6.1",
    "ERR10680128": "L6.2",
    "ERR3170426": "L6.2",
    "ERR3800960": "L6.2",
    "ERR3806601": "L6.3",
    "ERR439992": "L6.3",
    "ERR1023286": "L6.3",
    "ERR181315": "L9",
    "ERR2514391": "L9",
    "SRR24827312": "L9",
    "SRR3500411": "Mungi",
    "SRR3745458": "Dassie",
    "SRR3745456": "Dassie",
    "SRR3745457": "Dassie",
    "ERR713575": "Chimpanze",
    "ERR970412": "Suricattae",
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

def concat(fasta_file, motif_source='', motif_cible=''):
    print(motif_source)
    print(motif_cible)
    with open(fasta_file, 'r') as f:
        fasta = f.read()
    if len(motif_source)>0:
        fasta = fasta.replace(motif_source, motif_cible)
    with open('data/all_contigs.fasta', 'a') as f:
        f.write('\n')
        f.write(fasta)

stop = False
for SRR in sra_list:
    if f"{SRR}_mapped_contigs.fasta" not in os.listdir("contigs_mapped/"):    
        print(f" - On téléchage {SRR} ({sra_list[SRR]})")
        os.system(f"fasterq-dump -e 16 --split-files --force --outdir fastq {SRR}")
        print(f" - Nettoyage {SRR} ({sra_list[SRR]})")
        subprocess.run(
            [
                "python3",
                "preprocess_reads.py",
                "-1",
                f"fastq/{SRR}_1.fastq",
                "-2",
                f"fastq/{SRR}_2.fastq",
                "-o",
                f"fastq/cleaned/{SRR}",
                "--threads",
                "16",
            ],
            check=True,
        )
        print(f" - On aligne {SRR} ({sra_list[SRR]})")
        os.system(
            f"bwa mem -t 16 data/H37Rv.fasta fastq/cleaned/{SRR}_1.fastq fastq/cleaned/{SRR}_2.fastq > alignments/{SRR}.sam"
        )
        os.system(f"samtools view -bS alignments/{SRR}.sam | samtools sort -o alignments/{SRR}_sorted.bam")
        os.system(f"rm -f alignments/{SRR}.sam")
        os.system(f"samtools index alignments/{SRR}_sorted.bam")
        # Extraction en BAM des reads non mappés (les deux mates non mappés)
        os.system(f"samtools view -b -f 12 -F 256 alignments/{SRR}_sorted.bam > unmapped/{SRR}_unmapped.bam")
        # Extraction en fastq directement avec samtools
        os.system(f"samtools fastq -1 unmapped/{SRR}_unmapped_1.fastq -2 unmapped/{SRR}_unmapped_2.fastq -0 /dev/null -s /dev/null -n unmapped/{SRR}_unmapped.bam")
        # Assemblage unmapped
        outdir = f"assemblies_unmapped/{SRR}"
        contigs = f"contigs_unmapped/{SRR}_unmapped_contigs.fasta"
        k_values = choose_k_values(f"unmapped/{SRR}_unmapped_1.fastq")
        os.system(
            f"spades.py --careful -t 16 --cov-cutoff auto -k {k_values} -1 unmapped/{SRR}_unmapped_1.fastq -2 unmapped/{SRR}_unmapped_2.fastq -o {outdir}"
        )
        assembled = f"{outdir}/contigs.fasta"
        if os.path.exists(assembled):
            os.rename(assembled, contigs)
        os.system(f"rm -f fastq/{SRR}*fastq")
        stop = True
    elif f"{SRR}.sam" in os.listdir("alignments/"):    
        os.system(f"rm -f alignments/{SRR}.sam")

    sorted_bam = f"alignments/{SRR}_sorted.bam"
    
    # 1) extraction des reads mappés
    mapped_bam = f"mapped/{SRR}_mapped.bam"
    if not os.path.exists(mapped_bam):
        os.system(f"samtools view -b -F 4 {sorted_bam} > {mapped_bam}")
    
    # 2) conversion en FASTQ
    mapped_fastq1 = f"mapped/{SRR}_mapped_1.fastq"
    mapped_fastq2 = f"mapped/{SRR}_mapped_2.fastq"
    if not (os.path.exists(mapped_fastq1) and os.path.exists(mapped_fastq2)):
        os.system(
            f"samtools fastq -1 {mapped_fastq1} -2 {mapped_fastq2} "
            f"-0 /dev/null -s /dev/null -n {mapped_bam}"
        )
    
    # 3) assemblage SPAdes
    outdir_m = f"assemblies_mapped/{SRR}"
    contigs_m = f"contigs_mapped/{SRR}_mapped_contigs.fasta"
    if not os.path.exists(contigs_m):
        k_values = choose_k_values(mapped_fastq1)
        os.system(
            f"spades.py -k {k_values} -1 {mapped_fastq1} -2 {mapped_fastq2} -o {outdir_m}"
        )
        assembled_m = f"{outdir_m}/contigs.fasta"
        if os.path.exists(assembled_m):
            os.rename(assembled_m, contigs_m)


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
            print(debut)
            print(fin)
            print(nom)
            assert debut < fin
            with open(f'data/RD/{nom}.fasta', 'w') as f:
                f.write(f'>{nom}\n')
                f.write(H37Rv[debut-1:fin-1])    
            concat(f"data/RD/{nom}.fasta")
        

    os.system(f"makeblastdb -in data/all_contigs.fasta -dbtype nucl -out mydb")
    if stop:
        exit()    

# blastn -query data/sequences/TbD1.fasta -db mydb

