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

## Assemblage avec SPAdes

Le script `make_blastdb.py` assemble les lectures non mappées et mappées avec
**SPAdes**. Les valeurs passées à l'option `-k` sont déterminées
automatiquement à partir de la longueur des lectures du fichier FASTQ. Les
tailles de k-mers 21, 33, 55, 77, 99 et 127 sont filtrées pour ne conserver que
celles strictement inférieures à la longueur des lectures.

Par exemple, pour des lectures de 150 nt, tous les k-mers précédents sont
utilisés, alors que pour des lectures de 100 nt seuls `21,33,55,77,99` seront
passés à SPAdes.

Lors de la construction des bases de données de régions de délétions (RD),
les en-têtes FASTA sont nettoyés pour supprimer les caractères spéciaux.
Ceci évite les erreurs `makeblastdb` quand des symboles non ASCII comme `Δ`
apparaissent dans les noms de séquences.

## Calcul du taux de GC

Le script `analyse_seq.py` calcule le pourcentage global de GC d'un fichier FASTA.

### Usage

```bash
./analyse_seq.py sequences.fasta
```

Le pourcentage de GC est affiché sur la sortie standard avec deux décimales.

### Affichage de l'arbre des sous-lignées

Si vous disposez de bases BLAST par sous-lignée (dans un répertoire `bdd/` par
exemple), le script peut afficher un arbre phylogénétique ASCII indiquant pour
chaque feuille le pourcentage de la séquence retrouvée dans la base
correspondante :

```bash
./analyse_seq.py sequences.fasta --lineage-db-dir bdd
```

## Vérification de l'origine des contigs

Le même script peut rechercher si des contigs correspondent à des plasmides,
des phages ou des éléments transposables.

### Bases de données et outils requis

- **PLSDB** pour les plasmides. Téléchargez les séquences puis créez une base
  BLAST :

  ```bash
  makeblastdb -in plsdb.fasta -dbtype nucl -out plsdb_db
  ```

- **PHASTER** permet l'analyse de séquences phagiques en ligne sur
  <https://phaster.ca/>. Si vous disposez d'une base locale de phages,
  indiquez-la avec `--phage-db`.

- **ISEScan** et **TransposonPSI** sont utilisés pour détecter les éléments
  d'insertion. Installez-les séparément puis fournissez leur chemin au script.

### Exemple d'utilisation

```bash
./analyse_seq.py contigs.fasta --plasmid-db plsdb_db --phage-db phages_db \
    --isescan /path/to/isescan.py --transposonpsi /path/to/TransposonPSI.pl
```

## Lister les CDS prédits par Prodigal

`analyse_seq.py` peut également lancer Prodigal pour prédire les CDS. Les
annotations GFF et les protéines traduites sont écrites à l'emplacement
spécifié par `--prodigal-prefix` (dans `--tmpdir` par défaut). Les protéines
prédites sont toujours copiées dans `tmp/orfs/orfs.faa` pour une éventuelle
recherche BLAST ou HMMER.

### Usage

```bash
./analyse_seq.py genome.fasta --list-cds summary --prodigal-prefix resultat/prodigal
```

Les fichiers `resultat/prodigal.gff` et `resultat/prodigal.faa` contiendront
respectivement les coordonnées des CDS et leurs séquences protéiques. Utilisez
`--list-cds full` pour afficher toutes les entrées ou `--list-cds none` pour ne
pas afficher de résumé.

## Recherche de protéines par BLASTX

Spécifiez une ou plusieurs bases avec `--orf-db` pour interroger les protéines
prédites. Par défaut, seuls les hits contenant le mot clé "transposase" sont
résumés. Utilisez `--orf-detailed` pour obtenir la sortie BLAST complète ou
changez le mot-clé avec `--orf-keyword` ("none" pour désactiver le filtrage).

### Exemple

```bash
./analyse_seq.py genome.fasta --orf-db myco_proteins --orf-db isfinder_prot
```

Une recherche de domaines PFAM peut être ajoutée en indiquant `--hmmer` et la
base `--pfam-db`. Le seuil de recherche est contrôlé par `--evalue`, utilisé
aussi bien pour BLAST que pour HMMER. L'option `--orf-detailed` affiche alors
les lignes complètes des résultats.
