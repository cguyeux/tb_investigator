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

La base BLAST finale est écrite dans le répertoire `bdd/` sous le nom
`mydb` pour éviter de créer des fichiers à la racine du projet.

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
correspondante. Lorsque plusieurs séquences sont fournies via `--filename`,
les pourcentages sont la moyenne calculée sur l'ensemble de ces entrées :

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

## Recherche dans la base interne `mydb`

Une base BLAST de séquences spécifiques de *M. tuberculosis* peut être placée
dans `bdd/mydb`. Le script la consulte automatiquement pour la séquence
analysée et pour chaque ORF détecté :

```bash
./analyse_seq.py genome.fasta --tb-db bdd/mydb
```

Les identifiants des séquences correspondantes sont affichés si la similarité
est d'au moins 95 % sur la moitié de la longueur de la requête.

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

Pour une annotation fonctionnelle plus complète, activez `--eggnog` pour lancer
`eggnog-mapper` sur les protéines prédites. Utilisez `--eggnog-data` pour
spécifier le répertoire de données et `--eggnog-cpu` pour ajuster le nombre de
threads.
Vous pouvez également fournir `--kegg-db` ou `--go-db` pour réaliser une
annotation par BLASTp contre une base KEGG Orthology ou Gene Ontology.

## Analyse de sélection positive/négative

Ajoutez l'option `--dnds` pour calculer un ratio dN/dS entre chaque ORF
prédit et son meilleur orthologue trouvé dans la base spécifiée par
`--dnds-db` (ou `--tb-db` si elle n'est pas fournie). Ce test rapide indique
si la séquence semble conservée (`dN/dS < 1`) ou soumise à une sélection
positive (`dN/dS > 1`).

## Localiser une perte dans H37Rv

En fournissant `--h37rv-db` (base BLAST du génome H37Rv) et `--h37rv-gff`,
`analyse_seq.py` indique la position de la meilleure correspondance et les gènes
présents dans une fenêtre définie par `--context-window` autour de cette région.

### Exemple complet

```bash
python analyse_seq.py \
  --filename ../L8_investigations/contigs_ERR12115321_filtered.fasta --seqname NODE_3_length_3231_cov_125.486367 \
  --filename ../L8_investigations/contigs_SRR10828835_filtered.fasta --seqname NODE_2_length_3097_cov_99.353311 \
  --filename ../L8_investigations/contigs_SRR1173284_filtered.fasta --seqname NODE_1_length_5572_cov_25.084622 \
  --filename ../L10_investigations/contigs_ERR2516384_filtered.fasta --seqname NODE_3_length_3096_cov_16.757204 \
  --filename ../L10_investigations/contigs_ERR2707158_filtered.fasta --seqname NODE_3_length_3107_cov_57.167987 \
  --lineage-db-dir bdd/ \
  --plasmid-db bdd/plsdb \
  --phage-db bdd/phagedb \
  --isescan ISEScan/isescan.py \
  --list-cds full \
  --prodigal-mode meta \
  --orf-db bdd/myco_proteins \
  --orf-db bdd/isfinder_prot \
  --hmmer \
  --pfam-db Pfam-A.hmm \
  --orf-keyword none \
  --orf-detailed \
  --eggnog \
  --eggnog-data eggnog-mapper/data \
  --eggnog-cpu 16 \
  --h37rv-db bdd/H37Rv.fasta \
  --h37rv-gff data/sequences/CDS/Mycobacterium_tuberculosis_H37Rv_gff_v5.gff \
  --kegg-db bdd/kegg_prot \
  --go-db bdd/go_prot \
  --context-window 8000 \
  --trnascan /usr/bin/trnascan-1.4 \
  --check-circular
```
