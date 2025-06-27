#!/usr/bin/env python3
import argparse
import os
import subprocess
from typing import Dict, List

from ete3 import Tree


def compute_gc(fasta: str) -> float:
    """Calcule le % GC d'un fichier fasta."""
    total_len = 0
    gc_count = 0
    seq = []
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                if seq:
                    sequence = "".join(seq).upper()
                    gc_count += sequence.count("G") + sequence.count("C")
                    total_len += len(sequence)
                    seq = []
            else:
                seq.append(line.strip())
        if seq:
            sequence = "".join(seq).upper()
            gc_count += sequence.count("G") + sequence.count("C")
            total_len += len(sequence)
    if total_len == 0:
        return 0.0
    return gc_count / total_len * 100


def compute_total_length(fasta: str) -> int:
    """Retourne la longueur totale (en nt) des séquences d'un FASTA."""
    total = 0
    seq = []
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                if seq:
                    total += len("".join(seq))
                    seq = []
            else:
                seq.append(line.strip())
        if seq:
            total += len("".join(seq))
    return total


def blast_coverage(fasta: str, db: str, evalue: float = 1e-5) -> float:
    """Pourcentage de la séquence avec un hit BLAST dans ``db``."""
    total_len = compute_total_length(fasta)
    if total_len == 0:
        return 0.0
    cmd = [
        "blastn",
        "-query",
        fasta,
        "-db",
        db,
        "-max_target_seqs",
        "1",
        "-max_hsps",
        "1",
        "-evalue",
        str(evalue),
        "-outfmt",
        "6 qseqid qlen length",
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except FileNotFoundError as exc:
        raise RuntimeError(
            "blastn not found. Install BLAST+ to use this option."
        ) from exc
    covered = 0
    seen: set[str] = set()
    for line in result.stdout.strip().splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        qid, qlen, length = parts[:3]
        if qid not in seen:
            covered += int(length)
            seen.add(qid)
    return covered / total_len * 100


def blast_hits(fasta: str, db: str, evalue: float = 1e-5) -> List[str]:
    """Retourne les ID des queries avec hit contre *db* via BLASTn."""
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
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except FileNotFoundError as exc:
        raise RuntimeError(
            "blastn not found. Install BLAST+ to use this option."
        ) from exc
    hits = []
    for line in result.stdout.strip().splitlines():
        parts = line.split("\t")
        if parts:
            hits.append(parts[0])
    return hits


def run_isescan(fasta: str, isescan: str, outdir: str) -> bool:
    """Lance ISEScan sur *fasta*. Retourne True si IS trouvés."""
    os.makedirs(outdir, exist_ok=True)
    try:
        result = subprocess.run(
            [isescan, "--seqfile", fasta, "--output", outdir],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"{isescan} not found. Install ISEScan to use this option."
        ) from exc
    if result.returncode != 0:
        raise RuntimeError(
            f"ISEScan failed with status {result.returncode}:\n{result.stdout}"
        )
    gff = os.path.join(outdir, "prediction.gff3")
    return os.path.exists(gff)


def predict_orfs(fasta: str, outdir: str) -> str:
    """Prédis les ORFs avec prodigal (mode metagenomic, robustes pour contigs)."""
    os.makedirs(outdir, exist_ok=True)
    orfs_faa = os.path.join(outdir, "orfs.faa")
    try:
        subprocess.run(
            ["prodigal", "-i", fasta, "-a", orfs_faa, "-p", "meta", "-q"], check=True
        )
    except Exception as e:
        raise RuntimeError(f"Prodigal failed: {e}")
    return orfs_faa


def blastx_hits(
    faa: str, db: str, evalue: float = 1e-5, keyword="transposase"
) -> List[str]:
    """BLASTx des protéines contre *db* et recherche du mot-clé (e.g. transposase)."""
    cmd = [
        "blastp",  # blastx si tu pars du nucl, blastp si orfs.faa
        "-query",
        faa,
        "-db",
        db,
        "-evalue",
        str(evalue),
        "-outfmt",
        "6 qseqid sseqid pident length bitscore evalue stitle",
        "-max_target_seqs",
        "5",
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except FileNotFoundError as exc:
        raise RuntimeError(
            "blastp not found. Install BLAST+ to use this option."
        ) from exc
    hits = []
    for line in result.stdout.strip().splitlines():
        if keyword.lower() in line.lower():
            hits.append(line)
    return hits


def run_hmmer(faa: str, pfam_db: str, keyword="transposase") -> List[str]:
    """Recherche de domaines HMM PFAM transposase dans les protéines."""
    # pfam_db doit être une base hmmpressée
    cmd = ["hmmscan", "--tblout", "hmmer.tbl", pfam_db, faa]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except FileNotFoundError as exc:
        raise RuntimeError(
            "hmmscan not found. Install HMMER to use this option."
        ) from exc
    hits = []
    for line in result.stdout.strip().splitlines():
        if not line.startswith("#") and keyword.lower() in line.lower():
            hits.append(line)
    return hits


def main():
    parser = argparse.ArgumentParser(
        description="Analyse GC et recherche d'origine (plasmide, phage, IS, transposon, orf transposase)"
    )
    parser.add_argument("fasta", help="Fichier FASTA à analyser")
    parser.add_argument("--plasmid-db", help="Base BLAST de plasmides (PLSDB)")
    parser.add_argument("--phage-db", help="Base BLAST de phages")
    parser.add_argument("--isescan", help="Chemin vers isescan.py")
    parser.add_argument(
        "--tmpdir", default="tmp", help="Répertoire pour les fichiers temporaires"
    )
    parser.add_argument("--evalue", type=float, default=1e-5)
    parser.add_argument(
        "--orf-search",
        action="store_true",
        help="Prédire les ORFs et faire une recherche de transposase",
    )
    parser.add_argument(
        "--orf-db", help="Base BLAST protéines (e.g. NR ou ISFinder_proteins)"
    )
    parser.add_argument(
        "--hmmer", action="store_true", help="Faire aussi une recherche HMMER/PFAM"
    )
    parser.add_argument("--pfam-db", help="Base HMM PFAM (hmmscan)")
    parser.add_argument(
        "--lineage-db-dir",
        help="Répertoire des bases BLAST par sous-lignée pour afficher l'arbre",
    )
    args = parser.parse_args()

    if args.lineage_db_dir:
        newick = "(Canettii, (L8, ((L1, (L7, ((L4.1, (L4.2, (((L4.4, L4.13), (L4.17, (L4.3, L4.18))), (L4.14, (L4.5, ((L4.6.1, L4.6.2), (L4.11, (L4.12, (L4.16, (L4.15, (L4.7, ((L4.9, L4.9H37Rv), L4.8)))))))))))), (L3, (L2.1proto, L2.2))))), (L5, (((Pinipedii, Microti), (OrygisLa3, (BovisLa1, (CapraeLa2, La4)))), ((L10, ((L6.1, (L6.2, L6.3)), L9)), (Dassie, Suricattae)))))));"
        tree = Tree(newick)
        coverages: Dict[str, float] = {}
        for leaf in tree.iter_leaves():
            db_prefix = os.path.join(args.lineage_db_dir, leaf.name)
            try:
                pct = blast_coverage(args.fasta, db_prefix, args.evalue)
            except RuntimeError as err:
                print(err)
                pct = 0.0
            coverages[leaf.name] = pct
        for leaf in tree.iter_leaves():
            pct = coverages.get(leaf.name, 0.0)
            leaf.name = f"{leaf.name} ({pct:.1f}%)"
        print(tree.get_ascii(attributes=[]))

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
            found = run_isescan(
                args.fasta, args.isescan, os.path.join(args.tmpdir, "isescan")
            )
            msg = (
                "ISEScan detected IS elements."
                if found
                else "ISEScan found no IS elements."
            )
            print(msg)
        except RuntimeError as err:
            print(err)

    # Désactivation par défaut de TransposonPSI car obsolète/incompatible
    # if args.transposonpsi:
    #     print("ATTENTION : TransposonPSI nécessite BLAST legacy et est souvent incompatible BLAST+.")
    #     # Ajoute ici si besoin la gestion si BLAST legacy disponible

    # Ajout de la recherche ORF et transposase
    if args.orf_search:
        print("\nRecherche de transposases ou intégrases par annotation ORF/blastx :")
        faa = predict_orfs(args.fasta, os.path.join(args.tmpdir, "orfs"))
        if args.orf_db:
            try:
                hits = blastx_hits(faa, args.orf_db, args.evalue)
                if hits:
                    print("Hits transposase/integrase trouvés (blastp) :")
                    for h in hits:
                        print(h)
                else:
                    print("Aucun hit transposase/integrase (blastp).")
            except RuntimeError as err:
                print(err)
        if args.hmmer and args.pfam_db:
            try:
                hits = run_hmmer(faa, args.pfam_db)
                if hits:
                    print("Domaines transposase trouvés (HMMER) :")
                    for h in hits:
                        print(h)
                else:
                    print("Aucun domaine transposase (HMMER/PFAM).")
            except RuntimeError as err:
                print(err)


if __name__ == "__main__":
    main()
