#!/usr/bin/env python3
import argparse
import os
import shutil
import subprocess
import sys
from collections import Counter
from typing import Dict, List, Any, Tuple
from ete3 import Tree


# Newick tree describing the M. tuberculosis complex lineages
LINEAGE_NEWICK = (
    "(M.riyadhense, (M.shinjukuense, (M.lacus, (M.decipiens,(Canettii, (L8, ((L1, (L7, "
    "((L4.1, (L4.2, (((L4.4, L4.13), (L4.17, (L4.3, L4.18))), (L4.14, (L4.5, "
    "((L4.6.1, L4.6.2), (L4.11, (L4.12, (L4.16, (L4.15, (L4.7, ((L4.9, L4.9H37Rv), L4.8)))))))))))), "
    "(L3, (L2.1proto, L2.2))))), (L5, (((Pinipedii, Microti), (OrygisLa3, (BovisLa1, (CapraeLa2, La4)))), "
    "((L10, ((L6.1, (L6.2, L6.3)), L9)), (Chimpanze, (Mungi, (Dassie, Suricattae)))))))))))));"
)


def run_prodigal(fasta: str, prefix: str, mode: str = "meta") -> tuple[str, str]:
    """Run prodigal and return paths to GFF and protein FASTA."""
    gff = f"{prefix}.gff"
    faa = f"{prefix}.faa"
    subprocess.run(
        [
            "prodigal",
            "-i",
            fasta,
            "-p",
            mode,
            "-a",
            faa,
            "-f",
            "gff",
            "-o",
            gff,
            "-q",
        ],
        check=True,
    )
    return gff, faa


def parse_gff(gff: str) -> List[Dict[str, str]]:
    """Parse Prodigal GFF and return CDS info."""
    cds = []
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            seqid, _source, type_, start, end, _score, strand, _phase, attrs = parts
            if type_ != "CDS":
                continue
            attr_dict: Dict[str, str] = {}
            for item in attrs.split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    attr_dict[key] = value
            cds.append(
                {
                    "seqid": seqid,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "id": attr_dict.get("ID", ""),
                }
            )
    return cds


def summarize_cds(cds: List[Dict[str, str]]) -> None:
    """Print a short summary of CDS."""
    print(f"{len(cds)} CDS predicted")
    strand_counts: Dict[str, int] = {}
    for c in cds:
        strand_counts[c["strand"]] = strand_counts.get(c["strand"], 0) + 1
    for strand, count in strand_counts.items():
        print(f"{count} CDS on strand {strand}")


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


def count_fasta_seqs(fasta: str) -> int:
    """Retourne le nombre d'entrées dans un fichier FASTA."""
    count = 0
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                count += 1
    return count


def parse_fasta_sequences(fasta: str) -> Dict[str, str]:
    """Parse un fichier FASTA et retourne un dictionnaire id->sequence."""
    sequences: Dict[str, str] = {}
    current_id = None
    seq_lines: List[str] = []
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(seq_lines)
                current_id = line[1:].strip().split()[0]
                seq_lines = []
            else:
                seq_lines.append(line.strip())
    if current_id is not None:
        sequences[current_id] = "".join(seq_lines)
    return sequences


def parse_prodigal_faa(faa: str) -> Dict[str, Dict[str, Any]]:
    """Parse les sorties FASTA de Prodigal et retourne info par ORF."""
    info: Dict[str, Dict[str, Any]] = {}
    current_id: str | None = None
    seq_lines: List[str] = []
    with open(faa) as fh:
        for line in fh:
            if line.startswith(">"):
                if current_id is not None:
                    info.setdefault(current_id, {})["sequence"] = "".join(seq_lines)
                header = line[1:].strip()
                parts = header.split("#")
                current_id = parts[0].strip().split()[0]
                if len(parts) >= 4:
                    try:
                        start = int(parts[1].strip())
                        end = int(parts[2].strip())
                        strand_raw = parts[3].strip()
                        strand = "-" if strand_raw.startswith("-") else "+"
                        info[current_id] = {
                            "start": start,
                            "end": end,
                            "strand": strand,
                        }
                    except ValueError:
                        info[current_id] = {}
                else:
                    info[current_id] = {}
                seq_lines = []
            else:
                seq_lines.append(line.strip())
    if current_id is not None:
        info.setdefault(current_id, {})["sequence"] = "".join(seq_lines)
    return info


def reverse_complement(seq: str) -> str:
    """Retourne le complément inverse d'une séquence ADN."""
    trans = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(trans)[::-1]


def add_nuc_sequences(orf_info: Dict[str, Dict[str, Any]], dna_seq: str) -> None:
    """Ajoute la séquence nucléotidique à chaque ORF."""
    for info in orf_info.values():
        start = info.get("start")
        end = info.get("end")
        if start is None or end is None:
            continue
        sub = dna_seq[start - 1 : end]
        if info.get("strand") == "-":
            sub = reverse_complement(sub)
        info["nuc_sequence"] = sub


def parse_gff_dict(gff: str) -> Dict[str, Dict[str, Any]]:
    """Parse un GFF Prodigal et retourne un dictionnaire par ID."""
    info: Dict[str, Dict[str, Any]] = {}
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            seqid, _source, type_, start, end, _score, strand, _phase, attrs = parts
            if type_ != "CDS":
                continue
            attr_dict: Dict[str, str] = {}
            for item in attrs.split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    attr_dict[key] = value
            orf_id = attr_dict.get("ID", "")
            info[orf_id] = {
                "seqid": seqid,
                "start": int(start),
                "end": int(end),
                "strand": strand,
            }
    return info


def first_fasta_header(fasta: str) -> str:
    """Retourne le premier identifiant d'un fichier FASTA."""
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                return line[1:].strip().split()[0]
    return ""


def extract_sequence(
    filename: str, seqname: str | None, tmpdir: str, prefix: str = "selected"
) -> str:
    """Extrait une séquence d'un multi-FASTA dans ``tmpdir``.

    Le fichier FASTA resultants sera ``prefix``.fasta dans ``tmpdir`` et son
    chemin est renvoyé.
    """
    os.makedirs(tmpdir, exist_ok=True)
    out_fasta = os.path.join(tmpdir, f"{prefix}.fasta")

    current_id = None
    wanted = seqname
    seq_lines: list[str] = []
    found = False
    with open(filename) as fh:
        for line in fh:
            if line.startswith(">"):
                header = line[1:].strip()
                ident = header.split()[0]
                if found:
                    break
                if wanted is None or ident == wanted:
                    current_id = ident
                    seq_lines = []
                    found = True
                continue
            if found:
                seq_lines.append(line.strip())

    if not found:
        raise ValueError(f"Sequence '{seqname}' not found in {filename}")

    with open(out_fasta, "w") as out:
        out.write(f">{current_id}\n")
        out.write("\n".join(seq_lines) + "\n")
    return out_fasta


def consensus_from_alignment(aln_fasta: str) -> str:
    """Construit une séquence consensus à partir d'un alignement FASTA."""
    seqs = parse_fasta_sequences(aln_fasta)
    if not seqs:
        return ""
    sequences = list(seqs.values())
    length = len(sequences[0])
    cons = []
    for i in range(length):
        column = [s[i] for s in sequences if i < len(s)]
        column = [c for c in column if c != "-"]
        if not column:
            cons.append("N")
        else:
            counts = Counter(column)
            base = max(counts.items(), key=lambda x: x[1])[0]
            cons.append(base)
    return "".join(cons)


def align_sequences(fastas: List[str], tmpdir: str) -> str:
    """Aligne plusieurs séquences avec MAFFT et renvoie le fichier consensus."""
    if len(fastas) == 1:
        return fastas[0]

    os.makedirs(tmpdir, exist_ok=True)
    lengths = [compute_total_length(f) for f in fastas]
    idx_ref = lengths.index(max(lengths))
    ref = fastas[idx_ref]
    fragments = [f for i, f in enumerate(fastas) if i != idx_ref]

    frag_path = os.path.join(tmpdir, "fragments.fasta")
    with open(frag_path, "w") as out:
        for path in fragments:
            with open(path) as fh:
                out.write(fh.read())

    ref_aln = os.path.join(tmpdir, "ref_aligned.fasta")
    with open(ref_aln, "w") as out:
        try:
            subprocess.run(["mafft", ref], check=True, text=True, stdout=out)
        except FileNotFoundError as exc:
            raise RuntimeError("mafft not found") from exc

    final_aln = os.path.join(tmpdir, "final_alignment.fasta")
    with open(final_aln, "w") as out:
        try:
            subprocess.run(
                ["mafft", "--addfragments", frag_path, ref_aln],
                check=True,
                text=True,
                stdout=out,
            )
        except FileNotFoundError as exc:
            raise RuntimeError("mafft not found") from exc

    print_header("Alignement multiple (MAFFT)")
    with open(final_aln) as fh:
        print(fh.read())

    consensus_seq = consensus_from_alignment(final_aln)
    consensus_path = os.path.join(tmpdir, "consensus.fasta")
    with open(consensus_path, "w") as out:
        out.write(">consensus\n")
        out.write(consensus_seq + "\n")
    return consensus_path


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


def blast_significant_hits(
    fasta: str,
    db: str,
    evalue: float = 1e-5,
    min_pident: float = 95.0,
    min_cov: float = 0.5,
) -> Dict[str, List[Dict[str, float]]]:
    """Return detailed BLASTn hits grouped by query.

    Each hit dictionary contains ``sseqid`` (subject id), ``pident`` (percent
    identity) as well as query and subject coverage percentages (``qcov`` and
    ``scov``).

    A hit is kept if ``pident`` is at least ``min_pident`` and the alignment
    length covers at least ``min_cov`` of the query sequence.
    """

    cmd = [
        "blastn",
        "-query",
        fasta,
        "-db",
        db,
        "-evalue",
        str(evalue),
        "-outfmt",
        "6 qseqid sseqid pident length qlen slen",
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except FileNotFoundError as exc:
        raise RuntimeError(
            "blastn not found. Install BLAST+ to use this option."
        ) from exc

    hits: Dict[str, List[Dict[str, float]]] = {}
    for line in result.stdout.strip().splitlines():
        parts = line.split("\t")
        if len(parts) < 6:
            continue
        qid, sid, pid, length, qlen, slen = parts[:6]
        try:
            pid_f = float(pid)
            length_i = int(length)
            qlen_f = float(qlen)
            slen_f = float(slen)
        except ValueError:
            continue
        qcov = length_i / qlen_f
        scov = length_i / slen_f if slen_f else 0.0
        if pid_f >= min_pident and qcov >= min_cov:
            hits.setdefault(qid, []).append(
                {
                    "sseqid": sid,
                    "pident": pid_f,
                    "qcov": qcov * 100,
                    "scov": scov * 100,
                }
            )
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
    faa: str, db: str, evalue: float = 1e-5, keyword: str | None = "transposase"
) -> List[Dict[str, str]]:
    """BLASTp des protéines contre *db* et retourne les champs parsés.

    Si ``keyword`` est fourni, seuls les hits contenant ce mot-clé sont
    retournés. Avec ``keyword`` à ``None`` ou ``"none"`` tous les résultats sont
    conservés. Chaque élément retourné est un dictionnaire contenant les clés
    ``qseqid``, ``sseqid``, ``pident``, ``length``, ``bitscore``, ``evalue`` et
    ``stitle``.
    """
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
    hits: List[Dict[str, str]] = []
    for line in result.stdout.strip().splitlines():
        if not line:
            continue
        if keyword is None or keyword == "" or keyword.lower() == "none" or keyword.lower() in line.lower():
            parts = line.split("\t")
            if len(parts) >= 7:
                hit = {
                    "qseqid": parts[0],
                    "sseqid": parts[1],
                    "pident": parts[2],
                    "length": parts[3],
                    "bitscore": parts[4],
                    "evalue": parts[5],
                    "stitle": parts[6],
                }
                hits.append(hit)
    return hits


def summarize_orf_hits(hits: List[Dict[str, str]]) -> List[str]:
    """Résumé simple des hits BLASTp par ORF."""
    best: Dict[str, Dict[str, str]] = {}
    for h in hits:
        if h["qseqid"] not in best:
            best[h["qseqid"]] = h
    summary = []
    for qseqid, info in best.items():
        desc = info.get("stitle", "")
        summary.append(
            f"{qseqid}: {desc} (identité {info['pident']}%, longueur {info['length']})"
        )
    return summary


def run_hmmer(
    faa: str, pfam_db: str, keyword: str | None = "transposase", evalue: float = 1e-5
) -> List[str]:
    """Recherche de domaines HMM PFAM dans les protéines.

    Si ``keyword`` est fourni, seuls les domaines contenant ce mot-clé sont
    rapportés. Avec ``keyword`` à ``None`` ou ``"none"`` tous les résultats sont
    retournés. Le paramètre ``evalue`` définit le seuil ``--domE`` passé à
    ``hmmscan``.
    """
    # pfam_db doit être une base hmmpressée
    cmd = [
        "hmmscan",
        "--domtblout",
        "hmmer.tbl",
        "--domE",
        str(evalue),
        pfam_db,
        faa,
    ]
    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
    except FileNotFoundError as exc:
        raise RuntimeError(
            "hmmscan not found. Install HMMER to use this option."
        ) from exc
    hits = []
    with open("hmmer.tbl") as tbl:
        for line in tbl:
            if line.startswith("#"):
                continue
            if keyword is None or keyword == "" or keyword.lower() == "none":
                hits.append(line.strip())
            elif keyword.lower() in line.lower():
                hits.append(line.strip())
    return hits


def parse_emapper_annotations(path: str) -> Dict[str, Dict[str, str]]:
    """Parse an eggNOG-mapper annotation file."""
    annotations: Dict[str, Dict[str, str]] = {}
    header: List[str] | None = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                if line.startswith("#query"):
                    header = line[1:].strip().split("\t")
                continue
            if header is None:
                continue
            parts = line.rstrip("\n").split("\t")
            row = {k: parts[i] if i < len(parts) else "" for i, k in enumerate(header)}
            qid = row.get("query") or row.get("query_name") or row.get("#query")
            if not qid:
                continue
            annotations[qid] = {
                "eggnog_preferred": row.get("Preferred_name", ""),
                "eggnog_description": row.get("Description", ""),
            }
    return annotations


def run_eggnog_mapper(
    faa: str, outdir: str, cpu: int = 1, data_dir: str | None = None
) -> Dict[str, Dict[str, str]]:
    """Run eggNOG-mapper on ``faa`` and return parsed annotations."""
    os.makedirs(outdir, exist_ok=True)
    prefix = "eggnog"
    cmd = [
        "eggnog-mapper/emapper.py",
        "-i",
        faa,
        "-o",
        prefix,
        "--override",
        "--output_dir",
        outdir,
        "--cpu",
        str(cpu),
    ]
    if data_dir:
        cmd.extend(["--data_dir", data_dir])
    try:
        subprocess.run(cmd, check=True)
    except FileNotFoundError as exc:
        raise RuntimeError(
            "emapper.py not found. Install eggnog-mapper to use this option."
        ) from exc
    ann_path = os.path.join(outdir, f"{prefix}.emapper.annotations")
    return parse_emapper_annotations(ann_path)


def blast_first_hit(fasta: str, db: str, evalue: float = 1e-5) -> Dict[str, Any] | None:
    """Return first BLASTn hit with coordinates against ``db``."""
    cmd = [
        "blastn",
        "-query",
        fasta,
        "-db",
        db,
        "-max_target_seqs",
        "1",
        "-outfmt",
        "6 sseqid sstart send qstart qend pident qcovs sstrand",
        "-evalue",
        str(evalue),
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except FileNotFoundError as exc:
        raise RuntimeError(
            "blastn not found. Install BLAST+ to use this option."
        ) from exc
    line = result.stdout.strip().splitlines()
    if not line:
        return None
    parts = line[0].split("\t")
    if len(parts) < 8:
        return None
    sseqid, sstart, send, qstart, qend, pid, qcov, strand = parts[:8]
    try:
        return {
            "sseqid": sseqid,
            "sstart": int(sstart),
            "send": int(send),
            "qstart": int(qstart),
            "qend": int(qend),
            "pident": float(pid),
            "qcov": float(qcov),
            "strand": strand,
        }
    except ValueError:
        return None


def parse_general_gff(gff: str) -> List[Dict[str, Any]]:
    """Parse a GFF file and return a list of gene/CDS features."""
    features: List[Dict[str, Any]] = []
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            seqid, _source, type_, start, end, _score, strand, _phase, attrs = parts
            if type_ not in {"gene", "CDS"}:
                continue
            attr_dict: Dict[str, str] = {}
            for item in attrs.split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    attr_dict[key] = value
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue
            features.append(
                {
                    "seqid": seqid,
                    "start": start_i,
                    "end": end_i,
                    "strand": strand,
                    "attrs": attr_dict,
                }
            )
    return features


def features_in_window(features: List[Dict[str, Any]], seqid: str, start: int, end: int) -> List[Dict[str, Any]]:
    """Return features overlapping [start, end] on ``seqid``."""
    subset = []
    for feat in features:
        if feat["seqid"] != seqid:
            continue
        if feat["end"] < start or feat["start"] > end:
            continue
        subset.append(feat)
    subset.sort(key=lambda x: x["start"])
    return subset


def print_header(title: str) -> None:
    """Affiche un titre encadré de séparateurs."""
    bar = "=" * 60
    print(f"\n{bar}\n{title}\n{bar}")


def print_orf_details(
    orf_info: Dict[str, Dict[str, Any]],
    seqname: str,
    lineage_db_dir: str | None = None,
    evalue: float = 1e-5,
    tmpdir: str = "tmp",
) -> None:
    """Affiche les séquences ORF et leurs informations.

    Si ``lineage_db_dir`` est fourni, un arbre ASCII du complexe est affiché
    pour chaque ORF avec le pourcentage de couverture BLASTn dans chaque
    sous-lignée.
    """
    for oid, info in sorted(orf_info.items(), key=lambda x: x[1].get("start", 0)):
        print("-" * 60)
        strand = info.get("strand", "+")
        orientation = "brin complément" if strand == "-" else "brin direct"
        short_oid = oid
        if oid.startswith(seqname + "_"):
            short_oid = oid[len(seqname) + 1 :]
        header = f"{seqname}_{short_oid}_{info.get('start','?')}-{info.get('end','?')}_{orientation}"
        prot_seq = info.get("sequence", "")
        nuc_seq = info.get("nuc_sequence", "")
        print(f">{header}")
        if prot_seq:
            print(prot_seq)
        if nuc_seq:
            print(nuc_seq)
        details = []
        if lineage_db_dir and nuc_seq:
            os.makedirs(tmpdir, exist_ok=True)
            orf_fasta = os.path.join(tmpdir, f"{oid}.fasta")
            with open(orf_fasta, "w") as out:
                out.write(f">{oid}\n{nuc_seq}\n")
            tree = Tree(LINEAGE_NEWICK)
            coverages: Dict[str, float] = {}
            for leaf in tree.iter_leaves():
                db_prefix = os.path.join(lineage_db_dir, leaf.name)
                try:
                    pct = blast_coverage(orf_fasta, db_prefix, evalue)
                except RuntimeError as err:
                    print(err)
                    pct = 0.0
                coverages[leaf.name] = pct
            for leaf in tree.iter_leaves():
                pct = coverages.get(leaf.name, 0.0)
                leaf.name = f"{leaf.name} ({pct:.1f}%)"
            details.append("Sous-lignées:")
            details.append(tree.get_ascii(attributes=[]).rstrip())
        if "blast_hits" in info and info["blast_hits"]:
            details.append("BLASTp:")
            for h in info["blast_hits"]:
                desc = h.get("stitle", "")
                details.append(
                    f"  {h['sseqid']} identité {h['pident']}% longueur {h['length']} - {desc}"
                )
        if "hmmer_hits" in info and info["hmmer_hits"]:
            details.append("HMMER:")
            for line in info["hmmer_hits"]:
                details.append(f"  {line}")
        if "tb_hits" in info and info["tb_hits"]:
            details.append("mydb:")
            for h in info["tb_hits"]:
                note = (
                    " (unmapped)" if "unmapped" in h.get("sseqid", "") else ""
                )
                details.append(
                    "  "
                    + f"{h['sseqid']} identité {h['pident']:.1f}% qcov {h['qcov']:.1f}% scov {h['scov']:.1f}%{note}"
                )
        if "eggnog_description" in info or "eggnog_preferred" in info:
            desc = info.get("eggnog_description", "")
            pref = info.get("eggnog_preferred", "")
            joined = f"{pref} {desc}".strip()
            details.append(f"eggNOG: {joined}")
        for d in details:
            print(d)


def main():
    parser = argparse.ArgumentParser(
        description="Analyse GC et recherche d'origine (plasmide, phage, IS, transposon, annotation d'ORFs)"
    )
    parser.add_argument(
        "--filename",
        action="append",
        required=True,
        help="Fichier multi-FASTA contenant la séquence à analyser. Peut être spécifié plusieurs fois.",
    )
    parser.add_argument(
        "--seqname",
        action="append",
        help="Nom de la séquence à extraire juste après chaque --filename. Si omis, la première séquence sera utilisée.",
    )
    parser.add_argument(
        "--plasmid-db",
        help="Base BLAST de plasmides (PLSDB), située dans bdd/",
    )
    parser.add_argument(
        "--phage-db",
        help="Base BLAST de phages, située dans bdd/",
    )
    parser.add_argument(
        "--tb-db",
        default="bdd/mydb",
        help="Base BLASTn de séquences connues de M. tuberculosis",
    )
    parser.add_argument("--isescan", help="Chemin vers isescan.py")
    parser.add_argument(
        "--tmpdir", default="tmp", help="Répertoire pour les fichiers temporaires"
    )
    parser.add_argument(
        "--list-cds",
        choices=["none", "summary", "full"],
        default="none",
        help=(
            "Contrôle l'affichage des CDS prédits par Prodigal : "
            "'summary' pour un résumé, 'full' pour tout afficher, "
            "'none' pour ne rien afficher"
        ),
    )
    parser.add_argument(
        "--prodigal-prefix",
        default="prodigal",
        help="Préfixe des fichiers générés par Prodigal",
    )
    parser.add_argument(
        "--prodigal-mode",
        choices=["meta", "single"],
        default="meta",
        help="Mode de Prodigal (meta ou single)",
    )
    parser.add_argument(
        "--evalue",
        type=float,
        default=1e-5,
        help="Seuil d'e-value pour BLAST et HMMER",
    )
    parser.add_argument(
        "--orf-db",
        action="append",
        help=(
            "Base BLAST protéines (e.g. NR ou ISFinder_proteins). "
            "Cherchée dans bdd/ si aucun chemin n'est fourni. "
            "Peut être spécifiée plusieurs fois."
        ),
    )
    parser.add_argument(
        "--orf-detailed",
        action="store_true",
        help="Afficher les lignes BLAST/HMMER complètes pour chaque ORF",
    )
    parser.add_argument(
        "--hmmer", action="store_true", help="Faire aussi une recherche HMMER/PFAM"
    )
    parser.add_argument("--pfam-db", help="Base HMM PFAM (hmmscan)")
    parser.add_argument(
        "--orf-keyword",
        default="transposase",
        help=(
            "Mot-clé pour filtrer les hits BLAST/HMMER sur les ORFs. "
            'Utiliser "none" pour désactiver le filtrage.'
        ),
    )
    parser.add_argument(
        "--eggnog",
        action="store_true",
        help="Annoter les ORFs avec eggnog-mapper",
    )
    parser.add_argument(
        "--eggnog-data",
        help="Répertoire des données eggnog-mapper",
    )
    parser.add_argument(
        "--eggnog-cpu",
        type=int,
        default=1,
        help="Nombre de threads pour eggnog-mapper",
    )
    parser.add_argument(
        "--h37rv-db",
        help="Base BLASTn du génome H37Rv pour localiser la perte",
    )
    parser.add_argument(
        "--h37rv-gff",
        help="Annotations GFF de H37Rv pour décrire le contexte",
    )
    parser.add_argument(
        "--context-window",
        type=int,
        default=1000,
        help="Taille en bp de la fenêtre autour de la perte",
    )
    parser.add_argument(
        "--lineage-db-dir",
        help="Répertoire des bases BLAST par sous-lignée pour afficher l'arbre",
    )
    args = parser.parse_args()

    def add_bdd_prefix(db: str) -> str:
        """Return path with bdd/ prefix if no directory is provided."""
        if db and os.sep not in db:
            return os.path.join("bdd", db)
        return db

    if args.plasmid_db:
        args.plasmid_db = add_bdd_prefix(args.plasmid_db)
    if args.phage_db:
        args.phage_db = add_bdd_prefix(args.phage_db)
    if args.tb_db:
        args.tb_db = add_bdd_prefix(args.tb_db)
    if args.h37rv_db:
        args.h37rv_db = add_bdd_prefix(args.h37rv_db)
    if args.orf_db:
        args.orf_db = [add_bdd_prefix(db) for db in args.orf_db]

    print_header("Extraction de la séquence")

    # Reconstitue la liste (filename, seqname) en suivant l'ordre
    def parse_pairs(argv: List[str]) -> List[Tuple[str, str | None]]:
        pairs: List[Tuple[str, str | None]] = []
        i = 0
        while i < len(argv):
            if argv[i] == "--filename":
                if i + 1 >= len(argv):
                    parser.error("--filename doit être suivi d'un chemin")
                fname = argv[i + 1]
                i += 2
                sname: str | None = None
                if i < len(argv) and argv[i] == "--seqname":
                    if i + 1 >= len(argv):
                        parser.error("--seqname doit être suivi d'un identifiant")
                    sname = argv[i + 1]
                    i += 2
                pairs.append((fname, sname))
            else:
                i += 1
        return pairs

    pairs = parse_pairs(sys.argv[1:])

    selected_fastas: List[str] = []
    for idx, (fname, sname) in enumerate(pairs):
        try:
            fasta_path = extract_sequence(fname, sname, args.tmpdir, f"selected_{idx}")
        except Exception as err:
            parser.error(str(err))
        if sname is None:
            print(f"{fname} : pas de --seqname fourni, première séquence utilisée")
        else:
            print(f"{fname} : séquence {sname} extraite")
        print(f"Fichier temporaire : {fasta_path}")
        selected_fastas.append(fasta_path)

    args.fasta = align_sequences(selected_fastas, args.tmpdir)
    selected_seq_id = first_fasta_header(args.fasta)

    if args.lineage_db_dir:
        print_header("Recherche de la sous-lignée (BLASTn)")
        print(f"Bases dans {args.lineage_db_dir}")
        newick = "(M.riyadhense, (M.shinjukuense, (M.lacus, (M.decipiens,(Canettii, (L8, ((L1, (L7, ((L4.1, (L4.2, (((L4.4, L4.13), (L4.17, (L4.3, L4.18))), (L4.14, (L4.5, ((L4.6.1, L4.6.2), (L4.11, (L4.12, (L4.16, (L4.15, (L4.7, ((L4.9, L4.9H37Rv), L4.8)))))))))))), (L3, (L2.1proto, L2.2))))), (L5, (((Pinipedii, Microti), (OrygisLa3, (BovisLa1, (CapraeLa2, La4)))), ((L10, ((L6.1, (L6.2, L6.3)), L9)), (Chimpanze, (Mungi, (Dassie, Suricattae)))))))))))));"
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

    print_header("Calcul du contenu GC")
    print(f"Fichier analysé : {args.fasta}")
    gc = compute_gc(args.fasta)
    print(f"GC content: {gc:.2f}%")

    if args.plasmid_db:
        print_header("Recherche de plasmides (BLASTn)")
        print(
            f"Commande : blastn -query {args.fasta} -db {args.plasmid_db} -evalue {args.evalue}"
        )
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
        print_header("Recherche de phages (BLASTn)")
        print(
            f"Commande : blastn -query {args.fasta} -db {args.phage_db} -evalue {args.evalue}"
        )
        try:
            hits = blast_hits(args.fasta, args.phage_db, args.evalue)
        except RuntimeError as err:
            print(err)
        else:
            if hits:
                print("Possible phage origin for:", ", ".join(hits))
            else:
                print("No phage match found.")

    if args.tb_db:
        print_header("Recherche de séquences M. tuberculosis (BLASTn)")
        print(
            f"Commande : blastn -query {args.fasta} -db {args.tb_db} -evalue {args.evalue}"
        )
        if os.path.normpath(args.tb_db) == os.path.normpath("bdd/mydb"):
            print(
                "Base interne bdd/mydb contenant des contigs et regions specifiques de M. tuberculosis."
            )
        try:
            tb_hits_map = blast_significant_hits(
                args.fasta, args.tb_db, args.evalue
            )
            tb_hits = tb_hits_map.get(selected_seq_id, [])
        except RuntimeError as err:
            print(err)
        else:
            if tb_hits:
                for h in tb_hits:
                    note = (
                        " (match dans une zone ne mappant pas sur H37Rv)"
                        if "unmapped" in h["sseqid"]
                        else ""
                    )
                    print(
                        f"{h['sseqid']}: identite {h['pident']:.1f}% qcov {h['qcov']:.1f}% scov {h['scov']:.1f}%{note}"
                    )
            else:
                print("Aucune correspondance significative.")

    if args.h37rv_db:
        print_header("Position dans H37Rv (BLASTn)")
        try:
            hit = blast_first_hit(args.fasta, args.h37rv_db, args.evalue)
        except RuntimeError as err:
            print(err)
        else:
            if hit:
                start = min(hit["sstart"], hit["send"])
                end = max(hit["sstart"], hit["send"])
                print(
                    f"Correspondance principale : {hit['sseqid']}:{start}-{end} ({hit['strand']})"
                )
                if args.h37rv_gff:
                    try:
                        feats = parse_general_gff(args.h37rv_gff)
                    except FileNotFoundError:
                        print(f"Fichier {args.h37rv_gff} introuvable")
                    else:
                        win_start = max(1, start - args.context_window)
                        win_end = end + args.context_window
                        context = features_in_window(
                            feats, hit["sseqid"], win_start, win_end
                        )
                        if context:
                            print(
                                f"G\u00e8nes dans la fen\u00eatre {win_start}-{win_end} :"
                            )
                            for feat in context:
                                attrs = feat["attrs"]
                                name = (
                                    attrs.get("gene")
                                    or attrs.get("Name")
                                    or attrs.get("locus_tag", "")
                                )
                                product = attrs.get("product", "")
                                print(
                                    f"- {name} {feat['start']}-{feat['end']} {feat['strand']} {product}"
                                )
                        else:
                            print("Aucun g\u00e8ne dans la fen\u00eatre sp\u00e9cifi\u00e9e.")
            else:
                print("Pas de correspondance dans H37Rv.")

    if args.isescan:
        print_header("Recherche d'IS avec ISEScan")
        outdir = os.path.join(args.tmpdir, "isescan")
        print(
            f"Commande : {args.isescan} --seqfile {args.fasta} --output {outdir}"
        )
        try:
            found = run_isescan(args.fasta, args.isescan, outdir)
            msg = (
                "ISEScan detected IS elements."
                if found
                else "ISEScan found no IS elements."
            )
            print(msg)
        except RuntimeError as err:
            print(err)

    need_prodigal = args.list_cds != "none" or args.orf_db or (
        args.hmmer and args.pfam_db
    )

    orf_info: Dict[str, Dict[str, Any]] = {}

    if need_prodigal:
        prefix = os.path.join(args.tmpdir, args.prodigal_prefix)
        os.makedirs(os.path.dirname(prefix), exist_ok=True)
        if args.list_cds != "none":
            print_header("Prédiction des CDS avec Prodigal")
            print(
                f"Commande : prodigal -i {args.fasta} -p {args.prodigal_mode} -a {prefix}.faa -f gff -o {prefix}.gff -q"
            )
        try:
            gff, faa = run_prodigal(args.fasta, prefix, args.prodigal_mode)
            orf_dir = os.path.join(args.tmpdir, "orfs")
            os.makedirs(orf_dir, exist_ok=True)
            shutil.copy(faa, os.path.join(orf_dir, "orfs.faa"))
            if args.list_cds in ("summary", "full"):
                cds = parse_gff(gff)
                if args.list_cds == "full":
                    for c in cds:
                        print(
                            "\t".join([
                                c["seqid"],
                                c["start"],
                                c["end"],
                                c["strand"],
                                c["id"],
                            ])
                        )
                else:
                    summarize_cds(cds)
                print(f"CDS annotations written to {gff}")
                print(f"Protein translations written to {faa}")
            orf_info = parse_prodigal_faa(faa)
            dna_seq = parse_fasta_sequences(args.fasta).get(selected_seq_id, "")
            if dna_seq:
                add_nuc_sequences(orf_info, dna_seq)
        except Exception as err:
            print(f"Prodigal failed: {err}")


    # Recherche optionnelle d'ORFs et annotation par BLAST
    if args.orf_db or (args.hmmer and args.pfam_db):
        faa = os.path.join(args.tmpdir, "orfs", "orfs.faa")
        print_header("Recherche de protéines par ORF/BLASTp")
        orf_total = count_fasta_seqs(faa)
        print(f"{orf_total} ORFs prédits")
        if args.orf_db:
            for db in args.orf_db:
                print(f"\nRecherche dans la base {db} (BLASTp evalue {args.evalue}) :")
                try:
                    hits = blastx_hits(faa, db, args.evalue, args.orf_keyword)
                    if hits:
                        for h in hits:
                            oid = h["qseqid"]
                            orf_info.setdefault(oid, {}).setdefault("blast_hits", []).append(h)
                        summaries = summarize_orf_hits(hits)
                        print(f"{len(summaries)} ORFs avec hits :")
                        for s in summaries:
                            print(s)
                        if args.orf_detailed:
                            for h in hits:
                                print(
                                    "\t".join(
                                        [
                                            h["qseqid"],
                                            h["sseqid"],
                                            h["pident"],
                                            h["length"],
                                            h["bitscore"],
                                            h["evalue"],
                                            h["stitle"],
                                        ]
                                    )
                                )
                    else:
                        print("Aucun hit BLASTp trouvé.")
                except RuntimeError as err:
                    print(err)
        if args.hmmer and args.pfam_db:
            print_header("Recherche de domaines PFAM avec HMMER")
            print(
                f"Commande : hmmscan --domtblout hmmer.tbl --domE {args.evalue} {args.pfam_db} {faa}"
            )
            try:
                hits = run_hmmer(
                    faa,
                    args.pfam_db,
                    args.orf_keyword,
                    args.evalue,
                )
                if hits:
                    print(f"{len(hits)} domaines HMMER trouvés")
                    for line in hits:
                        oid = line.split()[0]
                        orf_info.setdefault(oid, {}).setdefault("hmmer_hits", []).append(line)
                    if args.orf_detailed:
                        for h in hits:
                            print(h)
                else:
                    print("Aucun domaine trouvé (HMMER/PFAM).")
            except RuntimeError as err:
                print(err)

        if args.eggnog:
            print_header("Annotation fonctionnelle avec eggnog-mapper")
            try:
                ann = run_eggnog_mapper(
                    faa,
                    os.path.join(args.tmpdir, "eggnog"),
                    args.eggnog_cpu,
                    args.eggnog_data,
                )
                if ann:
                    for oid, info in ann.items():
                        orf_info.setdefault(oid, {}).update(info)
                    print(f"{len(ann)} ORFs annotés par eggnog-mapper")
                else:
                    print("Aucune annotation eggnog-mapper trouvée.")
            except RuntimeError as err:
                print(err)

        if args.tb_db:
            print_header("Recherche des ORFs dans la base TB (BLASTn)")
            nuc_faa = os.path.join(args.tmpdir, "orfs", "orfs_nuc.fasta")
            with open(nuc_faa, "w") as out:
                for oid, info in orf_info.items():
                    seq = info.get("nuc_sequence")
                    if not seq:
                        continue
                    out.write(f">{oid}\n{seq}\n")
            try:
                hits_map = blast_significant_hits(nuc_faa, args.tb_db, args.evalue)
            except RuntimeError as err:
                print(err)
            else:
                if hits_map:
                    for oid, hits in hits_map.items():
                        if hits:
                            orf_info.setdefault(oid, {}).setdefault("tb_hits", []).extend(hits)
                    total = sum(1 for n in hits_map.values() if n)
                    if total:
                        print(f"{total} ORFs avec correspondance dans mydb")
                        for oid, hits in hits_map.items():
                            if hits:
                                for h in hits:
                                    note = (
                                        " (match dans une zone ne mappant pas sur H37Rv)"
                                        if "unmapped" in h["sseqid"]
                                        else ""
                                    )
                                    print(
                                        f"{oid}: {h['sseqid']} identite {h['pident']:.1f}% qcov {h['qcov']:.1f}% scov {h['scov']:.1f}%{note}"
                                    )
                    else:
                        print("Aucun ORF n'a de correspondance significative.")
                else:
                    print("Aucun ORF n'a de correspondance significative.")

    if need_prodigal:
        print_header("Détails des ORFs")
        print_orf_details(
            orf_info,
            selected_seq_id,
            args.lineage_db_dir,
            args.evalue,
            args.tmpdir,
        )


if __name__ == "__main__":
    main()
