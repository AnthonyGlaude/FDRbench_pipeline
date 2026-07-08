#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Étape 2 — Vérification I/L et préparation des fichiers finaux d'analyse.

Cette étape prend les sorties du cleanup basique :

    analysis_outdir/{cell}/01_basic_cleanup/
        canonical_final_merged_non_mutant.tsv
        canonical_final_merged_mutant_only.tsv
        nc_final_merged_non_mutant.tsv
        nc_final_merged_mutant_only.tsv

et produit :

    analysis_outdir/{cell}/02_IL_filtering/final_analysis_input/
        canonical.tsv
        canonical_variant.tsv
        nc_wt.tsv
        nc_variant.tsv

Règles
------
1) canonical.tsv
   - copie directe de canonical_final_merged_non_mutant.tsv
   - pas de filtre gene-mismatch
   - pas de filtre I/L

2) canonical_variant.tsv
   - source : canonical_final_merged_mutant_only.tsv
   - filtre gene-consistency
   - filtre I/L contre référence théorique non-mutante

3) nc_wt.tsv
   - source : nc_final_merged_non_mutant.tsv
   - filtre gene-consistency
   - filtre cousins I/L dans référence théorique non-mutante

4) nc_variant.tsv
   - source : nc_final_merged_mutant_only.tsv
   - filtre gene-consistency
   - filtre I/L contre référence théorique non-mutante

Conçu pour être appelé depuis main.py avec :

    from module.il_filtering import run_il_filtering
    il_outputs = run_il_filtering(config)
"""

from __future__ import annotations

from pathlib import Path
import ast
import re
import shutil
from typing import Optional

import pandas as pd


PEPTIDE_COL = "peptide"
PROTEIN_COL = "protein"

MISSING_GENE_VALUES = {
    "", "nan", "none", "na", "n/a", "null", "missing", "unknown", "-", "."
}


# ============================================================
# General helpers
# ============================================================

def normalize_il(seq: object) -> str:
    """Remplace I et L par J pour gérer l'ambiguïté isobarique I/L."""
    if pd.isna(seq):
        return ""
    return str(seq).strip().upper().replace("I", "J").replace("L", "J")


def clean_gene_value(x: object) -> Optional[str]:
    """Retourne un symbole/accession propre, ou None si valeur manquante."""
    if x is None:
        return None

    value = str(x).strip()
    if value.lower() in MISSING_GENE_VALUES:
        return None

    return value


def parse_protein_list(x: object) -> list[str]:
    """
    Parse la colonne protein de MSFragger/FDRBench.

    Elle ressemble souvent à :
        ['header1', 'header2']
    """
    if isinstance(x, list):
        return [str(v) for v in x if v is not None]

    if pd.isna(x):
        return []

    s = str(x).strip()
    if not s:
        return []

    try:
        value = ast.literal_eval(s)
        if isinstance(value, list):
            return [str(v) for v in value if v is not None]
        return [str(value)]
    except Exception:
        return [s]


def extract_gn_from_header(header: object) -> Optional[str]:
    """Extrait GN=... d'un header FASTA/protéine."""
    h = str(header)
    match = re.search(r"(?:^|\|)GN=([^|]+)", h)
    if not match:
        return None
    return clean_gene_value(match.group(1))


def extract_gene_accession_from_header(header: object) -> Optional[str]:
    """
    Extrait une accession ENSG si présente.

    Exemple :
        gene:ENSG00000277282.1 -> ENSG00000277282
    """
    h = str(header)
    match = re.search(r"gene:(ENSG[0-9]+)(?:\.[0-9]+)?", h)
    if not match:
        return None
    return clean_gene_value(match.group(1))


def has_gene_mismatch(protein_value: object) -> bool:
    """
    Retourne True seulement si un PSM mappe vers plusieurs vrais gènes.

    Les valeurs GN=nan, GN=None, GN=NA, etc. sont ignorées.
    """
    proteins = parse_protein_list(protein_value)

    if len(proteins) <= 1:
        return False

    gene_accessions = {
        g for g in (extract_gene_accession_from_header(p) for p in proteins)
        if g is not None
    }
    if len(gene_accessions) > 1:
        return True
    if len(gene_accessions) == 1:
        return False

    gene_symbols = {
        g for g in (extract_gn_from_header(p) for p in proteins)
        if g is not None
    }

    return len(gene_symbols) > 1


def apply_gene_consistency_filter(df: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """Retire les lignes dont la liste de protéines contient plusieurs vrais gènes."""
    if PROTEIN_COL not in df.columns:
        raise ValueError(
            f"Colonne '{PROTEIN_COL}' absente. Colonnes disponibles : {list(df.columns)}"
        )

    mismatch_mask = df[PROTEIN_COL].map(has_gene_mismatch)
    filtered = df.loc[~mismatch_mask].copy()

    return filtered, int(mismatch_mask.sum())


# ============================================================
# Référence théorique digérée
# ============================================================

def is_mutant_header(header: object) -> bool:
    """Détecte les entrées mutantes dans le header de la référence digérée."""
    h = str(header)
    h_upper = h.upper()
    return (
        "_mut_" in h
        or "|type=mut" in h
        or "type=mut" in h
        or "_MUT_" in h_upper
    )


def parse_kiwi_line(line: str) -> Optional[tuple[str, str]]:
    """
    Parse une ligne CSV KIWI sans header :
        peptide, protein_header, mass

    Le header peut contenir des virgules.
    On prend donc :
        peptide = avant la première virgule
        mass = après la dernière virgule
        protein_header = entre les deux
    """
    line = line.rstrip("\n")
    if not line.strip():
        return None

    left, _mass = line.rsplit(",", 1)
    peptide, protein_header = left.split(",", 1)

    return peptide.strip(), protein_header.strip()


def load_nonmutant_reference(digested_csv: Path) -> tuple[set[str], dict[str, set[str]]]:
    """
    Charge la référence théorique NON-MUTANTE depuis la digestion KIWI.

    Retourne :
        ref_il_set :
            ensemble des peptides non-mutants après collapse I/L

        ref_il_to_exact :
            dict peptide_IL -> peptides exacts observés
    """
    if not digested_csv.exists():
        raise FileNotFoundError(f"Fichier digéré introuvable : {digested_csv}")

    ref_il_set: set[str] = set()
    ref_il_to_exact: dict[str, set[str]] = {}

    n_lines = 0
    n_bad = 0
    n_mutant_skipped = 0
    n_nonmut_used = 0

    print(f"[LOAD NON-MUT REF] {digested_csv}")

    with open(digested_csv, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            n_lines += 1

            try:
                parsed = parse_kiwi_line(line)
                if parsed is None:
                    continue

                peptide, protein_header = parsed

                if is_mutant_header(protein_header):
                    n_mutant_skipped += 1
                    continue

                exact_peptide = str(peptide).strip().upper()
                pep_il = normalize_il(exact_peptide)

                if not pep_il:
                    continue

                ref_il_set.add(pep_il)
                ref_il_to_exact.setdefault(pep_il, set()).add(exact_peptide)
                n_nonmut_used += 1

            except Exception:
                n_bad += 1

    print(f"[REF] lines read             : {n_lines:,}")
    print(f"[REF] mutant headers skipped : {n_mutant_skipped:,}")
    print(f"[REF] non-mutant used        : {n_nonmut_used:,}")
    print(f"[REF] unique IL peptides     : {len(ref_il_set):,}")

    if n_bad:
        print(f"[WARNING] bad lines ignored  : {n_bad:,}")

    if not ref_il_set:
        raise RuntimeError(
            f"Référence non-mutante vide pour {digested_csv}. "
            "Vérifie le chemin ou le format du CSV digéré."
        )

    return ref_il_set, ref_il_to_exact


# ============================================================
# Filtres I/L
# ============================================================

def filter_variant_against_nonmut_ref(
    df: pd.DataFrame,
    ref_il_set: set[str],
) -> tuple[pd.DataFrame, int]:
    """
    Pour les fichiers variants/mutants :
    retire les peptides dont la version I/L-collapsed existe dans la référence non-mutante.
    """
    if PEPTIDE_COL not in df.columns:
        raise ValueError(
            f"Colonne '{PEPTIDE_COL}' absente. Colonnes disponibles : {list(df.columns)}"
        )

    peptide_il = df[PEPTIDE_COL].map(normalize_il)
    remove_mask = peptide_il.isin(ref_il_set)

    filtered = df.loc[~remove_mask].copy()

    return filtered, int(remove_mask.sum())


def has_different_il_cousin_in_ref(
    peptide: object,
    ref_il_to_exact: dict[str, set[str]],
) -> bool:
    """
    Pour nc_wt :
    retire seulement si la référence contient un peptide exact différent
    avec le même collapse I/L.
    """
    if pd.isna(peptide):
        return False

    exact = str(peptide).strip().upper()
    pep_il = normalize_il(exact)
    ref_exacts = ref_il_to_exact.get(pep_il, set())

    return any(ref_exact != exact for ref_exact in ref_exacts)


def filter_nonmutant_il_cousins(
    df: pd.DataFrame,
    ref_il_to_exact: dict[str, set[str]],
) -> tuple[pd.DataFrame, int]:
    """
    Pour nc_wt :
    garde les lignes qui n'ont pas de cousin I/L différent dans la référence non-mutante.
    """
    if PEPTIDE_COL not in df.columns:
        raise ValueError(
            f"Colonne '{PEPTIDE_COL}' absente. Colonnes disponibles : {list(df.columns)}"
        )

    remove_mask = df[PEPTIDE_COL].map(
        lambda pep: has_different_il_cousin_in_ref(pep, ref_il_to_exact)
    )

    filtered = df.loc[~remove_mask].copy()

    return filtered, int(remove_mask.sum())


# ============================================================
# File processing
# ============================================================

def read_tsv(path: Path) -> pd.DataFrame:
    """Lit un TSV avec vérification d'existence."""
    if not path.exists():
        raise FileNotFoundError(f"Fichier d'entrée introuvable : {path}")
    return pd.read_csv(path, sep="\t")


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    """Écrit un TSV en créant le dossier parent si nécessaire."""
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def copy_canonical_as_is(
    *,
    cell: str,
    input_root: Path,
    final_root: Path,
) -> dict[str, object]:
    """
    canonical.tsv est copié tel quel.
    Aucun filtre gene-mismatch.
    Aucun filtre I/L.
    """
    src = input_root / cell / "01_basic_cleanup" / "canonical_final_merged_non_mutant.tsv"
    dst = final_root / cell / "02_IL_filtering" / "canonical.tsv"

    if not src.exists():
        raise FileNotFoundError(f"Source canonical introuvable : {src}")

    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)

    n_rows = len(pd.read_csv(src, sep="\t"))

    print(f"[COPY AS IS] {src.name} -> {dst}")
    print(
        f"[COUNT] before={n_rows:,} | "
        f"removed_gene_mismatch=0 | removed_IL=0 | kept={n_rows:,}"
    )

    return {
        "cell": cell,
        "output_type": "canonical",
        "input": str(src),
        "output": str(dst),
        "filters_applied": "none_copy_as_is",
        "n_before": n_rows,
        "n_removed_gene_mismatch": 0,
        "n_removed_IL": 0,
        "n_after": n_rows,
    }


def process_file(
    *,
    cell: str,
    output_type: str,
    input_root: Path,
    final_root: Path,
    input_filename: str,
    output_filename: str,
    ref_il_set: set[str],
    ref_il_to_exact: dict[str, set[str]],
    il_mode: str,
) -> dict[str, object]:
    """
    Traite un fichier variant ou non-canonique.

    il_mode:
        variant_vs_nonmut_ref -> retire si peptide_IL existe en référence non-mutante
        nonmut_cousin         -> retire si cousin I/L exact différent existe
    """
    src = input_root / cell / "01_basic_cleanup" / input_filename
    dst = final_root / cell / "02_IL_filtering" / output_filename

    print(f"[INPUT] {src}")

    df = read_tsv(src)
    n_before = len(df)

    df_gene, n_removed_gene = apply_gene_consistency_filter(df)

    if il_mode == "variant_vs_nonmut_ref":
        df_il, n_removed_il = filter_variant_against_nonmut_ref(df_gene, ref_il_set)
    elif il_mode == "nonmut_cousin":
        df_il, n_removed_il = filter_nonmutant_il_cousins(df_gene, ref_il_to_exact)
    else:
        raise ValueError(f"Mode I/L invalide : {il_mode}")

    write_tsv(df_il, dst)

    n_after = len(df_il)
    filters = f"gene_consistency_ignore_missing_GN;IL_{il_mode}"

    print(f"[OUTPUT] {dst}")
    print(
        f"[COUNT] before={n_before:,} | "
        f"removed_gene_mismatch={n_removed_gene:,} | "
        f"removed_IL={n_removed_il:,} | kept={n_after:,}"
    )

    return {
        "cell": cell,
        "output_type": output_type,
        "input": str(src),
        "output": str(dst),
        "filters_applied": filters,
        "n_before": n_before,
        "n_removed_gene_mismatch": n_removed_gene,
        "n_removed_IL": n_removed_il,
        "n_after": n_after,
    }


def get_digested_file(config: dict, cell: str) -> Path:
    """
    Récupère le fichier digéré associé à une lignée.

    Format attendu dans config.yaml :

    path:
      digested_files: "/path/to/fasta_digéré"

    Le fichier attendu est :
      {cell}_contam_trypsin_l9_mc2_digested.csv
    """

    if "path" not in config:
        raise KeyError("Section 'path' absente du config.yaml.")

    if "digested_files" not in config["path"]:
        raise KeyError(
            "Chemin 'path: digested_files:' absent du config.yaml."
        )

    digested_files = config["path"]["digested_files"]

    # Cas 1 : digested_files est un dossier commun
    if isinstance(digested_files, str):
        digested_dir = Path(digested_files)
        digested_csv = digested_dir / f"{cell}_contam_trypsin_l9_mc2_digested.csv"

        if not digested_csv.exists():
            raise FileNotFoundError(f"Fichier digéré introuvable : {digested_csv}")

        return digested_csv

    # Cas 2 : digested_files est un dictionnaire par lignée
    if isinstance(digested_files, dict):
        if cell not in digested_files:
            raise KeyError(
                f"Aucun fichier digéré défini pour {cell}. "
                f"Clés disponibles : {list(digested_files.keys())}"
            )

        digested_csv = Path(digested_files[cell])

        if not digested_csv.exists():
            raise FileNotFoundError(f"Fichier digéré introuvable : {digested_csv}")

        return digested_csv

    raise TypeError(
        "config['path']['digested_files'] doit être soit un dossier string, "
        "soit un dictionnaire par lignée."
    )

def run_il_filtering(config: dict) -> dict[str, dict[str, Path]]:
    """
    Lance la vérification I/L pour toutes les lignées dans config.yaml.

    Paramètres attendus dans config.yaml :

    samples:
      SH-SY5Y:
        - run1
        - run2

    path:
      analysis_outdir: "/chemin/vers/analysis"

    digested_files:
      SH-SY5Y: "/chemin/vers/SH-SY5Y_contam_trypsin_l9_mc2_digested.csv"
      SK-N-Be2: "/chemin/vers/SK-N-Be2_contam_trypsin_l9_mc2_digested.csv"
    """
    analysis_outdir = Path(config["path"]["analysis_outdir"])
    input_root = analysis_outdir
    final_root = analysis_outdir
    
    final_root.mkdir(parents=True, exist_ok=True)

    summary: list[dict[str, object]] = []
    outputs: dict[str, dict[str, Path]] = {}

    for cell in config["samples"].keys():
        print("\n" + "=" * 100)
        print(f"[Étape 2] Vérification I/L — {cell}")
        print("=" * 100)

        digested_csv = get_digested_file(config, cell)
        ref_il_set, ref_il_to_exact = load_nonmutant_reference(digested_csv)

        outputs[cell] = {}

        # 1) canonical normal : copie directe.
        result = copy_canonical_as_is(
            cell=cell,
            input_root=input_root,
            final_root=final_root,
        )
        summary.append(result)
        outputs[cell]["canonical"] = Path(result["output"])

        # 2) canonical variant.
        result = process_file(
            cell=cell,
            output_type="canonical_variant",
            input_root=input_root,
            final_root=final_root,
            input_filename="canonical_final_merged_mutant_only.tsv",
            output_filename="canonical_variant.tsv",
            ref_il_set=ref_il_set,
            ref_il_to_exact=ref_il_to_exact,
            il_mode="variant_vs_nonmut_ref",
        )
        summary.append(result)
        outputs[cell]["canonical_variant"] = Path(result["output"])

        # 3) nc WT.
        result = process_file(
            cell=cell,
            output_type="nc_wt",
            input_root=input_root,
            final_root=final_root,
            input_filename="nc_final_merged_non_mutant.tsv",
            output_filename="nc_wt.tsv",
            ref_il_set=ref_il_set,
            ref_il_to_exact=ref_il_to_exact,
            il_mode="nonmut_cousin",
        )
        summary.append(result)
        outputs[cell]["nc_wt"] = Path(result["output"])

        # 4) nc variant.
        result = process_file(
            cell=cell,
            output_type="nc_variant",
            input_root=input_root,
            final_root=final_root,
            input_filename="nc_final_merged_mutant_only.tsv",
            output_filename="nc_variant.tsv",
            ref_il_set=ref_il_set,
            ref_il_to_exact=ref_il_to_exact,
            il_mode="variant_vs_nonmut_ref",
        )
        summary.append(result)
        outputs[cell]["nc_variant"] = Path(result["output"])

    summary_df = pd.DataFrame(summary)

    print(f"[FINAL_ROOT] {final_root}")

    return outputs


if __name__ == "__main__":
    raise SystemExit(
        "Ce module est conçu pour être appelé depuis main.py.\n"
        "Utilise plutôt : python main.py -c config.yaml"
    )
