#!/usr/bin/env python3

"""
Basic cleanup des résultats MSFragger/FDRBench.

Cette étape :
1. charge les fichiers canonical_final.tsv et nc_final.tsv pour chaque run;
2. fusionne les runs par lignée cellulaire;
3. ajoute les colonnes cell_line, run et sample;
4. sépare les identifications mutant-only et non-mutant;
5. écrit les fichiers merged, mutant_only et non_mutant.

Ce module est conçu pour être appelé depuis main.py.
"""

from pathlib import Path
import ast
import pandas as pd


FILES = {
    "canonical": "canonical_final.tsv",
    "nc": "nc_final.tsv",
}


def parse_protein_list(value):
    """
    Convertit la colonne protein en liste Python.

    La colonne protein peut être :
    - déjà une liste;
    - une chaîne représentant une liste;
    - une valeur vide / NaN;
    - une simple chaîne.
    """
    if isinstance(value, list):
        return value

    if pd.isna(value):
        return []

    try:
        parsed = ast.literal_eval(value)
        if isinstance(parsed, list):
            return parsed
        return [str(parsed)]
    except Exception:
        return [str(value)]


def keep_mutant_only(protein_value):
    """
    Retourne True si toutes les protéines associées sont mutantes.

    Deux formats sont reconnus :
    - '_mut_' dans l'identifiant;
    - '|type=mut' dans l'identifiant FASTA.
    """
    proteins = parse_protein_list(protein_value)

    proteins = [
        str(protein)
        for protein in proteins
        if protein is not None and str(protein).strip() not in ["", "None", "nan"]
    ]

    if len(proteins) == 0:
        return False

    return all(
        ("_mut_" in protein) or ("|type=mut" in protein)
        for protein in proteins
    )


def load_and_merge_msfragger(cell_line, runs, results_root, filename):
    """
    Charge et fusionne un fichier MSFragger/FDRBench pour tous les runs
    d'une lignée cellulaire.

    Les chemins attendus sont :
    results_root / cell_line / "ms_fragger" / "final" / run / filename
    """
    dfs = []

    for run in runs:
        path = (
            results_root
            / cell_line
            / "ms_fragger"
            / "final"
            / run
            / filename
        )

        if not path.exists():
            raise FileNotFoundError(f"Fichier introuvable : {path}")

        df = pd.read_csv(path, sep="\t")
        df["cell_line"] = cell_line
        df["run"] = run
        df["sample"] = run

        dfs.append(df)

    if not dfs:
        raise ValueError(f"Aucun fichier chargé pour {cell_line} avec {filename}")

    return pd.concat(dfs, ignore_index=True)


def split_mutant_only(merged_df, protein_col="protein"):
    """
    Sépare un DataFrame fusionné en :
    - mutant_only : toutes les protéines associées sont mutantes;
    - non_mutant : au moins une protéine associée n'est pas mutante.
    """
    if protein_col not in merged_df.columns:
        raise KeyError(
            f"Colonne '{protein_col}' absente du fichier. "
            f"Colonnes disponibles : {list(merged_df.columns)}"
        )

    mutant_mask = merged_df[protein_col].apply(keep_mutant_only)

    mutant_df = merged_df[mutant_mask].copy()
    non_mutant_df = merged_df[~mutant_mask].copy()

    return mutant_df, non_mutant_df


def process_one_file(cell_line, runs, results_root, outdir, label, filename):
    """
    Traite un seul type de fichier : canonical ou nc.

    Produit :
    - {label}_final_merged.tsv
    - {label}_final_merged_mutant_only.tsv
    - {label}_final_merged_non_mutant.tsv
    """
    print(f"\nProcessing {cell_line} — {label}: {filename}")

    merged_df = load_and_merge_msfragger(
        cell_line=cell_line,
        runs=runs,
        results_root=results_root,
        filename=filename,
    )

    mutant_df, non_mutant_df = split_mutant_only(merged_df)

    merged_out = outdir / f"{label}_final_merged.tsv"
    mutant_out = outdir / f"{label}_final_merged_mutant_only.tsv"
    non_mutant_out = outdir / f"{label}_final_merged_non_mutant.tsv"

    merged_df.to_csv(merged_out, sep="\t", index=False)
    mutant_df.to_csv(mutant_out, sep="\t", index=False)
    non_mutant_df.to_csv(non_mutant_out, sep="\t", index=False)

    print(f"[OK] merged      : {merged_out} {merged_df.shape}")
    print(f"[OK] mutant_only : {mutant_out} {mutant_df.shape}")
    print(f"[OK] non_mutant  : {non_mutant_out} {non_mutant_df.shape}")

    return {
        "merged": merged_out,
        "mutant_only": mutant_out,
        "non_mutant": non_mutant_out,
    }


def run_basic_cleanup(config):
    """
    Lance le cleanup basique pour toutes les lignées définies dans config.yaml.

    Paramètres attendus dans config.yaml :

    samples:
      SH-SY5Y:
        - run1
        - run2

    path:
      results_root: "/chemin/vers/results"
      analysis_outdir: "/chemin/vers/sorties"
    """
    results_root = Path(config["path"]["results_root"])
    analysis_outdir = Path(config["path"]["analysis_outdir"])

    cleanup_outputs = {}

    for cell_line, runs in config["samples"].items():
        print("\n========================================")
        print(f"[Étape 1] Basic cleanup — {cell_line}")
        print("========================================")

        cell_outdir = analysis_outdir / cell_line / "01_basic_cleanup"
        cell_outdir.mkdir(parents=True, exist_ok=True)

        cleanup_outputs[cell_line] = {}

        for label, filename in FILES.items():
            cleanup_outputs[cell_line][label] = process_one_file(
                cell_line=cell_line,
                runs=runs,
                results_root=results_root,
                outdir=cell_outdir,
                label=label,
                filename=filename,
            )

    return cleanup_outputs
