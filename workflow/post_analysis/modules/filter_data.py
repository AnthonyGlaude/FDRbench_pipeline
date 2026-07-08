#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Étape 3 — Filtrage / résumé des données peptidiques et PSM.

Cette étape prend les sorties de l'étape 2 :

    analysis_outdir/{cell}/02_IL_filtering/
        canonical.tsv
        canonical_variant.tsv
        nc_wt.tsv
        nc_variant.tsv

et les croise avec les fichiers PeptideProphet pepXML :

    results_root/{cell}/ms_fragger/search/{run}/{run}_prophet.pepXML

Elle produit :

    analysis_outdir/{cell}/03_filter_data/
        {cell}_{category}_psm_summary.tsv

et un fichier global :

    analysis_outdir/03_filter_data/all_peptide_psm_summary.tsv

Catégories produites :
    canonical
    canonical_variant
    nc_wt
    nc_variant

Ce module est conçu pour être appelé depuis main.py :

    from module.filter_data import run_filter_data
    filter_outputs = run_filter_data(config)
"""

from __future__ import annotations

from pathlib import Path
import ast
import xml.etree.ElementTree as ET

import pandas as pd


FINAL_INPUT_FILES = {
    "canonical": "canonical.tsv",
    "canonical_variant": "canonical_variant.tsv",
    "nc_wt": "nc_wt.tsv",
    "nc_variant": "nc_variant.tsv",
}


# ============================================================
# Helpers pepXML
# ============================================================

def strip_namespace(tag: str) -> str:
    """Retire le namespace XML d'un tag."""
    return tag.split("}", 1)[-1] if "}" in tag else tag


def read_prophet_pepxml(pepxml_path: Path, sample_name: str) -> pd.DataFrame:
    """
    Lit un fichier PeptideProphet pepXML et retourne les PSM top-rank.

    Colonnes retournées :
        sample
        spectrum
        peptide
        modified_peptide
        charge
        peptideprophet_probability
    """
    rows = []

    if not pepxml_path.exists():
        raise FileNotFoundError(f"pepXML introuvable : {pepxml_path}")

    context = ET.iterparse(pepxml_path, events=("end",))

    for _event, elem in context:
        if strip_namespace(elem.tag) != "spectrum_query":
            continue

        spectrum = elem.attrib.get("spectrum")
        charge = elem.attrib.get("assumed_charge")

        for child in elem.iter():
            if strip_namespace(child.tag) != "search_hit":
                continue

            if str(child.attrib.get("hit_rank")) != "1":
                continue

            peptide = child.attrib.get("peptide")
            modified_peptide = peptide
            probability = None

            for subchild in child.iter():
                tag = strip_namespace(subchild.tag)

                if tag == "modification_info":
                    modified_peptide = subchild.attrib.get("modified_peptide", peptide)

                elif tag == "peptideprophet_result":
                    probability = subchild.attrib.get("probability")

            rows.append(
                {
                    "sample": sample_name,
                    "spectrum": spectrum,
                    "peptide": peptide,
                    "modified_peptide": modified_peptide,
                    "charge": int(charge) if charge is not None else None,
                    "peptideprophet_probability": (
                        float(probability) if probability is not None else None
                    ),
                }
            )

        elem.clear()

    return pd.DataFrame(rows)


def get_pepxml_path(config: dict, cell_line: str, run: str) -> Path:
    """
    Retourne le chemin pepXML pour un run.

    Priorité :
    1. config["prophet_files"][cell_line][run], si défini explicitement;
    2. path["pepxml_root"], si défini;
    3. path["results_root"] par défaut.

    Structure attendue par défaut :
        root/{cell_line}/ms_fragger/search/{run}/{run}_prophet.pepXML
    """
    if "prophet_files" in config:
        cell_prophets = config["prophet_files"].get(cell_line, {})
        if run in cell_prophets:
            return Path(cell_prophets[run])

    path_cfg = config.get("path", {})

    if "pepxml_root" in path_cfg:
        root = Path(path_cfg["pepxml_root"])
    else:
        root = Path(path_cfg["results_root"])

    return (
        root
        / cell_line
        / "ms_fragger"
        / "search"
        / run
        / f"{run}_prophet.pepXML"
    )


def summarize_prophet_psms(
    *,
    config: dict,
    cell_line: str,
    runs: list[str],
    min_probability: float | None = None,
) -> pd.DataFrame:
    """
    Résume les PSM PeptideProphet pour tous les runs d'une lignée.
    """
    dfs = []

    for run in runs:
        pepxml_path = get_pepxml_path(config, cell_line, run)

        if not pepxml_path.exists():
            print(f"[MISSING pepXML] {pepxml_path}")
            continue

        print(f"[READ pepXML] {cell_line} — {run}")
        df = read_prophet_pepxml(pepxml_path, run)

        if min_probability is not None:
            df = df[df["peptideprophet_probability"] >= min_probability].copy()

        dfs.append(df)

    if not dfs:
        print(f"[WARNING] Aucun pepXML lu pour {cell_line}")
        return pd.DataFrame(
            columns=[
                "peptide",
                "modified_peptide",
                "charge",
                "n_samples_prophet",
                "n_PSM_prophet",
                "samples_prophet",
                "max_peptideprophet_probability",
                "mean_peptideprophet_probability",
            ]
        )

    all_psms = pd.concat(dfs, ignore_index=True)

    summary = (
        all_psms
        .groupby(["peptide", "modified_peptide", "charge"], dropna=False)
        .agg(
            n_samples_prophet=("sample", "nunique"),
            n_PSM_prophet=("spectrum", "size"),
            samples_prophet=("sample", lambda x: ";".join(sorted(set(map(str, x))))),
            max_peptideprophet_probability=("peptideprophet_probability", "max"),
            mean_peptideprophet_probability=("peptideprophet_probability", "mean"),
        )
        .reset_index()
    )

    return summary


# ============================================================
# Helpers final TSV
# ============================================================

def parse_list_column(x: object) -> list[object]:
    """Parse une colonne écrite comme liste Python."""
    if isinstance(x, list):
        return x

    if pd.isna(x):
        return []

    try:
        value = ast.literal_eval(str(x))
        if isinstance(value, list):
            return value
        return [str(value)]
    except Exception:
        return [str(x)]


def clean_list(values: list[object]) -> list[str]:
    """Nettoie une liste de valeurs protéiques."""
    return [
        str(v)
        for v in values
        if v is not None and str(v).strip() not in ["", "None", "nan"]
    ]


def get_charge_column(df: pd.DataFrame, path: Path) -> str:
    """Trouve la colonne de charge dans un fichier final."""
    if "assumed_charge" in df.columns:
        return "assumed_charge"
    if "charge" in df.columns:
        return "charge"

    raise ValueError(f"Aucune colonne charge trouvée dans {path}")


def get_sample_column(df: pd.DataFrame, path: Path) -> str:
    """Trouve la colonne de run/sample dans un fichier final."""
    if "sample_name" in df.columns:
        return "sample_name"
    if "replicate" in df.columns:
        return "replicate"
    if "run" in df.columns:
        return "run"
    if "sample" in df.columns:
        return "sample"

    raise ValueError(f"Aucune colonne sample/run trouvée dans {path}")


def summarize_final_file(path: Path, cell_line: str, category: str) -> pd.DataFrame:
    """
    Résume un fichier final au niveau peptide/modification/charge.

    Ajoute aussi un indicateur de support protéique :
        protein_supported_by_2plus_peptides
    """
    if not path.exists():
        raise FileNotFoundError(f"Fichier final introuvable : {path}")

    df = pd.read_csv(path, sep="\t")

    if "peptide" not in df.columns:
        raise ValueError(f"Colonne peptide absente dans {path}")

    if "protein" not in df.columns:
        raise ValueError(f"Colonne protein absente dans {path}")

    if "modified_peptide" not in df.columns:
        df["modified_peptide"] = df["peptide"]

    charge_col = get_charge_column(df, path)
    sample_col = get_sample_column(df, path)

    df["protein_list"] = df["protein"].apply(parse_list_column).apply(clean_list)
    df["protein_clean"] = df["protein_list"].apply(
        lambda x: ";".join(sorted(set(x)))
    )

    summary = (
        df.groupby(["peptide", "modified_peptide", charge_col], dropna=False)
        .agg(
            protein=(
                "protein_clean",
                lambda x: ";".join(sorted(set(";".join(x).split(";")) - {""})),
            ),
            n_proteins=(
                "protein_clean",
                lambda x: len(set(";".join(x).split(";")) - {""}),
            ),
            n_samples_final=(sample_col, "nunique"),
            n_PSM_final=("peptide", "size"),
            samples_final=(sample_col, lambda x: ";".join(sorted(set(map(str, x))))),
        )
        .reset_index()
        .rename(columns={charge_col: "charge"})
    )

    summary["cell_line"] = cell_line
    summary["category"] = category

    exploded = summary.copy()
    exploded["protein_single"] = exploded["protein"].str.split(";")
    exploded = exploded.explode("protein_single")

    exploded = exploded[
        exploded["protein_single"].notna()
        & (exploded["protein_single"].astype(str).str.strip() != "")
    ].copy()

    if exploded.empty:
        summary["protein_n_unique_peptides"] = 0
        summary["protein_supported_by_2plus_peptides"] = False
    else:
        protein_support = (
            exploded
            .groupby("protein_single", dropna=False)
            .agg(protein_n_unique_peptides=("modified_peptide", "nunique"))
            .reset_index()
        )

        exploded = exploded.merge(
            protein_support,
            on="protein_single",
            how="left",
        )

        protein_support_per_peptide = (
            exploded
            .groupby(["peptide", "modified_peptide", "charge"], dropna=False)
            .agg(protein_n_unique_peptides=("protein_n_unique_peptides", "max"))
            .reset_index()
        )

        summary = summary.merge(
            protein_support_per_peptide,
            on=["peptide", "modified_peptide", "charge"],
            how="left",
        )

        summary["protein_n_unique_peptides"] = (
            summary["protein_n_unique_peptides"]
            .fillna(0)
            .astype(int)
        )

        summary["protein_supported_by_2plus_peptides"] = (
            summary["protein_n_unique_peptides"] >= 2
        )

    return summary[
        [
            "cell_line",
            "category",
            "peptide",
            "modified_peptide",
            "charge",
            "protein",
            "n_proteins",
            "protein_n_unique_peptides",
            "protein_supported_by_2plus_peptides",
            "n_samples_final",
            "n_PSM_final",
            "samples_final",
        ]
    ]


# ============================================================
# Main function for pipeline
# ============================================================

def run_filter_data(config: dict) -> dict[str, dict[str, Path]]:
    """
    Lance l'étape 3 pour toutes les lignées dans config.yaml.

    Utilise les fichiers de l'étape 2 :
        analysis_outdir/{cell}/02_IL_filtering/*.tsv

    Produit les résumés :
        analysis_outdir/{cell}/03_filter_data/*_psm_summary.tsv
        analysis_outdir/03_filter_data/all_peptide_psm_summary.tsv
    """
    analysis_outdir = Path(config["path"]["analysis_outdir"])

    min_probability = config.get("filter_data", {}).get(
        "min_peptideprophet_probability",
        None,
    )

    outputs: dict[str, dict[str, Path]] = {}

    for cell_line, runs in config["samples"].items():
        print("\n" + "=" * 100)
        print(f"[Étape 3] Filtrage/résumé des données — {cell_line}")
        print("=" * 100)

        cell_input_dir = analysis_outdir / cell_line / "02_IL_filtering"
        cell_outdir = analysis_outdir / cell_line / "03_filter_data"
        cell_outdir.mkdir(parents=True, exist_ok=True)

        outputs[cell_line] = {}

        prophet_summary = summarize_prophet_psms(
            config=config,
            cell_line=cell_line,
            runs=runs,
            min_probability=min_probability,
        )

        for category, filename in FINAL_INPUT_FILES.items():
            final_path = cell_input_dir / filename

            if not final_path.exists():
                print(f"[MISSING FINAL] {final_path}")
                continue

            print(f"[PROCESS FINAL] {cell_line} — {category}")

            final_summary = summarize_final_file(
                path=final_path,
                cell_line=cell_line,
                category=category,
            )

            merged_summary = final_summary.merge(
                prophet_summary,
                on=["peptide", "modified_peptide", "charge"],
                how="left",
            )

            merged_summary["n_PSM_prophet"] = (
                merged_summary["n_PSM_prophet"]
                .fillna(0)
                .astype(int)
            )

            merged_summary["n_samples_prophet"] = (
                merged_summary["n_samples_prophet"]
                .fillna(0)
                .astype(int)
            )

            out = cell_outdir / f"{cell_line}_{category}_psm_summary.tsv"
            merged_summary.to_csv(out, sep="\t", index=False)

            print(f"[OK] {out} {merged_summary.shape}")

            outputs[cell_line][category] = out


    print("\n[DONE] Étape 3 terminée.")
    return outputs


if __name__ == "__main__":
    raise SystemExit(
        "Ce module est conçu pour être appelé depuis main.py.\n"
        "Utilise plutôt : python main.py -c config.yaml"
    )
