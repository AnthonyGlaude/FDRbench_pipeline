#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Étape 7 — QC final des peptides.

Cette étape combine :
    1. le support peptidique calculé à l'étape 6
    2. le scoring expérimental MS2PIP calculé à l'étape 5

Entrées attendues :

    analysis_outdir/{sample}/06_peptide_support/
        {sample}_canonical_variant_support.tsv
        {sample}_nc_variant_support.tsv
        {sample}_nc_wt_support.tsv

    analysis_outdir/{sample}/05_ms2pip_experimental_scoring/
        {sample}_canonical_variant_ndp_scores.tsv
        {sample}_nc_variant_ndp_scores.tsv
        {sample}_nc_wt_ndp_scores.tsv

Sorties :

    analysis_outdir/{sample}/07_final_QC/
        canonical_variant_final_QC_all.tsv
        canonical_variant_final_QC_kept.tsv
        nc_variant_final_QC_all.tsv
        nc_variant_final_QC_kept.tsv
        nc_wt_final_QC_all.tsv
        nc_wt_final_QC_kept.tsv
        ALL_categories_final_QC_all.tsv
        ALL_categories_final_QC_kept.tsv
        ALL_categories_final_QC_summary.tsv

Config attendue :

    final_qc:
      ndp_threshold: 0.70

Conçu pour être appelé depuis main.py :

    from modules.final_qc import run_final_qc
    final_qc_outputs = run_final_qc(config)
"""

from pathlib import Path
import re

import pandas as pd


# ============================================================
# CATEGORIES
# ============================================================

CATEGORIES = [
    "canonical_variant",
    "nc_variant",
    "nc_wt",
]


AA_ONLY = re.compile(r"[A-Z]")


# ============================================================
# PATH HELPERS
# ============================================================

def get_analysis_outdir(config: dict) -> Path:
    """
    Dossier principal des analyses modulaires.

    Exemple config :
        path:
          analysis_outdir: results
    """
    path_cfg = config.get("path", {})

    if "analysis_outdir" not in path_cfg:
        raise KeyError(
            "Chemin analysis_outdir absent du config.yaml. Ajoute :\n"
            "path:\n"
            "  analysis_outdir: results"
        )

    return Path(path_cfg["analysis_outdir"])


def get_ndp_threshold(config: dict) -> float:
    """
    Seuil NDP utilisé pour le QC final.

    Exemple config :
        final_qc:
          ndp_threshold: 0.70
    """
    final_qc_cfg = config.get("final_qc", {})

    return float(final_qc_cfg.get("ndp_threshold", 0.70))


def get_support_dir(config: dict, sample: str) -> Path:
    """
    Dossier de sortie de l'étape 6.

    Attendu :
        analysis_outdir/{sample}/06_peptide_support/
    """
    analysis_outdir = get_analysis_outdir(config)

    return analysis_outdir / sample / "06_peptide_support"


def get_ndp_dir(config: dict, sample: str) -> Path:
    """
    Dossier de sortie de l'étape 5.

    Attendu :
        analysis_outdir/{sample}/05_ms2pip_experimental_scoring/
    """
    analysis_outdir = get_analysis_outdir(config)

    return analysis_outdir / sample / "05_ms2pip_experimental_scoring"


def get_final_qc_outdir(config: dict, sample: str) -> Path:
    """
    Dossier de sortie du QC final.

    Sortie :
        analysis_outdir/{sample}/07_final_QC/
    """
    analysis_outdir = get_analysis_outdir(config)

    return analysis_outdir / sample / "07_final_QC"


def support_path(config: dict, sample: str, category: str) -> Path:
    """
    Fichier support peptidique.

    Exemple :
        results/SH-SY5Y/06_peptide_support/SH-SY5Y_nc_variant_support.tsv
    """
    return get_support_dir(config, sample) / f"{sample}_{category}_support.tsv"


def ndp_path(config: dict, sample: str, category: str) -> Path:
    """
    Fichier NDP.

    Exemple :
        results/SH-SY5Y/05_ms2pip_experimental_scoring/SH-SY5Y_nc_variant_ndp_scores.tsv
    """
    return get_ndp_dir(config, sample) / f"{sample}_{category}_ndp_scores.tsv"


# ============================================================
# HELPERS
# ============================================================

def pep_key(x) -> str:
    """
    Normalise les peptides :
    - garde seulement A-Z
    - uppercase
    - I/L -> J
    """
    if x is None:
        return ""

    try:
        if pd.isna(x):
            return ""
    except Exception:
        pass

    seq = "".join(AA_ONLY.findall(str(x).upper()))

    return seq.replace("I", "J").replace("L", "J")


def load_tsv(path: Path) -> pd.DataFrame:
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"Fichier introuvable: {path}")

    if path.stat().st_size == 0:
        raise ValueError(f"Fichier vide: {path}")

    return pd.read_csv(path, sep="\t")


def prepare_support(df: pd.DataFrame, path: Path) -> pd.DataFrame:
    """
    Prépare le fichier support pour le merge.
    """
    df = df.copy()

    required = [
        "cell_line",
        "category",
        "spectrum",
        "peptide",
        "support_decision",
    ]

    missing = [c for c in required if c not in df.columns]

    if missing:
        raise ValueError(
            f"Colonnes manquantes dans support file {path}: {missing}"
        )

    df["support_pep_key"] = df["peptide"].map(pep_key)

    return df


def prepare_ndp(df: pd.DataFrame, category: str, path: Path) -> pd.DataFrame:
    """
    Prépare le fichier NDP pour le merge.
    """
    df = df.copy()

    required = [
        "cell_line",
        "spectrum",
        "modified_peptide",
        "ndp_mean",
    ]

    missing = [c for c in required if c not in df.columns]

    if missing:
        raise ValueError(
            f"Colonnes manquantes dans NDP file {path}: {missing}"
        )

    df["category"] = category
    df["ndp_pep_key"] = df["modified_peptide"].map(pep_key)

    numeric_cols = [
        "ndp_exp_ref",
        "ndp_theo_ref",
        "ndp_mean",
        "n_exp_annotated",
        "n_theo_ions",
        "n_shared_ions",
        "fraction_theo_found",
        "peptideprophet_probability",
        "hyperscore",
        "massdiff",
    ]

    for c in numeric_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    return df


def call_final_decision(row, ndp_threshold: float) -> str:
    """
    Décision finale :
    - KEEP si support OK et NDP OK
    - REJECT_SUPPORT si support échoue seulement
    - REJECT_NDP si NDP échoue seulement
    - REJECT_BOTH si les deux échouent
    """
    support_ok = row["support_decision"] == "KEEP"

    ndp = row.get("ndp_mean")
    ndp_ok = pd.notna(ndp) and ndp >= ndp_threshold

    if support_ok and ndp_ok:
        return "KEEP"

    if not support_ok and not ndp_ok:
        return "REJECT_BOTH"

    if not support_ok:
        return "REJECT_SUPPORT"

    return "REJECT_NDP"


def call_final_reject_reason(row, ndp_threshold: float) -> str:
    """
    Raison du rejet final.
    """
    final_decision = row["final_decision"]

    if final_decision == "KEEP":
        return ""

    if final_decision == "REJECT_SUPPORT":
        return str(row.get("support_reject_reason", "support_qc_failed"))

    if final_decision == "REJECT_NDP":
        return f"ndp_mean_below_{ndp_threshold}"

    if final_decision == "REJECT_BOTH":
        support_reason = str(
            row.get("support_reject_reason", "support_qc_failed")
        )
        return f"{support_reason};ndp_mean_below_{ndp_threshold}"

    return "unknown"


# ============================================================
# MERGE
# ============================================================

def merge_support_ndp(
    config: dict,
    sample: str,
    category: str,
    ndp_threshold: float,
) -> pd.DataFrame:
    print(f"\n=== {sample} | {category} ===")

    spath = support_path(config, sample, category)
    npath = ndp_path(config, sample, category)

    print(f"[1] Support: {spath}")
    support = prepare_support(load_tsv(spath), path=spath)

    print(f"[2] NDP: {npath}")
    ndp = prepare_ndp(load_tsv(npath), category=category, path=npath)

    print("[3] Merge support + NDP on cell_line + category + spectrum")

    merged = support.merge(
        ndp,
        on=["cell_line", "category", "spectrum"],
        how="left",
        suffixes=("_support", "_ndp"),
    )

    merged["peptide_match_support_ndp"] = (
        merged["support_pep_key"] == merged["ndp_pep_key"]
    )

    merged["ndp_missing"] = merged["ndp_mean"].isna()
    merged["ndp_pass"] = merged["ndp_mean"] >= ndp_threshold

    merged["final_decision"] = merged.apply(
        lambda row: call_final_decision(row, ndp_threshold),
        axis=1,
    )

    merged["final_reject_reason"] = merged.apply(
        lambda row: call_final_reject_reason(row, ndp_threshold),
        axis=1,
    )

    merged["final_keep"] = merged["final_decision"] == "KEEP"

    print("--- merge check ---")
    print("support rows:", len(support))
    print("merged rows :", len(merged))
    print("with NDP    :", merged["ndp_mean"].notna().sum())
    print("missing NDP :", merged["ndp_mean"].isna().sum())
    print("pep matches :", merged["peptide_match_support_ndp"].sum(), "/", len(merged))

    print("--- final decision ---")
    print(merged["final_decision"].value_counts(dropna=False))

    return merged


# ============================================================
# MAIN PIPELINE FUNCTION
# ============================================================

def run_final_qc(config: dict) -> dict[str, Path]:
    """
    Étape 7 — QC final.

    Combine :
        - support peptidique
        - score NDP MS2PIP

    Sorties par sample :
        analysis_outdir/{sample}/07_final_QC/

    Retourne :
        dictionnaire des fichiers écrits.
    """
    ndp_threshold = get_ndp_threshold(config)

    print("\n" + "#" * 100)
    print("[Étape 7] QC final des peptides")
    print("#" * 100)
    print(f"[CONFIG] ndp_threshold: {ndp_threshold}")

    outputs: dict[str, Path] = {}
    all_results = []

    for sample in config["samples"]:
        print("\n" + "=" * 100)
        print(f"[SAMPLE] {sample}")
        print("=" * 100)

        cell_outdir = get_final_qc_outdir(config, sample)
        cell_outdir.mkdir(parents=True, exist_ok=True)

        print(f"[OUTPUT] {cell_outdir}")

        sample_results = []

        for category in CATEGORIES:
            merged = merge_support_ndp(
                config=config,
                sample=sample,
                category=category,
                ndp_threshold=ndp_threshold,
            )

            out_full = cell_outdir / f"{category}_final_QC_all.tsv"
            out_keep = cell_outdir / f"{category}_final_QC_kept.tsv"

            merged.to_csv(out_full, sep="\t", index=False)

            kept = merged[merged["final_decision"] == "KEEP"].copy()
            kept.to_csv(out_keep, sep="\t", index=False)

            print(f"[WRITE] {out_full}")
            print(f"[WRITE] {out_keep}")

            outputs[f"{sample}_{category}_all"] = out_full
            outputs[f"{sample}_{category}_kept"] = out_keep

            sample_results.append(merged)
            all_results.append(merged)

        if not sample_results:
            print(f"[WARN] Aucun résultat produit pour {sample}")
            continue

        sample_df = pd.concat(sample_results, ignore_index=True)

        sample_all_full = cell_outdir / "ALL_categories_final_QC_all.tsv"
        sample_all_keep = cell_outdir / "ALL_categories_final_QC_kept.tsv"
        sample_summary = cell_outdir / "ALL_categories_final_QC_summary.tsv"

        sample_df.to_csv(sample_all_full, sep="\t", index=False)

        sample_kept = sample_df[sample_df["final_decision"] == "KEEP"].copy()
        sample_kept.to_csv(sample_all_keep, sep="\t", index=False)

        summary = (
            sample_df
            .groupby(
                [
                    "cell_line",
                    "category",
                    "support_decision",
                    "ndp_pass",
                    "final_decision",
                ],
                dropna=False,
            )
            .size()
            .reset_index(name="n")
        )

        summary.to_csv(sample_summary, sep="\t", index=False)

        print(f"[WRITE] {sample_all_full}")
        print(f"[WRITE] {sample_all_keep}")
        print(f"[WRITE] {sample_summary}")

        outputs[f"{sample}_all_full"] = sample_all_full
        outputs[f"{sample}_all_keep"] = sample_all_keep
        outputs[f"{sample}_summary"] = sample_summary

    if not all_results:
        raise RuntimeError("Aucun résultat produit par le QC final.")

    print("\n[DONE] Étape 7 — QC final terminé.")

    return outputs


if __name__ == "__main__":
    ()