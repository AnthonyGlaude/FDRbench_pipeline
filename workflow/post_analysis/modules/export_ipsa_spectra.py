#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Étape 9 — Export de spectres pour IPSA.

Le script lit un fichier d'entrée par sample :

    ipsa_input_root/{sample}/{input_name}

Exemple :
    data/Spectra_list/SH-SY5Y/ipsa_spectra_input.tsv
    data/Spectra_list/SK-N-Be2/ipsa_spectra_input.tsv

Ce fichier doit contenir au minimum de quoi identifier le spectre :
    - soit spectrum
    - soit run_name + scan
    - soit run + scan

Colonnes utiles si présentes :
    sample / cell_line
    category
    peptide / modified_peptide_ndp / modified_peptide
    charge / assumed_charge
    run_name / run
    scan / start_scan
    spectrum
    final_keep
    final_decision

Sorties :
    analysis_outdir/{sample}/09_spectra/
        {category}/*.tsv
        plots/{category}/*.png
        extraction_summary.tsv

Comportement important :
    Si un sample n'a pas de fichier d'entrée IPSA, il est ignoré.
    Le pipeline ne crash pas.
"""

from pathlib import Path
import re
import hashlib

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyteomics import mzml


# ============================================================
# CONFIG HELPERS
# ============================================================

def get_analysis_outdir(config: dict) -> Path:
    path_cfg = config.get("path", {})

    if "analysis_outdir" not in path_cfg:
        raise KeyError(
            "Chemin analysis_outdir absent du config.yaml. Ajoute :\n"
            "path:\n"
            "  analysis_outdir: results"
        )

    return Path(path_cfg["analysis_outdir"])


def get_mzml_root(config: dict) -> Path:
    path_cfg = config.get("path", {})

    if "mzml_root" in path_cfg:
        return Path(path_cfg["mzml_root"])

    if "mzml" in path_cfg:
        return Path(path_cfg["mzml"])

    raise KeyError(
        "Chemin mzML absent du config.yaml. Ajoute :\n"
        "path:\n"
        "  mzml_root: ../data/mzml/"
    )


def get_ipsa_input_root(config: dict) -> Path:
    """
    Racine des fichiers d'entrée IPSA.

    Défaut :
        data/Spectra_list

    Donc :
        data/Spectra_list/{sample}/ipsa_spectra_input.tsv
    """
    path_cfg = config.get("path", {})
    return Path(path_cfg.get("ipsa_input_root", "data/Spectra_list"))


def get_ipsa_config(config: dict) -> dict:
    return config.get("ipsa_export", {})


def ipsa_enabled(config: dict) -> bool:
    return bool(get_ipsa_config(config).get("enabled", False))


def get_ipsa_input_name(config: dict) -> str:
    return str(get_ipsa_config(config).get("input_name", "ipsa_spectra_input.tsv"))


def get_ipsa_input_file(config: dict, sample: str) -> Path:
    """
    Fichier d'entrée par sample.

    Exemple :
        data/Spectra_list/SH-SY5Y/ipsa_spectra_input.tsv
    """
    return get_ipsa_input_root(config) / sample / get_ipsa_input_name(config)


def get_make_plots(config: dict) -> bool:
    return bool(get_ipsa_config(config).get("make_plots", True))


def get_normalize_to_basepeak(config: dict) -> bool:
    return bool(get_ipsa_config(config).get("normalize_to_basepeak", True))


def get_min_percent_basepeak(config: dict) -> float:
    return float(get_ipsa_config(config).get("min_percent_basepeak", 1.0))


def get_figsize(config: dict) -> tuple[float, float]:
    cfg = get_ipsa_config(config)

    return (
        float(cfg.get("fig_w", 12)),
        float(cfg.get("fig_h", 5)),
    )


def get_ipsa_outdir(config: dict, sample: str) -> Path:
    """
    Sortie :
        results/{sample}/09_spectra/
    """
    analysis_outdir = get_analysis_outdir(config)
    return analysis_outdir / sample / "09_spectra"


# ============================================================
# BASIC HELPERS
# ============================================================

def find_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    lower_map = {c.lower(): c for c in df.columns}

    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]

    return None


def parse_spectrum_string(s):
    """
    Parse un identifiant de spectre du type :
        raw.scan.scan.charge

    Exemple :
        20201027_LP_MB_SHSY5Y_control_R1.62158.62158.3

    Retourne :
        raw_name, scan, charge
    """
    s = str(s)

    m_scan = re.search(r"scan=(\d+)", s)

    if m_scan:
        scan = int(m_scan.group(1))
    else:
        scan = None

    m = re.match(r"(.+)\.(\d+)\.(\d+)\.(\d+)$", s)

    if m:
        raw_name = m.group(1)
        scan = int(m.group(2))
        charge = int(m.group(4))
        return raw_name, scan, charge

    return None, scan, None


def safe_name(text, max_len=80):
    text = str(text)
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text)

    if len(text) > max_len:
        h = hashlib.md5(text.encode()).hexdigest()[:8]
        text = text[:max_len] + "_" + h

    return text


def get_row_value(row, cols, default=None):
    for c in cols:
        if c in row.index:
            v = row[c]

            if pd.notna(v) and str(v).strip() != "":
                return v

    return default


def resolve_row_spectrum_info(row):
    """
    Résout raw_name, scan, charge, peptide et category à partir d'une ligne.

    Priorité :
        1. spectrum si disponible
        2. run_name/run + scan
    """

    spectrum = get_row_value(
        row,
        ["spectrum", "Spectrum", "spectrum_id", "SpecId", "spec_id"],
        default=None,
    )

    raw_name = None
    scan = None
    charge = None

    if spectrum is not None:
        raw_name, scan, charge = parse_spectrum_string(spectrum)

    if raw_name is None:
        raw_name = get_row_value(
            row,
            ["run_name", "run", "raw_name", "RawName"],
            default=None,
        )

    if scan is None:
        scan = get_row_value(
            row,
            ["scan", "start_scan", "end_scan"],
            default=None,
        )

        if scan is not None:
            scan = int(float(scan))

    if charge is None:
        charge = get_row_value(
            row,
            ["charge", "assumed_charge"],
            default=None,
        )

        if charge is not None:
            charge = int(float(charge))

    peptide = get_row_value(
        row,
        [
            "peptide",
            "modified_peptide_ndp",
            "modified_peptide",
            "modified_peptide_support",
            "Peptide",
            "sequence",
        ],
        default=f"row_{row.name}",
    )

    category = get_row_value(
        row,
        ["category", "dataset", "class_ndp", "class_support"],
        default="unknown",
    )

    return {
        "spectrum": spectrum,
        "raw_name": str(raw_name) if raw_name is not None else None,
        "scan": scan,
        "charge": charge,
        "peptide": peptide,
        "category": str(category),
    }


def read_input_file(path: Path, sample: str) -> pd.DataFrame | None:
    """
    Lit le fichier d'entrée IPSA pour un sample.

    Si le fichier est absent, on retourne None au lieu de crasher.
    Ça permet de traiter seulement les lignées qui ont un fichier d'intérêt.
    """
    if not path.exists():
        print(f"[SKIP] Aucun fichier d'entrée IPSA pour {sample}: {path}")
        return None

    df = pd.read_csv(path, sep="\t")

    if df.empty:
        print(f"[SKIP] Fichier d'entrée IPSA vide pour {sample}: {path}")
        return None

    # Si le fichier contient final_keep, on garde uniquement True.
    # Si tu veux tout inclure, ne mets pas cette colonne ou mets-la déjà filtrée.
    if "final_keep" in df.columns:
        before = len(df)
        df = df[
            df["final_keep"]
            .astype(str)
            .str.lower()
            .isin(["true", "1", "yes"])
        ].copy()
        print(f"[FILTER] final_keep TRUE : {before} -> {len(df)}")

    # Si le fichier contient final_decision, on garde KEEP.
    elif "final_decision" in df.columns:
        before = len(df)
        df = df[df["final_decision"].astype(str) == "KEEP"].copy()
        print(f"[FILTER] final_decision KEEP : {before} -> {len(df)}")

    # Si le fichier contient sample/cell_line, on garde le sample courant.
    sample_col = find_column(df, ["sample", "cell_line"])

    if sample_col is not None:
        before = len(df)
        df = df[df[sample_col].astype(str) == str(sample)].copy()
        print(f"[FILTER] {sample_col} == {sample} : {before} -> {len(df)}")

    if df.empty:
        print(f"[SKIP] Aucun spectre restant après filtres pour {sample}")
        return None

    return df


# ============================================================
# mzML EXTRACTION
# ============================================================

def find_mzml_file(config: dict, sample: str, raw_name: str) -> Path | None:
    mzml_root = get_mzml_root(config)
    mzml_file = mzml_root / sample / f"{raw_name}.mzML"

    if mzml_file.exists():
        return mzml_file

    return None


def extract_scan_from_mzml(mzml_file: Path, scan_target: int, out_file: Path):
    """
    Extrait les pics mz/intensity d'un scan donné dans un mzML.
    """
    with mzml.read(str(mzml_file)) as spectra:
        for spec in spectra:
            spec_id = spec.get("id", "")

            if spec_id.endswith(f"scan={scan_target}") or f"scan={scan_target}" in spec_id:
                mz = spec["m/z array"]
                intensity = spec["intensity array"]

                with open(out_file, "w", encoding="utf-8") as f:
                    f.write("mz\tintensity\n")

                    for m, i in zip(mz, intensity):
                        f.write(f"{m:.4f}\t{i:.6f}\n")

                return spec_id

    return None


# ============================================================
# PLOTTING
# ============================================================

def load_spectrum_tsv(tsv_file: Path) -> pd.DataFrame:
    df = pd.read_csv(tsv_file, sep="\t")

    if "mz" not in df.columns or "intensity" not in df.columns:
        raise ValueError(f"Colonnes mz/intensity absentes dans {tsv_file}")

    df = df.copy()
    df["mz"] = pd.to_numeric(df["mz"], errors="coerce")
    df["intensity"] = pd.to_numeric(df["intensity"], errors="coerce")
    df = df.dropna(subset=["mz", "intensity"])

    return df


def process_spectrum_for_plot(
    df: pd.DataFrame,
    normalize_to_basepeak: bool,
    min_percent_basepeak: float,
) -> pd.DataFrame:
    """
    Normalise les intensités et filtre les petits pics pour la figure.
    Le fichier TSV brut n'est pas modifié.
    """
    df = df.copy()

    if len(df) == 0:
        return df

    max_int = df["intensity"].max()

    if normalize_to_basepeak and max_int > 0:
        df["plot_intensity"] = df["intensity"] / max_int * 100.0
        threshold = min_percent_basepeak
    else:
        df["plot_intensity"] = df["intensity"]
        threshold = max_int * (min_percent_basepeak / 100.0) if max_int > 0 else 0

    if min_percent_basepeak > 0:
        df = df[df["plot_intensity"] >= threshold].copy()

    return df


def plot_spectrum(
    tsv_file: Path,
    out_png: Path,
    title: str,
    *,
    normalize_to_basepeak: bool,
    min_percent_basepeak: float,
    figsize: tuple[float, float],
):
    """
    Génère une figure PNG simple du spectre.
    """
    df = load_spectrum_tsv(tsv_file)

    raw_n_peaks = len(df)
    max_intensity = df["intensity"].max() if raw_n_peaks > 0 else np.nan

    plot_df = process_spectrum_for_plot(
        df,
        normalize_to_basepeak=normalize_to_basepeak,
        min_percent_basepeak=min_percent_basepeak,
    )

    kept_n_peaks = len(plot_df)

    if kept_n_peaks == 0:
        return {
            "plot_status": "empty_after_filter",
            "raw_n_peaks": raw_n_peaks,
            "plotted_n_peaks": kept_n_peaks,
            "max_intensity": max_intensity,
        }

    fig_w, fig_h = figsize

    plt.figure(figsize=(fig_w, fig_h))
    plt.vlines(plot_df["mz"], 0, plot_df["plot_intensity"], linewidth=0.8)

    plt.xlabel("m/z")

    if normalize_to_basepeak:
        plt.ylabel("Intensité relative (%)")
    else:
        plt.ylabel("Intensité")

    plt.title(title, fontsize=9)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

    return {
        "plot_status": "plotted",
        "raw_n_peaks": raw_n_peaks,
        "plotted_n_peaks": kept_n_peaks,
        "max_intensity": max_intensity,
    }


# ============================================================
# PROCESSING
# ============================================================

def process_sample(config: dict, sample: str) -> dict[str, Path]:
    input_file = get_ipsa_input_file(config, sample)
    out_dir = get_ipsa_outdir(config, sample)
    plots_dir = out_dir / "plots"

    out_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 100)
    print(f"[IPSA] {sample}")
    print("=" * 100)
    print(f"[INPUT] {input_file}")
    print(f"[OUTPUT] {out_dir}")

    df = read_input_file(input_file, sample)

    if df is None:
        return {}

    make_plots = get_make_plots(config)
    normalize_to_basepeak = get_normalize_to_basepeak(config)
    min_percent_basepeak = get_min_percent_basepeak(config)
    figsize = get_figsize(config)

    expected_raw_names = set(config["samples"][sample])
    summary_rows = []

    for idx, row in df.iterrows():
        info = resolve_row_spectrum_info(row)

        raw_name = info["raw_name"]
        scan = info["scan"]
        charge = info["charge"]
        peptide = info["peptide"]
        category = safe_name(info["category"], max_len=40)
        spectrum = info["spectrum"]

        if raw_name is None:
            print(f"[WARN] raw_name introuvable ligne {idx}")
            continue

        if scan is None:
            print(f"[WARN] scan introuvable ligne {idx}")
            continue

        if raw_name not in expected_raw_names:
            print(f"[WARN] Raw non attendu pour {sample} : {raw_name}")
            continue

        mzml_file = find_mzml_file(config, sample, raw_name)

        if mzml_file is None:
            print(f"[WARN] mzML introuvable pour {sample} | {raw_name}")
            continue

        category_out = out_dir / category
        category_plot_out = plots_dir / category

        category_out.mkdir(parents=True, exist_ok=True)
        category_plot_out.mkdir(parents=True, exist_ok=True)

        peptide_safe = safe_name(peptide)
        charge_label = f"z{charge}" if charge is not None else "zNA"

        out_file = category_out / f"{raw_name}_scan_{scan}_{charge_label}_{peptide_safe}.tsv"
        out_png = category_plot_out / f"{raw_name}_scan_{scan}_{charge_label}_{peptide_safe}.png"

        if out_file.exists():
            status = "already_exists"
            found_spec_id = "NA"
        else:
            found_spec_id = extract_scan_from_mzml(mzml_file, scan, out_file)
            status = "extracted" if found_spec_id is not None else "not_found"

        plot_info = {
            "plot_status": "not_plotted",
            "raw_n_peaks": np.nan,
            "plotted_n_peaks": np.nan,
            "max_intensity": np.nan,
        }

        if make_plots and status in ["extracted", "already_exists"] and out_file.exists():
            if out_png.exists():
                plot_info["plot_status"] = "already_exists"
            else:
                title = f"{sample} | {category} | {peptide} | scan {scan} | {charge_label}"

                plot_info = plot_spectrum(
                    tsv_file=out_file,
                    out_png=out_png,
                    title=title,
                    normalize_to_basepeak=normalize_to_basepeak,
                    min_percent_basepeak=min_percent_basepeak,
                    figsize=figsize,
                )

        summary_rows.append({
            "sample": sample,
            "cell_line": sample,
            "category": category,
            "input_row_index": idx,
            "peptide": peptide,
            "spectrum": spectrum,
            "raw_name": raw_name,
            "scan": scan,
            "charge": charge,
            "mzml_file": str(mzml_file),
            "peak_file": str(out_file) if status in ["extracted", "already_exists"] else "",
            "plot_file": str(out_png) if plot_info["plot_status"] in ["plotted", "already_exists"] else "",
            "status": status,
            "mzml_spectrum_id": found_spec_id,
            **plot_info,
        })

        if status == "extracted":
            print(f"[OK] {sample} | {category} | scan {scan} | {charge_label} | {peptide}")
        elif status == "already_exists":
            print(f"[SKIP] déjà extrait | {sample} | {category} | scan {scan} | {charge_label} | {peptide}")
        elif status == "not_found":
            print(f"[WARN] Scan {scan} non trouvé dans {mzml_file.name}")

    if not summary_rows:
        print(f"[SKIP] Aucun spectre extrait pour {sample}")
        return {}

    summary = pd.DataFrame(summary_rows)
    summary_file = out_dir / "extraction_summary.tsv"
    summary.to_csv(summary_file, sep="\t", index=False)

    print(f"\n[DONE] Extraction terminée pour {sample}.")
    print(f"[WRITE] Résumé : {summary_file}")

    if make_plots:
        print(f"[INFO] Figures : {plots_dir}")

    return {
        f"{sample}_spectra_summary": summary_file,
    }


# ============================================================
# MAIN RUNNER
# ============================================================

def run_ipsa_spectrum_export(config: dict) -> dict[str, Path]:
    """
    Étape 9 — Export optionnel des spectres pour IPSA.

    Contrôle :
        ipsa_export.enabled

    Input :
        ipsa_input_root/{sample}/{input_name}

    Output :
        analysis_outdir/{sample}/09_spectra/
    """
    if not ipsa_enabled(config):
        print("\n[SKIP] Étape 9 — Export IPSA désactivé dans config.yaml.")
        return {}

    print("\n" + "#" * 100)
    print("[Étape 9] Export des spectres pour IPSA")
    print("#" * 100)

    outputs = {}

    for sample in config["samples"]:
        sample_outputs = process_sample(config, sample)
        outputs.update(sample_outputs)

    if not outputs:
        print("\n[SKIP] Aucun spectre IPSA exporté. Aucun fichier d'entrée trouvé ou aucun spectre valide.")

    print("\n[DONE] Étape 9 — Export IPSA terminé.")

    return outputs


if __name__ == "__main__":
    raise SystemExit(
        "Ce module est conçu pour être appelé depuis main.py.\n"
        "Utilise plutôt : python main.py -c config.yaml"
    )