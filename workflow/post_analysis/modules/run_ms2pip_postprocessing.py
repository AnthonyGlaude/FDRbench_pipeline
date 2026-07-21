"""
Étape 5 — Scoring expérimental des prédictions MS2PIP.

Cette étape est lancée après l'exécution externe de MS2PIP.

Entrées attendues :

    analysis_outdir/{cell}/02_IL_filtering/
        canonical.tsv
        canonical_variant.tsv
        nc_wt.tsv
        nc_variant.tsv

    ms2pip_output_root/{cell}/
        canonical_wt_ms2pip.spectronaut.tsv
        c_mut_ms2pip.spectronaut.tsv
        nc_wt_ms2pip.spectronaut.tsv
        nc_mut_ms2pip.spectronaut.tsv

    mzml_root/{cell}/
        fichiers *.mzML

Sorties :

    analysis_outdir/{cell}/05_ms2pip_experimental_scoring/
        {cell}_{category}_ndp_scores.tsv
        {cell}_{category}_ndp_mean_distribution.png
        {cell}_{category}_ndp_mean_stats.tsv

Conçu pour être appelé depuis main.py :

    from modules.run_ms2pip_postprocessing import run_ms2pip_experimental_scoring
    ndp_outputs = run_ms2pip_experimental_scoring(config)
"""

from __future__ import annotations

from pathlib import Path
import re

import numpy as np
import pandas as pd
from pyteomics import mzml

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


FILE_PAIRS = {
    "canonical": {
        "psm": "canonical.tsv",
        "ms2pip_candidates": [
            "canonical_wt_ms2pip.spectronaut.tsv",
        ],
    },
    "canonical_variant": {
        "psm": "canonical_variant.tsv",
        "ms2pip_candidates": [
            "c_mut_ms2pip.spectronaut.tsv",
        ],
    },
    "nc_variant": {
        "psm": "nc_variant.tsv",
        "ms2pip_candidates": [
            "nc_mut_ms2pip.spectronaut.tsv",
        ],
    },
    "nc_wt": {
        "psm": "nc_wt.tsv",
        "ms2pip_candidates": [
            "nc_wt_ms2pip.spectronaut.tsv",
        ],
    },
}


# ============================================================
# Path helpers
# ============================================================

def get_analysis_outdir(config: dict) -> Path:
    """
    Récupère le dossier principal des analyses.

    Format attendu :
        path:
          analysis_outdir: results
    """
    path_cfg = config.get("path", {})

    if "analysis_outdir" not in path_cfg:
        raise KeyError(
            "Chemin analysis_outdir absent du config.yaml. Ajoute par exemple :\n"
            "path:\n"
            "  analysis_outdir: results"
        )

    return Path(path_cfg["analysis_outdir"])


def get_ms2pip_output_root(config: dict) -> Path:
    """
    Récupère le dossier racine où déposer les sorties MS2PIP du serveur.

    Format attendu :
        path:
          ms2pip_output_root: ../data/MS2PIP_intensities
    """
    path_cfg = config.get("path", {})

    if "ms2pip_output_root" in path_cfg:
        return Path(path_cfg["ms2pip_output_root"])

    if "ms2pip_outdir" in path_cfg:
        return Path(path_cfg["ms2pip_outdir"])

    raise KeyError(
        "Chemin MS2PIP absent du config.yaml. Ajoute par exemple :\n"
        "path:\n"
        "  ms2pip_output_root: ../data/MS2PIP_intensities"
    )


def get_ms2pip_dir(config: dict, cell_line: str) -> Path:
    """
    Retourne le dossier contenant les sorties MS2PIP pour une lignée.

    Structure attendue :
        ms2pip_output_root/{cell_line}/

    Exemple :
        ../data/MS2PIP_intensities/SK-N-Be2/
        ../data/MS2PIP_intensities/SH-SY5Y/
    """
    path_cfg = config.get("path", {})

    # Option avancée : chemin spécifique par lignée.
    if "ms2pip_output_files" in path_cfg:
        files_cfg = path_cfg["ms2pip_output_files"]

        if isinstance(files_cfg, dict) and cell_line in files_cfg:
            return Path(files_cfg[cell_line])

    ms2pip_output_root = get_ms2pip_output_root(config)

    return ms2pip_output_root / cell_line


def get_mzml_root(config: dict) -> Path:
    """
    Récupère le dossier racine mzML.

    Format attendu :
        path:
          mzml_root: ../data/mzml/
    """
    path_cfg = config.get("path", {})

    if "mzml_root" in path_cfg:
        return Path(path_cfg["mzml_root"])

    if "mzml" in path_cfg:
        return Path(path_cfg["mzml"])

    raise KeyError(
        "Chemin mzML absent du config.yaml. Ajoute par exemple :\n"
        "path:\n"
        "  mzml_root: ../data/mzml/"
    )


# ============================================================
# Basic utils
# ============================================================

def get_scan_number(spec_id: object) -> int | None:
    match = re.search(r"scan=(\d+)", str(spec_id))
    return int(match.group(1)) if match else None


def extract_run_name_from_spectrum(spectrum: object) -> str:
    """
    Exemple :
        20201027_LP_MB_SHSY5Y_control_R1.43742.43742.2
        -> 20201027_LP_MB_SHSY5Y_control_R1
    """
    spectrum = str(spectrum)
    parts = spectrum.rsplit(".", 3)

    if len(parts) == 4:
        return parts[0]

    return spectrum


def find_existing_ms2pip_file(ms2pip_dir: Path, candidates: list[str]) -> Path:
    """
    Retourne le premier fichier MS2PIP existant et non vide.
    """
    ms2pip_dir = Path(ms2pip_dir)

    for name in candidates:
        p = ms2pip_dir / name
        if p.exists() and p.stat().st_size > 0:
            return p

    existing = [
        f"{p.name} ({p.stat().st_size} bytes)"
        for p in sorted(ms2pip_dir.glob("*.tsv"))
    ]

    raise FileNotFoundError(
        f"Aucun fichier MS2PIP valide trouvé dans {ms2pip_dir}\n"
        f"Candidats testés : {candidates}\n"
        f"Fichiers présents : {existing}"
    )


def build_mzml_index(mzml_dir: Path) -> dict[str, Path]:
    """
    Construit un index :
        run_name -> mzML path
    """
    mzml_dir = Path(mzml_dir)
    mzml_index: dict[str, Path] = {}

    for mzml_path in sorted(mzml_dir.glob("*.mzML")):
        mzml_index[mzml_path.stem] = mzml_path

    return mzml_index


def load_mzml_ms2_scans(mzml_path: Path) -> dict[int, pd.DataFrame]:
    """
    Charge tous les scans MS2 d'un mzML.

    Retourne :
        scan_number -> DataFrame(mz, intensity)
    """
    scan_dict: dict[int, pd.DataFrame] = {}

    with mzml.read(str(mzml_path)) as reader:
        for spec in reader:
            if spec.get("ms level") != 2:
                continue

            scan_number = get_scan_number(spec.get("id", ""))

            if scan_number is None:
                continue

            scan_dict[scan_number] = pd.DataFrame({
                "mz": spec["m/z array"],
                "intensity": spec["intensity array"],
            })

    return scan_dict


def load_ms2pip_library(ms2pip_file: Path) -> dict[tuple[str, int], pd.DataFrame]:
    """
    Charge un fichier MS2PIP Spectronaut.

    Retourne :
        (ModifiedPeptide, PrecursorCharge) -> fragments DataFrame
    """
    ms2pip_file = Path(ms2pip_file)

    if not ms2pip_file.exists():
        raise FileNotFoundError(f"MS2PIP file not found: {ms2pip_file}")

    if ms2pip_file.stat().st_size == 0:
        raise ValueError(f"MS2PIP file is empty: {ms2pip_file}")

    df = pd.read_csv(ms2pip_file, sep="\t")

    if df.empty:
        raise ValueError(f"MS2PIP file has header but no rows: {ms2pip_file}")

    required_cols = {
        "ModifiedPeptide",
        "PrecursorCharge",
        "RelativeFragmentIntensity",
        "FragmentMz",
        "FragmentType",
        "FragmentNumber",
        "FragmentCharge",
    }

    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Colonnes manquantes dans {ms2pip_file}: {missing}")

    lib: dict[tuple[str, int], pd.DataFrame] = {}

    for (pep, charge), group in df.groupby(
        ["ModifiedPeptide", "PrecursorCharge"],
        sort=False,
    ):
        lib[(str(pep), int(charge))] = group.copy()

    return lib


def get_charge_column(df: pd.DataFrame, path: Path) -> str:
    if "assumed_charge" in df.columns:
        return "assumed_charge"

    if "charge" in df.columns:
        return "charge"

    raise ValueError(f"Aucune colonne de charge trouvée dans {path}")


def get_scan_column(df: pd.DataFrame, path: Path) -> str:
    if "start_scan" in df.columns:
        return "start_scan"

    if "scan" in df.columns:
        return "scan"

    raise ValueError(f"Aucune colonne scan/start_scan trouvée dans {path}")


def get_modified_peptide_column(df: pd.DataFrame) -> str:
    if "modified_peptide" in df.columns:
        return "modified_peptide"

    if "peptide" in df.columns:
        return "peptide"

    raise ValueError("Aucune colonne modified_peptide ou peptide trouvée.")


# ============================================================
# NDP functions
# ============================================================

def match_theoretical_to_experimental(
    theo_df: pd.DataFrame,
    exp_df: pd.DataFrame,
    tolerance: float = 0.02,
) -> pd.DataFrame:
    """
    Pour chaque fragment théorique MS2PIP, cherche le pic expérimental
    le plus proche dans le scan mzML.
    """
    rows = []

    exp_mz = exp_df["mz"].to_numpy(float)
    exp_int = exp_df["intensity"].to_numpy(float)

    if len(exp_mz) == 0:
        return pd.DataFrame()

    for _, row in theo_df.iterrows():
        theo_mz = float(row["FragmentMz"])

        diffs = np.abs(exp_mz - theo_mz)
        best_idx = int(np.argmin(diffs))
        best_diff = float(diffs[best_idx])

        ion = (
            f"{row['FragmentType']}"
            f"{int(row['FragmentNumber'])}"
            f"^{int(row['FragmentCharge'])}"
        )

        if best_diff <= tolerance:
            rows.append({
                "ion": ion,
                "theoretical_mz": theo_mz,
                "theoretical_intensity": float(row["RelativeFragmentIntensity"]),
                "experimental_mz": float(exp_mz[best_idx]),
                "experimental_intensity": float(exp_int[best_idx]),
                "mz_error": best_diff,
                "matched": True,
            })
        else:
            rows.append({
                "ion": ion,
                "theoretical_mz": theo_mz,
                "theoretical_intensity": float(row["RelativeFragmentIntensity"]),
                "experimental_mz": np.nan,
                "experimental_intensity": 0.0,
                "mz_error": np.nan,
                "matched": False,
            })

    return pd.DataFrame(rows)


def build_dicts_from_matched(
    matched: pd.DataFrame,
) -> tuple[dict[str, list[float]], dict[str, list[float]]]:
    """
    theo_dict:
        tous les fragments théoriques
        ion -> [theoretical_mz, theoretical_intensity]

    fragments_dict:
        seulement les fragments expérimentaux matchés
        ion -> [experimental_intensity, experimental_mz]
    """
    theo_dict = {}
    fragments_dict = {}

    for _, row in matched.iterrows():
        ion = row["ion"]

        theo_dict[ion] = [
            float(row["theoretical_mz"]),
            float(row["theoretical_intensity"]),
        ]

        if bool(row["matched"]):
            fragments_dict[ion] = [
                float(row["experimental_intensity"]),
                float(row["experimental_mz"]),
            ]

    return fragments_dict, theo_dict


def safe_max_normalize(values: list[float] | np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)

    if len(values) == 0:
        return values

    max_val = np.max(values)

    if max_val == 0:
        return values

    return np.array([val / max_val * 100 for val in values], dtype=float)


def normalized_dot_product(
    exp_intensity: list[float],
    theo_intensity: list[float],
) -> float:
    """
    Normalise les intensités au pic maximal, puis calcule le produit scalaire
    entre vecteurs unitaires.
    """
    norm_exp = safe_max_normalize(exp_intensity)
    norm_theo = safe_max_normalize(theo_intensity)

    norm1 = np.linalg.norm(norm_exp)
    norm2 = np.linalg.norm(norm_theo)

    if norm1 == 0 or norm2 == 0:
        return 0.0

    unit_vec1 = norm_exp / norm1
    unit_vec2 = norm_theo / norm2

    return float(np.dot(unit_vec1, unit_vec2))


def ndp_exp_reference(
    fragments_dict: dict[str, list[float]],
    theo_dict: dict[str, list[float]],
) -> float:
    """
    Analyse 1 :
        référence = fragments expérimentaux annotés.
    """
    exp_intensity = [float(val[0]) for val in fragments_dict.values()]

    theo_intensity = []

    for ion in fragments_dict:
        if ion in theo_dict:
            theo_intensity.append(float(theo_dict[ion][1]))
        else:
            theo_intensity.append(0.0)

    return normalized_dot_product(exp_intensity, theo_intensity)


def ndp_theo_reference(
    fragments_dict: dict[str, list[float]],
    theo_dict: dict[str, list[float]],
) -> float:
    """
    Analyse 2 :
        référence = fragments théoriques.
    """
    theo_intensity = [float(val[1]) for val in theo_dict.values()]

    exp_intensity = []

    for ion in theo_dict:
        if ion in fragments_dict:
            exp_intensity.append(float(fragments_dict[ion][0]))
        else:
            exp_intensity.append(0.0)

    return normalized_dot_product(exp_intensity, theo_intensity)


def compute_two_way_ndp(
    fragments_dict: dict[str, list[float]],
    theo_dict: dict[str, list[float]],
) -> dict[str, float | int]:
    ndp_exp_ref = ndp_exp_reference(fragments_dict, theo_dict)
    ndp_theo_ref = ndp_theo_reference(fragments_dict, theo_dict)

    ndp_mean = (ndp_exp_ref + ndp_theo_ref) / 2

    shared_ions = set(fragments_dict).intersection(set(theo_dict))

    return {
        "ndp_exp_ref": ndp_exp_ref,
        "ndp_theo_ref": ndp_theo_ref,
        "ndp_mean": ndp_mean,
        "n_exp_annotated": len(fragments_dict),
        "n_theo_ions": len(theo_dict),
        "n_shared_ions": len(shared_ions),
        "fraction_theo_found": (
            len(shared_ions) / len(theo_dict)
            if len(theo_dict) > 0
            else 0.0
        ),
    }


# ============================================================
# Process one PSM file
# ============================================================

def process_one_dataset(
    *,
    cell_line: str,
    dataset_name: str,
    psm_file: Path,
    ms2pip_file: Path,
    mzml_dir: Path,
    out_dir: Path,
    tolerance: float,
) -> Path | None:
    print("\n" + "=" * 100)
    print(f"[START] {cell_line} | {dataset_name}")
    print(f"[PSM]    {psm_file}")
    print(f"[MS2PIP] {ms2pip_file}")
    print(f"[MZML]   {mzml_dir}")

    psm_file = Path(psm_file)
    ms2pip_file = Path(ms2pip_file)
    mzml_dir = Path(mzml_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not psm_file.exists():
        print(f"[SKIP] Missing PSM file: {psm_file}")
        return None

    if psm_file.stat().st_size == 0:
        print(f"[SKIP] Empty PSM file: {psm_file}")
        return None

    psms = pd.read_csv(psm_file, sep="\t")

    peptide_col = get_modified_peptide_column(psms)
    charge_col = get_charge_column(psms, psm_file)
    scan_col = get_scan_column(psms, psm_file)

    if "spectrum" not in psms.columns:
        raise ValueError(f"Colonne spectrum absente dans {psm_file}")

    ms2pip_lib = load_ms2pip_library(ms2pip_file)
    mzml_index = build_mzml_index(mzml_dir)

    print(f"[INFO] PSMs: {len(psms)}")
    print(f"[INFO] MS2PIP peptide+charge entries: {len(ms2pip_lib)}")
    print(f"[INFO] mzML indexed: {len(mzml_index)}")

    mzml_cache: dict[str, dict[int, pd.DataFrame]] = {}
    results = []

    for idx, row in psms.iterrows():
        peptide = str(row[peptide_col])
        charge = int(row[charge_col])
        spectrum = str(row["spectrum"])
        scan = int(row[scan_col])

        run_name = extract_run_name_from_spectrum(spectrum)

        base = {
            "cell_line": cell_line,
            "dataset": dataset_name,
            "row_index": idx,
            "spectrum": spectrum,
            "run_name": run_name,
            "scan": scan,
            "modified_peptide": peptide,
            "charge": charge,
            "class": row.get("class", np.nan),
            "protein": row.get("protein", np.nan),
            "peptideprophet_probability": row.get(
                "peptideprophet_probability",
                np.nan,
            ),
            "hyperscore": row.get("hyperscore", np.nan),
            "massdiff": row.get("massdiff", np.nan),
            "num_matched_ions_original": row.get("num_matched_ions", np.nan),
            "tot_num_ions_original": row.get("tot_num_ions", np.nan),
        }

        theo_key = (peptide, charge)

        if theo_key not in ms2pip_lib:
            base.update({
                "status": "missing_ms2pip_prediction",
                "ndp_exp_ref": np.nan,
                "ndp_theo_ref": np.nan,
                "ndp_mean": np.nan,
                "n_exp_annotated": 0,
                "n_theo_ions": 0,
                "n_shared_ions": 0,
                "fraction_theo_found": np.nan,
            })
            results.append(base)
            continue

        if run_name not in mzml_index:
            base.update({
                "status": "missing_mzml_file",
                "ndp_exp_ref": np.nan,
                "ndp_theo_ref": np.nan,
                "ndp_mean": np.nan,
                "n_exp_annotated": 0,
                "n_theo_ions": len(ms2pip_lib[theo_key]),
                "n_shared_ions": 0,
                "fraction_theo_found": 0.0,
            })
            results.append(base)
            continue

        mzml_path = mzml_index[run_name]

        if run_name not in mzml_cache:
            print(f"[LOAD mzML] {cell_line} | {mzml_path.name}")
            mzml_cache[run_name] = load_mzml_ms2_scans(mzml_path)

        scan_dict = mzml_cache[run_name]

        if scan not in scan_dict:
            base.update({
                "status": "missing_scan_in_mzml",
                "ndp_exp_ref": np.nan,
                "ndp_theo_ref": np.nan,
                "ndp_mean": np.nan,
                "n_exp_annotated": 0,
                "n_theo_ions": len(ms2pip_lib[theo_key]),
                "n_shared_ions": 0,
                "fraction_theo_found": 0.0,
            })
            results.append(base)
            continue

        theo_df = ms2pip_lib[theo_key]
        exp_df = scan_dict[scan]

        matched = match_theoretical_to_experimental(
            theo_df=theo_df,
            exp_df=exp_df,
            tolerance=tolerance,
        )

        fragments_dict, theo_dict = build_dicts_from_matched(matched)
        scores = compute_two_way_ndp(fragments_dict, theo_dict)

        base.update({
            "status": "ok",
            "ndp_exp_ref": scores["ndp_exp_ref"],
            "ndp_theo_ref": scores["ndp_theo_ref"],
            "ndp_mean": scores["ndp_mean"],
            "n_exp_annotated": scores["n_exp_annotated"],
            "n_theo_ions": scores["n_theo_ions"],
            "n_shared_ions": scores["n_shared_ions"],
            "fraction_theo_found": scores["fraction_theo_found"],
        })

        results.append(base)

        if (idx + 1) % 500 == 0:
            print(
                f"[PROGRESS] {cell_line} | {dataset_name}: "
                f"{idx + 1}/{len(psms)}"
            )

    out = pd.DataFrame(results)

    out_file = out_dir / f"{cell_line}_{dataset_name}_ndp_scores.tsv"
    out.to_csv(out_file, sep="\t", index=False)

    print(f"[DONE] {cell_line} | {dataset_name}")
    print(out["status"].value_counts(dropna=False))
    print(f"[WRITE] {out_file}")

    return out_file


# ============================================================
# Plotting + stats
# ============================================================

def plot_distribution(
    *,
    df: pd.DataFrame,
    out_png: Path,
    threshold: float,
    title: str | None = None,
    out_stats_tsv: Path | None = None,
) -> None:
    ok = df.loc[
        (df["status"] == "ok") &
        df["ndp_mean"].notna(),
        "ndp_mean",
    ].astype(float)

    if len(ok) == 0:
        print(f"[WARN] No ok ndp_mean values for plot: {out_png}")
        return

    n_total = len(ok)
    n_kept = int((ok >= threshold).sum())
    pct_kept = 100 * n_kept / n_total

    plt.figure(figsize=(8, 4.8))

    plt.hist(
        ok,
        bins=50,
        edgecolor="black",
        alpha=0.85,
    )

    plt.axvline(
        threshold,
        linestyle="--",
        linewidth=2,
        color="black",
    )

    pct_str = f"{pct_kept:.1f}".replace(".", ",")
    n_kept_str = f"{n_kept:,}".replace(",", " ")

    plt.text(
        0.03,
        0.92,
        f"Seuil retenu = {threshold:.2f}\n"
        f"PSM conservés : {n_kept_str} ({pct_str} %)",
        transform=plt.gca().transAxes,
        ha="left",
        va="top",
        fontsize=11,
        bbox=dict(
            boxstyle="round",
            facecolor="white",
            edgecolor="gray",
            alpha=0.9,
        ),
    )

    plt.xlabel("Score moyen de produit scalaire normalisé")
    plt.ylabel("Nombre de PSM")

    if title is not None:
        plt.title(title)

    plt.xlim(0, 1.02)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

    print(f"[WRITE] {out_png}")

    quantiles = ok.quantile([
        0.001,
        0.005,
        0.01,
        0.025,
        0.05,
        0.10,
        0.25,
        0.50,
        0.75,
    ])

    threshold_rows = []

    for thr in [0.70, 0.75, 0.80, 0.85, 0.90]:
        n_below_thr = int((ok < thr).sum())
        n_kept_thr = int((ok >= thr).sum())

        threshold_rows.append({
            "threshold": thr,
            "n_total": n_total,
            "n_below": n_below_thr,
            "pct_below": 100 * n_below_thr / n_total,
            "n_kept": n_kept_thr,
            "pct_kept": 100 * n_kept_thr / n_total,
        })

    threshold_df = pd.DataFrame(threshold_rows)

    print("\n[DESCRIBE]")
    print(ok.describe())

    print("\n[QUANTILES]")
    print(quantiles)

    print("\n[THRESHOLD SUMMARY]")
    print(threshold_df)

    if out_stats_tsv is not None:
        stats_rows = []

        for q, val in quantiles.items():
            stats_rows.append({
                "stat_type": "quantile",
                "metric": f"q{q}",
                "threshold": np.nan,
                "value": val,
                "n_total": n_total,
                "n_below": np.nan,
                "pct_below": np.nan,
                "n_kept": np.nan,
                "pct_kept": np.nan,
            })

        for _, r in threshold_df.iterrows():
            stats_rows.append({
                "stat_type": "threshold",
                "metric": f"threshold_{r['threshold']}",
                "threshold": r["threshold"],
                "value": np.nan,
                "n_total": r["n_total"],
                "n_below": r["n_below"],
                "pct_below": r["pct_below"],
                "n_kept": r["n_kept"],
                "pct_kept": r["pct_kept"],
            })

        stats_df = pd.DataFrame(stats_rows)
        stats_df.to_csv(out_stats_tsv, sep="\t", index=False)
        print(f"[WRITE] {out_stats_tsv}")


# ============================================================
# Main pipeline function
# ============================================================

def run_ms2pip_experimental_scoring(
    config: dict,
) -> dict[str, dict[str, Path]]:
    """
    Lance l'étape 5 pour toutes les lignées dans config.yaml.

    Entrées :
        analysis_outdir/{cell}/02_IL_filtering/
        ms2pip_output_root/{cell}/
        mzml_root/{cell}/

    Sorties :
        analysis_outdir/{cell}/05_ms2pip_experimental_scoring/
    """
    analysis_outdir = get_analysis_outdir(config)
    ms2pip_output_root = get_ms2pip_output_root(config)
    mzml_root = get_mzml_root(config)

    scoring_cfg = config.get("ms2pip_scoring", {})
    tolerance = float(scoring_cfg.get("tolerance", 0.02))
    threshold = float(scoring_cfg.get("ndp_threshold", 0.70))

    print("\n[CONFIG] Chemins utilisés par l'étape 5")
    print(f"[CONFIG] analysis_outdir     : {analysis_outdir}")
    print(f"[CONFIG] ms2pip_output_root  : {ms2pip_output_root}")
    print(f"[CONFIG] mzml_root           : {mzml_root}")
    print(f"[CONFIG] tolerance           : {tolerance}")
    print(f"[CONFIG] ndp_threshold       : {threshold}")

    outputs: dict[str, dict[str, Path]] = {}

    for cell_line, _runs in config["samples"].items():
        print("\n" + "#" * 100)
        print(f"[Étape 5] Scoring expérimental MS2PIP — {cell_line}")
        print("#" * 100)

        input_dir = analysis_outdir / cell_line / "02_IL_filtering"
        ms2pip_dir = get_ms2pip_dir(config, cell_line)
        mzml_dir = mzml_root / cell_line
        out_dir = analysis_outdir / cell_line / "05_ms2pip_experimental_scoring"
        out_dir.mkdir(parents=True, exist_ok=True)

        print(f"[PATH] PSM input dir : {input_dir}")
        print(f"[PATH] MS2PIP dir    : {ms2pip_dir}")
        print(f"[PATH] mzML dir      : {mzml_dir}")
        print(f"[PATH] Output dir    : {out_dir}")

        outputs[cell_line] = {}

        if not input_dir.exists():
            print(f"[SKIP] Missing input dir: {input_dir}")
            continue

        if not ms2pip_dir.exists():
            print(f"[SKIP] Missing MS2PIP dir: {ms2pip_dir}")
            continue

        if not mzml_dir.exists():
            print(f"[SKIP] Missing mzML dir: {mzml_dir}")
            continue

        for dataset_name, cfg in FILE_PAIRS.items():
            psm_file = input_dir / cfg["psm"]

            try:
                ms2pip_file = find_existing_ms2pip_file(
                    ms2pip_dir=ms2pip_dir,
                    candidates=cfg["ms2pip_candidates"],
                )
            except Exception as e:
                print(f"[SKIP] {cell_line} | {dataset_name}: {e}")
                continue

            try:
                out_file = process_one_dataset(
                    cell_line=cell_line,
                    dataset_name=dataset_name,
                    psm_file=psm_file,
                    ms2pip_file=ms2pip_file,
                    mzml_dir=mzml_dir,
                    out_dir=out_dir,
                    tolerance=tolerance,
                )
            except Exception as e:
                print(f"[ERROR] {cell_line} | {dataset_name}: {e}")
                continue

            if out_file is None:
                continue

            outputs[cell_line][dataset_name] = out_file

            df_one = pd.read_csv(out_file, sep="\t")

            plot_distribution(
                df=df_one,
                out_png=out_dir / f"{cell_line}_{dataset_name}_ndp_mean_distribution.png",
                title=f"NDP mean distribution | {cell_line} | {dataset_name}",
                threshold=threshold,
                out_stats_tsv=out_dir / f"{cell_line}_{dataset_name}_ndp_mean_stats.tsv",
            )

    if not any(outputs[cell] for cell in outputs):
        raise RuntimeError("Aucun fichier NDP produit par l'étape 5.")

    print("\n[DONE] Étape 5 — scoring expérimental MS2PIP terminé.")

    return outputs


if __name__ == "__main__":
    raise SystemExit(
        "Ce module est conçu pour être appelé depuis main.py.\n"
        "Utilise plutôt : python main.py -c config.yaml"
    )
