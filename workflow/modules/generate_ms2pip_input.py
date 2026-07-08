#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Étape 4 — Génération des fichiers PEPREC pour MS2PIP.

Cette étape prend les fichiers propres produits par l'étape 2 :

    analysis_outdir/{cell}/02_IL_filtering/
        canonical.tsv
        canonical_variant.tsv
        nc_wt.tsv
        nc_variant.tsv

et produit les fichiers PEPREC utilisés par MS2PIP :

    analysis_outdir/{cell}/04_ms2pip_input/
        canonical.peprec
        canonical_variant.peprec
        nc_wt.peprec
        nc_variant.peprec
        ms2pip_config.txt
        unknown_modifications.tsv
        peprec_generation_summary.tsv

Important :
- aucun dossier global analysis_outdir/04_ms2pip_input n'est créé;
- un seul fichier unknown_modifications.tsv est produit par lignée;
- la colonne "group" indique canonical, canonical_variant, nc_wt ou nc_variant.

Ce module est conçu pour être appelé depuis main.py :

    from modules.generate_ms2pip_input import run_generate_ms2pip_input
    peprec_outputs = run_generate_ms2pip_input(config)
"""

from __future__ import annotations

from pathlib import Path
import ast
import re
from typing import Any

import pandas as pd


GROUPS = ["canonical", "canonical_variant", "nc_wt", "nc_variant"]

PEPREC_COLUMNS = ["spec_id", "modifications", "peptide", "charge"]

UNKNOWN_MOD_COLUMNS = [
    "cell",
    "group",
    "original_row_index",
    "peptide",
    "original_peptide",
    "modified_peptide",
    "modifications_original",
    "assumed_charge",
    "charge",
    "spectrum",
    "protein",
]


# ============================================================
# Helpers
# ============================================================

def parse_list_like(x: Any) -> list[Any]:
    """
    Convertit une chaîne de type liste Python en vraie liste.

    Exemples :
        []
        ['160.031@3']
        ['160.031@3', '147.035@10']
    """
    if pd.isna(x):
        return []

    if isinstance(x, list):
        return x

    s = str(x).strip()

    if s in ["", "None", "nan", "[]"]:
        return []

    try:
        value = ast.literal_eval(s)
        if isinstance(value, list):
            return value
        return [value]
    except Exception:
        return [s]


def clean_peptide_sequence(peptide: Any) -> str | None:
    """
    Nettoie la séquence peptidique pour garder seulement la séquence non modifiée.

    Exemples :
        A.PEPTIDE.K      -> PEPTIDE
        PEPTIDE          -> PEPTIDE
        M[15.9949]PEP    -> MPEP
        M(15.9949)PEP    -> MPEP
    """
    if pd.isna(peptide):
        return None

    s = str(peptide).strip()

    parts = s.split(".")
    if len(parts) == 3:
        s = parts[1]

    s = re.sub(r"\[.*?\]", "", s)
    s = re.sub(r"\(.*?\)", "", s)
    s = re.sub(r"[^A-Z]", "", s)

    return s if s else None


def valid_ms2pip_peptide(seq: str | None) -> bool:
    """
    Filtre basique de compatibilité MS2PIP.
    """
    if seq is None:
        return False

    if len(seq) <= 2:
        return False

    if len(seq) >= 100:
        return False

    forbidden = set("BJOUXZ")
    if any(aa in forbidden for aa in seq):
        return False

    return True


def parse_charge(x: Any) -> int | None:
    """
    Convertit la charge en int.
    """
    if pd.isna(x):
        return None

    try:
        return int(float(str(x).replace("+", "").strip()))
    except Exception:
        return None


def get_charge_value(row: pd.Series) -> Any:
    """
    Récupère la charge depuis assumed_charge ou charge.
    """
    if "assumed_charge" in row.index:
        return row.get("assumed_charge")

    if "charge" in row.index:
        return row.get("charge")

    return None


def msfragger_mods_to_ms2pip(
    modifications: Any,
    peptide: str | None = None,
) -> str | None:
    """
    Convertit les modifications MSFragger/Philosopher vers le format PEPREC MS2PIP.

    Entrées attendues :
        []
        ['160.031@3']
        ['160.031@3', '147.035@10']

    Sorties :
        -
        3|Carbamidomethyl
        3|Carbamidomethyl|10|Oxidation

    Conversions supportées :
        C avec masse ~160.031 -> Carbamidomethyl
        M avec masse ~147.035 -> Oxidation
        N-term avec masse ~43.018 à pos 0 -> Acetyl

    Retourne None si une modification est inconnue.
    """
    mods = parse_list_like(modifications)

    if len(mods) == 0:
        return "-"

    if peptide is None:
        return None

    ms2pip_mods: list[str] = []

    for mod in mods:
        mod = str(mod).strip()

        if "@" not in mod:
            return None

        mass_str, pos_str = mod.split("@", 1)

        try:
            mass = float(mass_str)
            pos = int(pos_str)
        except Exception:
            return None

        if pos == 0:
            # N-terminal acetylation selon la convention d'export.
            if abs(mass - 43.018) < 0.02:
                mod_name = "Acetyl"
            else:
                return None

        elif 1 <= pos <= len(peptide):
            aa = peptide[pos - 1]

            # C[160.031] = C + 57.021 ≈ Carbamidomethyl.
            if aa == "C" and abs(mass - 160.031) < 0.02:
                mod_name = "Carbamidomethyl"

            # M[147.035] = M + 15.995 ≈ Oxidation.
            elif aa == "M" and abs(mass - 147.035) < 0.02:
                mod_name = "Oxidation"

            else:
                return None

        else:
            return None

        ms2pip_mods.extend([str(pos), mod_name])

    return "|".join(ms2pip_mods) if ms2pip_mods else "-"


def make_spec_id(peptide: str, charge: int, modifications: str) -> str:
    """
    Crée un spec_id stable.

    La modification est incluse parce qu'un même peptide/charge peut exister
    sous différentes formes modifiées.
    """
    mod_tag = modifications.replace("|", "_").replace("-", "unmod")
    return f"{peptide}_{charge}_{mod_tag}"


def write_ms2pip_config(outdir: Path) -> Path:
    """
    Écrit un fichier de configuration minimal pour MS2PIP.
    """
    config_path = outdir / "ms2pip_config.txt"

    text = "\n".join(
        [
            "model=HCD",
            "frag_error=0.02",
            "out=csv",
            "ptm=Carbamidomethyl,57.021464,opt,C",
            "ptm=Oxidation,15.994915,opt,M",
            "ptm=Acetyl,42.010565,opt,N-term",
            "",
        ]
    )

    config_path.write_text(text)

    return config_path


# ============================================================
# Processing
# ============================================================

def generate_peprec(
    input_path: Path,
    output_path: Path,
    *,
    cell: str,
    group: str,
    keep_unknown_mods_as_unmodified: bool = False,
) -> tuple[dict[str, Any] | None, pd.DataFrame]:
    """
    Génère un fichier PEPREC à partir d'un TSV final.

    Retourne :
        result:
            résumé de génération pour ce fichier

        unknown_mods_df:
            DataFrame des modifications inconnues pour ce groupe.
            Il contient déjà les colonnes cell et group.
    """
    if not input_path.exists():
        print(f"[SKIP] Fichier introuvable : {input_path}")
        return None, pd.DataFrame(columns=UNKNOWN_MOD_COLUMNS)

    print("\n" + "=" * 100)
    print(f"[INPUT]  {input_path}")
    print(f"[OUTPUT] {output_path}")

    df = pd.read_csv(input_path, sep="\t")
    print(f"[READ] shape={df.shape}")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    if df.empty:
        pd.DataFrame(columns=PEPREC_COLUMNS).to_csv(output_path, sep=" ", index=False)

        result = {
            "cell": cell,
            "group": group,
            "input": str(input_path),
            "output": str(output_path),
            "n_input_rows": 0,
            "n_invalid_rows": 0,
            "n_unknown_mod_rows": 0,
            "n_peprec_before_dedup": 0,
            "n_dedup_removed": 0,
            "n_peprec_rows": 0,
            "status": "empty",
        }

        return result, pd.DataFrame(columns=UNKNOWN_MOD_COLUMNS)

    if "peptide" not in df.columns:
        raise ValueError(
            f"Colonne 'peptide' absente dans {input_path}. "
            f"Colonnes disponibles : {df.columns.tolist()}"
        )

    if "assumed_charge" not in df.columns and "charge" not in df.columns:
        raise ValueError(
            f"Aucune colonne de charge trouvée dans {input_path}. "
            "Colonnes attendues : 'assumed_charge' ou 'charge'. "
            f"Colonnes disponibles : {df.columns.tolist()}"
        )

    rows_peprec: list[dict[str, Any]] = []
    rows_unknown_mods: list[dict[str, Any]] = []

    n_invalid = 0
    n_unknown_mods = 0

    for idx, row in df.iterrows():
        peptide = clean_peptide_sequence(row.get("peptide", None))
        charge = parse_charge(get_charge_value(row))

        if not valid_ms2pip_peptide(peptide):
            n_invalid += 1
            continue

        if charge is None:
            n_invalid += 1
            continue

        modifications = msfragger_mods_to_ms2pip(
            row.get("modifications", "[]"),
            peptide=peptide,
        )

        if modifications is None:
            n_unknown_mods += 1

            rows_unknown_mods.append(
                {
                    "cell": cell,
                    "group": group,
                    "original_row_index": idx,
                    "peptide": peptide,
                    "original_peptide": row.get("peptide", pd.NA),
                    "modified_peptide": row.get("modified_peptide", pd.NA),
                    "modifications_original": row.get("modifications", pd.NA),
                    "assumed_charge": row.get("assumed_charge", pd.NA),
                    "charge": row.get("charge", pd.NA),
                    "spectrum": row.get("spectrum", pd.NA),
                    "protein": row.get("protein", pd.NA),
                }
            )

            if keep_unknown_mods_as_unmodified:
                modifications = "-"
            else:
                continue

        spec_id = make_spec_id(peptide, charge, modifications)

        rows_peprec.append(
            {
                "spec_id": spec_id,
                "modifications": modifications,
                "peptide": peptide,
                "charge": charge,
            }
        )

    peprec_df = pd.DataFrame(rows_peprec, columns=PEPREC_COLUMNS)
    unknown_mods_df = pd.DataFrame(rows_unknown_mods, columns=UNKNOWN_MOD_COLUMNS)

    before_dedup = len(peprec_df)

    # MS2PIP prédit la même chose pour peptide + charge + modification identiques.
    peprec_df = peprec_df.drop_duplicates(
        subset=["modifications", "peptide", "charge"]
    ).copy()

    n_dedup_removed = before_dedup - len(peprec_df)

    peprec_df.to_csv(output_path, sep=" ", index=False)

    print(f"[COUNT] input_rows={len(df):,}")
    print(f"[COUNT] invalid_rows={n_invalid:,}")
    print(f"[COUNT] unknown_mod_rows={n_unknown_mods:,}")
    print(f"[COUNT] peprec_before_dedup={before_dedup:,}")
    print(f"[COUNT] dedup_removed={n_dedup_removed:,}")
    print(f"[COUNT] peprec_final={len(peprec_df):,}")

    result = {
        "cell": cell,
        "group": group,
        "input": str(input_path),
        "output": str(output_path),
        "n_input_rows": len(df),
        "n_invalid_rows": n_invalid,
        "n_unknown_mod_rows": n_unknown_mods,
        "n_peprec_before_dedup": before_dedup,
        "n_dedup_removed": n_dedup_removed,
        "n_peprec_rows": len(peprec_df),
        "status": "ok",
    }

    return result, unknown_mods_df


def run_generate_ms2pip_input(config: dict) -> dict[str, dict[str, Path]]:
    """
    Génère les fichiers PEPREC pour toutes les lignées définies dans config.yaml.

    Paramètres optionnels dans config.yaml :

    ms2pip:
      keep_unknown_mods_as_unmodified: false
    """
    analysis_outdir = Path(config["path"]["analysis_outdir"])

    keep_unknown_mods_as_unmodified = config.get("ms2pip", {}).get(
        "keep_unknown_mods_as_unmodified",
        False,
    )

    outputs: dict[str, dict[str, Path]] = {}

    for cell in config["samples"].keys():
        print("\n" + "=" * 100)
        print(f"[Étape 4] Génération PEPREC MS2PIP — {cell}")
        print("=" * 100)

        input_dir = analysis_outdir / cell / "02_IL_filtering"
        cell_outdir = analysis_outdir / cell / "04_ms2pip_input"
        cell_outdir.mkdir(parents=True, exist_ok=True)

        config_path = write_ms2pip_config(cell_outdir)
        print(f"[CONFIG] {config_path}")

        outputs[cell] = {}

        cell_summary_rows: list[dict[str, Any]] = []
        cell_unknown_mods: list[pd.DataFrame] = []

        for group in GROUPS:
            input_path = input_dir / f"{group}.tsv"
            output_path = cell_outdir / f"{group}.peprec"

            result, unknown_mods_df = generate_peprec(
                input_path=input_path,
                output_path=output_path,
                cell=cell,
                group=group,
                keep_unknown_mods_as_unmodified=keep_unknown_mods_as_unmodified,
            )

            if not unknown_mods_df.empty:
                cell_unknown_mods.append(unknown_mods_df)

            if result is None:
                continue

            cell_summary_rows.append(result)
            outputs[cell][group] = Path(result["output"])

        # Un seul fichier unknown_modifications.tsv par lignée.
        unknown_out = cell_outdir / "unknown_modifications.tsv"

        if cell_unknown_mods:
            unknown_df = pd.concat(cell_unknown_mods, ignore_index=True)
        else:
            unknown_df = pd.DataFrame(columns=UNKNOWN_MOD_COLUMNS)

        unknown_df.to_csv(unknown_out, sep="\t", index=False)
        print(f"[OK] Unknown modifications : {unknown_out} {unknown_df.shape}")

        outputs[cell]["unknown_modifications"] = unknown_out

        # Un seul résumé par lignée.
        summary_path = cell_outdir / "peprec_generation_summary.tsv"

        if cell_summary_rows:
            summary_df = pd.DataFrame(cell_summary_rows)
            summary_df.to_csv(summary_path, sep="\t", index=False)

            print("\n" + "=" * 100)
            print(f"[SUMMARY] {summary_path}")
            print(
                summary_df[
                    [
                        "cell",
                        "group",
                        "n_input_rows",
                        "n_invalid_rows",
                        "n_unknown_mod_rows",
                        "n_dedup_removed",
                        "n_peprec_rows",
                        "status",
                    ]
                ].to_string(index=False)
            )
        else:
            summary_df = pd.DataFrame(
                columns=[
                    "cell",
                    "group",
                    "input",
                    "output",
                    "n_input_rows",
                    "n_invalid_rows",
                    "n_unknown_mod_rows",
                    "n_peprec_before_dedup",
                    "n_dedup_removed",
                    "n_peprec_rows",
                    "status",
                ]
            )
            summary_df.to_csv(summary_path, sep="\t", index=False)
            print(f"[WARNING] Aucun PEPREC produit pour {cell}. Summary vide : {summary_path}")

        outputs[cell]["summary"] = summary_path

    print("\n[DONE] Génération PEPREC terminée.")

    return outputs


if __name__ == "__main__":
    raise SystemExit(
        "Ce module est conçu pour être appelé depuis main.py.\n"
        "Utilise plutôt : python main.py -c config.yaml"
    )
