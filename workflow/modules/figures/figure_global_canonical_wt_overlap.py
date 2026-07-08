#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Figure globale — Recouvrement des peptides/protéines/gènes canoniques WT détectés.

Cette figure compare les lignées présentes dans config["samples"].

Entrées attendues :
    analysis_outdir/{sample}/02_IL_filtering/canonical.tsv

    ../data/fasta/{sample}/{sample}_contam.fasta

Optionnel :
    path:
      test_proteome_fasta: /chemin/vers/proteome.fasta

Si test_proteome_fasta est fourni, le script génère aussi les figures après filtre
sur les ENSP présents dans ce protéome.

Sorties :
    outdir/canonical_wt_overlap_levels/
        ALL_canonical_WT_peptide_to_protein_to_gene_long.tsv
        canonical_WT_counts_by_level_per_sample.tsv
        canonical_WT_peptide_level_venn.png/svg/pdf
        canonical_WT_protein_accession_level_venn.png/svg/pdf
        canonical_WT_gene_level_venn.png/svg/pdf
        ...
"""

from pathlib import Path
import ast
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    from matplotlib_venn import venn2
    HAS_VENN = True
except ImportError:
    HAS_VENN = False


# ============================================================
# CONSTANTES
# ============================================================

MIN_PEPTIDE_LENGTH = 9

COLOR_LEFT = "#F2B6B6"
COLOR_RIGHT = "#B7C9DD"
COLOR_SHARED = "#D9D9D9"
EDGE_COLOR = "#4D4D4D"


# ============================================================
# PATH HELPERS
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


def get_fasta_root(config: dict) -> Path:
    """
    Racine des FASTA custom.

    Défaut :
        ../data/fasta
    """
    path_cfg = config.get("path", {})
    return Path(path_cfg.get("fasta_root", "../data/fasta"))


def get_reference_fasta(config: dict, sample: str) -> Path:
    """
    FASTA attendue :
        ../data/fasta/{sample}/{sample}_contam.fasta
    """
    fasta_root = get_fasta_root(config)
    return fasta_root / sample / f"{sample}_contam.fasta"


def get_canonical_tsv(config: dict, sample: str) -> Path:
    """
    Fichier canonical.tsv après I/L filtering.

    Attendu :
        analysis_outdir/{sample}/02_IL_filtering/canonical.tsv
    """
    analysis_outdir = get_analysis_outdir(config)
    return analysis_outdir / sample / "02_IL_filtering" / "canonical.tsv"


def get_test_proteome_fasta(config: dict) -> Path | None:
    """
    Protéome optionnel utilisé comme whitelist ENSP.

    Exemple config :
        path:
          test_proteome_fasta: "../analysis/5_figures/test/proteome.fasta"

    Si absent, les figures filtrées sont ignorées.
    """
    path_cfg = config.get("path", {})
    p = path_cfg.get("test_proteome_fasta")

    if not p:
        return None

    return Path(p)


def get_config_samples(config: dict) -> list[str]:
    if "samples" not in config:
        raise KeyError("Section 'samples' absente du config.yaml.")

    samples = list(config["samples"].keys())

    if len(samples) != 2:
        raise ValueError(
            "Cette figure globale est conçue pour exactement deux samples.\n"
            f"Samples reçus : {samples}"
        )

    return samples


def sanitize_name(x: str) -> str:
    """
    Nom safe pour fichiers.
    """
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(x))


# ============================================================
# BASIC HELPERS
# ============================================================

def fmt_int_space(n):
    return f"{int(n):,}".replace(",", " ")


def parse_list_like(x):
    """
    Parser robuste pour des colonnes du type :
        ['ENSP...', 'ENSP...']
        ENSP...
        ENSP...;ENSP...
    """
    if isinstance(x, list):
        return x

    if x is None:
        return []

    try:
        if pd.isna(x):
            return []
    except Exception:
        pass

    x = str(x).strip()

    if not x or x.lower() == "nan":
        return []

    try:
        parsed = ast.literal_eval(x)

        if isinstance(parsed, list):
            return parsed

        if isinstance(parsed, tuple):
            return list(parsed)

        return [parsed]

    except Exception:
        if ";" in x:
            return [v.strip() for v in x.split(";") if v.strip()]

        if "," in x:
            return [v.strip() for v in x.split(",") if v.strip()]

        return [x]


def strip_version(identifier):
    identifier = str(identifier).strip().strip("'").strip('"')
    return identifier.split(".")[0]


def get_peptide(row):
    peptide = str(row.get("peptide", "")).strip()

    if not peptide or peptide.lower() == "nan":
        return None

    return peptide


def clean_canonical_ensp_accessions(protein_field):
    """
    Extrait les accessions ENSP canoniques depuis la colonne protein.

    Garde :
        ENSP00000...

    Exclut :
        - entrées mutées
        - headers custom contenant |
        - contaminants
        - entrées non ENSP
    """
    proteins = parse_list_like(protein_field)
    cleaned = []

    for p in proteins:
        if p is None:
            continue

        try:
            if pd.isna(p):
                continue
        except Exception:
            pass

        p = str(p).strip().strip("'").strip('"')

        if not p:
            continue

        if "_mut" in p:
            continue

        if "|" in p:
            continue

        if not re.fullmatch(r"ENSP\d+(?:\.\d+)?", p):
            continue

        cleaned.append(strip_version(p))

    return sorted(set(cleaned))


# ============================================================
# FASTA MAPPING
# ============================================================

def extract_ids_and_gene_from_fasta_header(header):
    """
    Parse les headers FASTA Ensembl.

    Exemple :
        >ENSP000003... pep chromosome:... gene:ENSG000001...
        transcript:ENST00000... gene_symbol:TP53 ...

    Retourne :
        ids, gene_key, gene_symbol, ensg

    gene_key :
        ENSG si disponible, sinon gene_symbol.
    """
    header = header.strip()

    ensp = None
    enst = None
    ensg = None
    gene_symbol = None

    m = re.search(r"(ENSP\d+(?:\.\d+)?)", header)
    if m:
        ensp = strip_version(m.group(1))

    m = re.search(r"transcript:(ENST\d+(?:\.\d+)?)", header)
    if m:
        enst = strip_version(m.group(1))

    m = re.search(r"gene:(ENSG\d+(?:\.\d+)?)", header)
    if m:
        ensg = strip_version(m.group(1))

    m = re.search(r"gene_symbol:([^\s]+)", header)
    if m:
        gene_symbol = m.group(1).strip()

    gene_key = ensg if ensg else gene_symbol

    ids = []
    for x in [ensp, enst, ensg]:
        if x:
            ids.append(x)

    return ids, gene_key, gene_symbol, ensg


def build_identifier_to_gene_mapping(reference_fastas: dict[str, Path], extra_fastas=None):
    """
    Construit :
        ENSP/ENST/ENSG -> gene_key
        ENSP/ENST/ENSG -> gene_symbol
        ENSP/ENST/ENSG -> ENSG
    """
    if extra_fastas is None:
        extra_fastas = []

    fasta_paths = list(reference_fastas.values()) + list(extra_fastas)

    id_to_gene_key = {}
    id_to_gene_symbol = {}
    id_to_ensg = {}

    for fasta_path in fasta_paths:
        fasta_path = Path(fasta_path)

        if not fasta_path.exists():
            print(f"[WARN] FASTA introuvable pour mapping : {fasta_path}")
            continue

        print(f"[LOAD] FASTA mapping : {fasta_path}")

        with open(fasta_path, "r") as handle:
            for line in handle:
                if not line.startswith(">"):
                    continue

                header = line[1:].strip()
                ids, gene_key, gene_symbol, ensg = extract_ids_and_gene_from_fasta_header(header)

                if not ids or not gene_key:
                    continue

                for identifier in ids:
                    id_to_gene_key[identifier] = gene_key

                    if gene_symbol:
                        id_to_gene_symbol[identifier] = gene_symbol

                    if ensg:
                        id_to_ensg[identifier] = ensg

    print(f"[INFO] Identifiants mappés à un gene_key : {len(id_to_gene_key)}")

    return id_to_gene_key, id_to_gene_symbol, id_to_ensg


def read_ensp_whitelist_from_fasta(path: Path):
    """
    Lit un protéome FASTA et retourne les ENSP présents.
    """
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"test_proteome_fasta introuvable : {path}")

    ensp_set = set()

    with open(path, "r") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue

            header = line[1:].strip()
            matches = re.findall(r"ENSP\d+(?:\.\d+)?", header)

            for m in matches:
                ensp_set.add(strip_version(m))

    return ensp_set


# ============================================================
# BUILD DETECTION TABLES
# ============================================================

def load_canonical_wt_detections(
    canonical_tsv: Path,
    sample: str,
    id_to_gene_key: dict,
    id_to_gene_symbol: dict,
    id_to_ensg: dict,
    allowed_ensp=None,
):
    """
    Construit une table longue peptide -> ENSP -> gène depuis canonical.tsv.

    Pour niveau peptide :
        set(peptide)

    Pour niveau protéine :
        set(protein_id)

    Pour niveau gène :
        set(gene_key)

    Si allowed_ensp est fourni :
        la liste ENSP est filtrée avant construction.
    """
    canonical_tsv = Path(canonical_tsv)

    if not canonical_tsv.exists():
        raise FileNotFoundError(f"canonical.tsv introuvable : {canonical_tsv}")

    df = pd.read_csv(canonical_tsv, sep="\t")

    required = {"protein", "peptide"}
    missing = required - set(df.columns)

    if missing:
        raise ValueError(f"Colonnes manquantes dans {canonical_tsv}: {missing}")

    rows = []
    unmapped_rows = []

    for idx, row in df.iterrows():
        peptide = get_peptide(row)

        if peptide is None:
            continue

        if len(peptide) < MIN_PEPTIDE_LENGTH:
            continue

        ensp_list = clean_canonical_ensp_accessions(row["protein"])

        if allowed_ensp is not None:
            ensp_list = sorted(set(p for p in ensp_list if p in allowed_ensp))

        if len(ensp_list) == 0:
            continue

        for ensp in ensp_list:
            gene_key = id_to_gene_key.get(ensp, f"UNMAPPED_{ensp}")
            gene_symbol = id_to_gene_symbol.get(ensp, "")
            ensg = id_to_ensg.get(ensp, "")

            if gene_key.startswith("UNMAPPED_"):
                unmapped_rows.append({
                    "sample": sample,
                    "cell_line": sample,
                    "row_index": idx,
                    "peptide": peptide,
                    "protein_id": ensp,
                })

            rows.append({
                "sample": sample,
                "cell_line": sample,
                "row_index": idx,
                "peptide": peptide,
                "protein_id": ensp,
                "gene_key": gene_key,
                "gene_symbol": gene_symbol,
                "ensg": ensg,
                "spectrum": row.get("spectrum", np.nan),
                "run": row.get("run", np.nan),
                "input_sample": row.get("sample", np.nan),
                "hyperscore": row.get("hyperscore", np.nan),
                "peptideprophet_probability": row.get("peptideprophet_probability", np.nan),
            })

    long_df = pd.DataFrame(rows)
    unmapped_df = pd.DataFrame(unmapped_rows)

    if long_df.empty:
        return long_df, unmapped_df

    long_df = long_df.drop_duplicates(
        subset=["sample", "peptide", "protein_id", "gene_key"]
    ).reset_index(drop=True)

    return long_df, unmapped_df


def summarize_counts(long_df: pd.DataFrame, samples: list[str], label: str):
    """
    Résumé par sample aux niveaux peptide/protéine/gène.
    """
    rows = []

    for sample in samples:
        sub = long_df[long_df["sample"] == sample].copy()

        rows.append({
            "analysis": label,
            "sample": sample,
            "n_peptides": sub["peptide"].nunique() if not sub.empty else 0,
            "n_protein_accessions_ENSP": sub["protein_id"].nunique() if not sub.empty else 0,
            "n_genes": sub["gene_key"].nunique() if not sub.empty else 0,
        })

    return pd.DataFrame(rows)


# ============================================================
# VENN PLOTS AND TABLE EXPORTS
# ============================================================

def style_venn_patches(v):
    patch_specs = {
        "10": (COLOR_LEFT, 0.85),
        "01": (COLOR_RIGHT, 0.85),
        "11": (COLOR_SHARED, 0.90),
    }

    for region_id, (color, alpha) in patch_specs.items():
        patch = v.get_patch_by_id(region_id)

        if patch is not None:
            patch.set_color(color)
            patch.set_alpha(alpha)
            patch.set_edgecolor(EDGE_COLOR)
            patch.set_linewidth(2.5)


def export_overlap_sets(
    long_df: pd.DataFrame,
    samples: list[str],
    id_col: str,
    outdir: Path,
    outprefix: str,
    level_label: str,
):
    """
    Exporte les listes left/shared/right et les tables contextuelles.
    """
    left_sample = samples[0]
    right_sample = samples[1]

    left_set = set(
        long_df.loc[long_df["sample"] == left_sample, id_col]
        .dropna()
        .astype(str)
    )

    right_set = set(
        long_df.loc[long_df["sample"] == right_sample, id_col]
        .dropna()
        .astype(str)
    )

    left_only = left_set - right_set
    right_only = right_set - left_set
    shared = left_set & right_set
    union = left_set | right_set

    left_safe = sanitize_name(left_sample)
    right_safe = sanitize_name(right_sample)

    summary = pd.DataFrame([{
        "level": level_label,
        "id_column": id_col,
        f"{left_sample}_specific": len(left_only),
        "shared": len(shared),
        f"{right_sample}_specific": len(right_only),
        "total_union": len(union),
        f"{left_sample}_total": len(left_set),
        f"{right_sample}_total": len(right_set),
        "pct_shared_over_union": 100 * len(shared) / len(union) if union else 0,
        f"pct_shared_over_{left_sample}": 100 * len(shared) / len(left_set) if left_set else 0,
        f"pct_shared_over_{right_sample}": 100 * len(shared) / len(right_set) if right_set else 0,
    }])

    summary_out = outdir / f"{outprefix}_summary.tsv"
    summary.to_csv(summary_out, sep="\t", index=False)

    pd.DataFrame({id_col: sorted(left_only)}).to_csv(
        outdir / f"{outprefix}_{left_safe}_specific_ids.tsv",
        sep="\t",
        index=False,
    )

    pd.DataFrame({id_col: sorted(shared)}).to_csv(
        outdir / f"{outprefix}_shared_ids.tsv",
        sep="\t",
        index=False,
    )

    pd.DataFrame({id_col: sorted(right_only)}).to_csv(
        outdir / f"{outprefix}_{right_safe}_specific_ids.tsv",
        sep="\t",
        index=False,
    )

    for label, ids in [
        (f"{left_safe}_specific", left_only),
        ("shared", shared),
        (f"{right_safe}_specific", right_only),
    ]:
        sub = long_df[long_df[id_col].astype(str).isin(ids)].copy()

        if not sub.empty:
            context = (
                sub.groupby(["sample", id_col], as_index=False)
                .agg({
                    "peptide": lambda x: ";".join(sorted(set(str(v) for v in x if str(v) != "nan"))),
                    "protein_id": lambda x: ";".join(sorted(set(str(v) for v in x if str(v) != "nan"))),
                    "gene_key": lambda x: ";".join(sorted(set(str(v) for v in x if str(v) != "nan"))),
                    "gene_symbol": lambda x: ";".join(sorted(set(str(v) for v in x if str(v) not in {"nan", ""}))),
                    "ensg": lambda x: ";".join(sorted(set(str(v) for v in x if str(v) not in {"nan", ""}))),
                })
            )

            context["n_peptides"] = context["peptide"].apply(
                lambda x: len([v for v in str(x).split(";") if v])
            )

            context["n_protein_accessions"] = context["protein_id"].apply(
                lambda x: len([v for v in str(x).split(";") if v])
            )

            context["n_genes"] = context["gene_key"].apply(
                lambda x: len([v for v in str(x).split(";") if v])
            )
        else:
            context = pd.DataFrame()

        context.to_csv(
            outdir / f"{outprefix}_{label}_context.tsv",
            sep="\t",
            index=False,
        )

    print(f"[WRITE] {summary_out}")

    return left_set, right_set, summary


def make_venn(
    long_df: pd.DataFrame,
    samples: list[str],
    id_col: str,
    title: str,
    outdir: Path,
    outprefix: str,
    level_label: str,
):
    """
    Exporte les tables et génère un Venn.
    """
    left_set, right_set, summary = export_overlap_sets(
        long_df=long_df,
        samples=samples,
        id_col=id_col,
        outdir=outdir,
        outprefix=outprefix,
        level_label=level_label,
    )

    if not HAS_VENN:
        print("[WARN] matplotlib-venn non installé. Venn ignoré.")
        return summary

    left_sample = samples[0]
    right_sample = samples[1]

    left_only = len(left_set - right_set)
    right_only = len(right_set - left_set)
    shared = len(left_set & right_set)
    total = len(left_set | right_set)

    fig, ax = plt.subplots(figsize=(9, 7))

    v = venn2(
        subsets=(left_only, right_only, shared),
        set_labels=(left_sample, right_sample),
        ax=ax,
    )

    style_venn_patches(v)

    region_counts = {
        "10": left_only,
        "01": right_only,
        "11": shared,
    }

    for region_id, n in region_counts.items():
        label = v.get_label_by_id(region_id)

        if label is not None:
            if id_col == "gene_key":
                label.set_text(fmt_int_space(n))
            else:
                pct = 100 * n / total if total > 0 else 0
                label.set_text(f"{fmt_int_space(n)}\n({pct:.1f}%)")

            label.set_fontsize(16)
            label.set_color("black")
            label.set_fontweight("bold")

    for label in v.set_labels:
        if label is not None:
            label.set_fontsize(22)
            label.set_fontweight("bold")

    ax.set_title(title, fontsize=17, fontweight="bold", pad=20)

    plt.tight_layout()

    for ext in ["png", "svg", "pdf"]:
        out = outdir / f"{outprefix}.{ext}"
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(f"[WRITE] {out}")

    plt.close(fig)

    return summary


def write_comparison_summary(summary_tables: dict, outdir: Path):
    rows = []

    for analysis_name, summary in summary_tables.items():
        if summary is None or summary.empty:
            continue

        row = summary.iloc[0].to_dict()
        row["analysis"] = analysis_name
        rows.append(row)

    if not rows:
        return None

    comp = pd.DataFrame(rows)
    cols = ["analysis"] + [c for c in comp.columns if c != "analysis"]
    comp = comp[cols]

    out = outdir / "canonical_WT_peptide_protein_gene_overlap_comparison_summary.tsv"
    comp.to_csv(out, sep="\t", index=False)

    print(f"[WRITE] {out}")

    return out


# ============================================================
# MAIN RUNNER
# ============================================================

def run_global_canonical_wt_overlap_figures(
    config: dict,
    outdir: Path,
) -> dict[str, Path]:
    """
    Génère les figures globales de recouvrement canonical WT.

    Sortie dans :
        outdir/canonical_wt_overlap_levels/
    """
    samples = get_config_samples(config)

    output_dir = Path(outdir) / "canonical_wt_overlap_levels"
    output_dir.mkdir(parents=True, exist_ok=True)

    print("\n[FIGURE] Recouvrement canonical WT peptide/protéine/gène")
    print(f"[SAMPLES] {samples}")
    print(f"[OUTPUT] {output_dir}")

    canonical_inputs = {
        sample: get_canonical_tsv(config, sample)
        for sample in samples
    }

    reference_fastas = {
        sample: get_reference_fasta(config, sample)
        for sample in samples
    }

    for sample, p in canonical_inputs.items():
        if not p.exists():
            raise FileNotFoundError(f"canonical.tsv introuvable pour {sample}: {p}")

    for sample, p in reference_fastas.items():
        if not p.exists():
            raise FileNotFoundError(f"FASTA introuvable pour {sample}: {p}")

    test_proteome_fasta = get_test_proteome_fasta(config)

    extra_fastas = []
    allowed_ensp = None
    has_test_proteome = False

    if test_proteome_fasta is not None:
        if test_proteome_fasta.exists():
            extra_fastas.append(test_proteome_fasta)
            allowed_ensp = read_ensp_whitelist_from_fasta(test_proteome_fasta)
            has_test_proteome = True
            print(f"[INFO] ENSP dans test proteome : {len(allowed_ensp)}")
        else:
            print(f"[WARN] test_proteome_fasta introuvable, figures filtrées ignorées : {test_proteome_fasta}")

    print("[LOAD] ENSP/ENST/ENSG -> gene mapping")

    id_to_gene_key, id_to_gene_symbol, id_to_ensg = build_identifier_to_gene_mapping(
        reference_fastas=reference_fastas,
        extra_fastas=extra_fastas,
    )

    all_long = []
    all_long_after_filter = []
    all_unmapped = []
    all_unmapped_after_filter = []

    for sample in samples:
        canonical_tsv = canonical_inputs[sample]

        print(f"\n[START] {sample}")
        print(f"[INPUT] {canonical_tsv}")

        long_df, unmapped_df = load_canonical_wt_detections(
            canonical_tsv=canonical_tsv,
            sample=sample,
            id_to_gene_key=id_to_gene_key,
            id_to_gene_symbol=id_to_gene_symbol,
            id_to_ensg=id_to_ensg,
            allowed_ensp=None,
        )

        print(f"[INFO] {sample} long rows : {len(long_df)}")

        all_long.append(long_df)

        if not unmapped_df.empty:
            all_unmapped.append(unmapped_df)

        if has_test_proteome:
            long_after_filter_df, unmapped_after_filter_df = load_canonical_wt_detections(
                canonical_tsv=canonical_tsv,
                sample=sample,
                id_to_gene_key=id_to_gene_key,
                id_to_gene_symbol=id_to_gene_symbol,
                id_to_ensg=id_to_ensg,
                allowed_ensp=allowed_ensp,
            )

            print(
                f"[INFO] {sample} long rows after test proteome filter : "
                f"{len(long_after_filter_df)}"
            )

            all_long_after_filter.append(long_after_filter_df)

            if not unmapped_after_filter_df.empty:
                all_unmapped_after_filter.append(unmapped_after_filter_df)

    all_df = pd.concat(all_long, ignore_index=True)

    all_long_out = output_dir / "ALL_canonical_WT_peptide_to_protein_to_gene_long.tsv"
    all_df.to_csv(all_long_out, sep="\t", index=False)
    print(f"[WRITE] {all_long_out}")

    if all_unmapped:
        unmapped_df = pd.concat(all_unmapped, ignore_index=True)
        unmapped_out = output_dir / "ALL_canonical_WT_unmapped_ENSP.tsv"
        unmapped_df.to_csv(unmapped_out, sep="\t", index=False)
        print(f"[WARN] Unmapped ENSP écrits : {unmapped_out}")

    count_tables = [
        summarize_counts(all_df, samples, "unfiltered")
    ]

    if has_test_proteome and all_long_after_filter:
        all_after_filter_df = pd.concat(all_long_after_filter, ignore_index=True)

        all_after_filter_out = (
            output_dir
            / "ALL_canonical_WT_peptide_to_protein_to_gene_long_after_test_proteome_filter.tsv"
        )

        all_after_filter_df.to_csv(all_after_filter_out, sep="\t", index=False)
        print(f"[WRITE] {all_after_filter_out}")

        if all_unmapped_after_filter:
            unmapped_after_filter_df = pd.concat(all_unmapped_after_filter, ignore_index=True)
            unmapped_after_filter_out = (
                output_dir
                / "ALL_canonical_WT_unmapped_ENSP_after_test_proteome_filter.tsv"
            )
            unmapped_after_filter_df.to_csv(unmapped_after_filter_out, sep="\t", index=False)
            print(f"[WARN] Unmapped ENSP après filtre écrits : {unmapped_after_filter_out}")

        count_tables.append(
            summarize_counts(all_after_filter_df, samples, "after_test_proteome_filter")
        )
    else:
        all_after_filter_df = None

    count_summary = pd.concat(count_tables, ignore_index=True)

    count_summary_out = output_dir / "canonical_WT_counts_by_level_per_sample.tsv"
    count_summary.to_csv(count_summary_out, sep="\t", index=False)
    print(f"[WRITE] {count_summary_out}")

    summary_tables = {}

    # ============================================================
    # UNFILTERED OVERLAPS
    # ============================================================

    summary_tables["peptide"] = make_venn(
        long_df=all_df,
        samples=samples,
        id_col="peptide",
        title="Peptides canoniques WT",
        outdir=output_dir,
        outprefix="canonical_WT_peptide_level_venn",
        level_label="peptide",
    )

    summary_tables["protein_accession_ENSP"] = make_venn(
        long_df=all_df,
        samples=samples,
        id_col="protein_id",
        title="Accessions protéiques canoniques WT",
        outdir=output_dir,
        outprefix="canonical_WT_protein_accession_level_venn",
        level_label="protein_accession_ENSP",
    )

    summary_tables["gene"] = make_venn(
        long_df=all_df,
        samples=samples,
        id_col="gene_key",
        title="Gènes canoniques WT",
        outdir=output_dir,
        outprefix="canonical_WT_gene_level_venn",
        level_label="gene",
    )

    # ============================================================
    # FILTERED OVERLAPS
    # ============================================================

    if all_after_filter_df is not None:
        summary_tables["peptide_after_test_proteome_filter"] = make_venn(
            long_df=all_after_filter_df,
            samples=samples,
            id_col="peptide",
            title="Peptides canoniques WT après filtre protéome",
            outdir=output_dir,
            outprefix="canonical_WT_peptide_level_after_test_proteome_filter_venn",
            level_label="peptide_after_test_proteome_filter",
        )

        summary_tables["protein_accession_ENSP_after_test_proteome_filter"] = make_venn(
            long_df=all_after_filter_df,
            samples=samples,
            id_col="protein_id",
            title="Accessions protéiques WT après filtre protéome",
            outdir=output_dir,
            outprefix="canonical_WT_protein_accession_level_after_test_proteome_filter_venn",
            level_label="protein_accession_ENSP_after_test_proteome_filter",
        )

        summary_tables["gene_after_test_proteome_filter"] = make_venn(
            long_df=all_after_filter_df,
            samples=samples,
            id_col="gene_key",
            title="Gènes canoniques WT après filtre protéome",
            outdir=output_dir,
            outprefix="canonical_WT_gene_level_after_test_proteome_filter_venn",
            level_label="gene_after_test_proteome_filter",
        )

    comparison_summary_out = write_comparison_summary(summary_tables, output_dir)

    outputs = {
        "canonical_wt_long": all_long_out,
        "canonical_wt_counts_summary": count_summary_out,
    }

    if comparison_summary_out is not None:
        outputs["canonical_wt_overlap_comparison_summary"] = comparison_summary_out

    print("\n[DONE] Figure globale canonical WT overlap terminée.")
    print(f"[OUTPUT] {output_dir}")

    return outputs


if __name__ == "__main__":
    raise SystemExit(
        "Ce module est conçu pour être appelé depuis run_all_figures.py.\n"
        "Utilise plutôt : python main.py -c config.yaml"
    )