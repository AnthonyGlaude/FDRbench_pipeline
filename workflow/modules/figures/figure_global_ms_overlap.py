#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Figures globales — Recouvrement des identifications MS.

Ce module génère :

1. Figure 4 panneaux :
    - Gènes canoniques détectés
    - Gènes associés à des variants canoniques
    - Gènes associés à des ncORF non mutés
    - Gènes associés à des ncORF variants

2. Figure Venn :
    - Séquences protéiques mutées

Les noms de lignées viennent directement de config["samples"].

Sorties :
    outdir/
        figure_4_7_overlap_MS_identifications_4panels.png
        figure_4_7_overlap_MS_identifications_4panels.svg
        figure_4_7_overlap_MS_identifications_4panels.pdf

        venn_mutated_sequences_overlap.png
        venn_mutated_sequences_overlap.svg
        venn_mutated_sequences_overlap.pdf
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2


# ============================================================
# COLORS
# ============================================================

COLOR_LEFT = "#F2B6B6"
COLOR_RIGHT = "#B7C9DD"
COLOR_SHARED = "#D9D9D9"
EDGE_COLOR = "#4D4D4D"


# ============================================================
# GLOBAL DATA
# ============================================================
# Important :
# Les clés doivent correspondre exactement aux noms dans config["samples"].
#
# Format des tuples :
#     (left_only, shared, right_only)
#
# Le "left" et "right" sont déterminés par l'ordre dans config["samples"].
# Si l'ordre est inversé, le script inverse automatiquement les valeurs.

MS_OVERLAP_COUNTS = {
    "canonical_genes": {
        "panel": "A",
        "title": "Gènes canoniques détectés",
        "counts": (625, 3244, 2373),
    },
    "canonical_variant_genes": {
        "panel": "B",
        "title": "Gènes associés à des variants canoniques",
        "counts": (22, 12, 46),
    },
    "nc_wt_genes": {
        "panel": "C",
        "title": "Gènes associés à des ncORF non mutés",
        "counts": (5, 10, 7),
    },
    "nc_variant_genes": {
        "panel": "D",
        "title": "Gènes associés à des ncORF variants",
        "counts": (2, 1, 1),
    },
}


MUTATED_SEQUENCE_OVERLAP = {
    "title": "Séquences protéiques mutées",
    "counts": (89359, 43624, 62049),
}

MUTATED_SEQUENCE_OVERLAP = {
    "title": "Séquences protéiques mutées",
    "counts": (89359, 43624, 62049),
}


# ============================================================
# HELPERS
# ============================================================

def fmt_int_space(n):
    return f"{int(n):,}".replace(",", " ")


def fmt_decimal_comma(x, ndigits=2):
    return f"{x:.{ndigits}f}".replace(".", ",")


def compute_oc(left_only, shared, right_only):
    """
    Overlap coefficient :
        OC = intersection / min(total_left, total_right)
    """
    left_total = left_only + shared
    right_total = right_only + shared
    denom = min(left_total, right_total)

    if denom == 0:
        return np.nan

    return shared / denom


def get_config_samples(config: dict) -> list[str]:
    """
    Retourne les samples dans l'ordre du config.yaml.
    """
    if "samples" not in config:
        raise KeyError("Section 'samples' absente du config.yaml.")

    samples = list(config["samples"].keys())

    if len(samples) != 2:
        raise ValueError(
            "Ces figures globales sont conçues pour exactement deux samples.\n"
            f"Samples reçus : {samples}"
        )

    return samples


def get_pair_data(pair_dict: dict, samples: list[str]) -> dict:
    """
    Récupère les données associées à la paire de samples.

    Si l'ordre dans config est inversé, les valeurs left_only/right_only
    sont inversées automatiquement.
    """
    pair = tuple(samples)
    reverse_pair = tuple(reversed(samples))

    if pair in pair_dict:
        return pair_dict[pair]

    if reverse_pair not in pair_dict:
        raise KeyError(
            "Aucune donnée trouvée pour la paire de samples : "
            f"{pair}. Ajoute cette paire dans le dictionnaire de données."
        )

    original_data = pair_dict[reverse_pair]
    reversed_data = {}

    for key, value in original_data.items():
        if "counts" not in value:
            reversed_data[key] = value
            continue

        left_only, shared, right_only = value["counts"]

        reversed_entry = dict(value)
        reversed_entry["counts"] = (right_only, shared, left_only)
        reversed_data[key] = reversed_entry

    return reversed_data


def style_venn(v):
    patch_specs = {
        "10": (COLOR_LEFT, 0.75),
        "01": (COLOR_RIGHT, 0.75),
        "11": (COLOR_SHARED, 0.85),
    }

    for region_id, (color, alpha) in patch_specs.items():
        patch = v.get_patch_by_id(region_id)

        if patch is not None:
            patch.set_color(color)
            patch.set_alpha(alpha)
            patch.set_edgecolor(EDGE_COLOR)
            patch.set_linewidth(1.8)


def draw_overlap_panel(
    ax,
    *,
    panel_label: str,
    title: str,
    left_label: str,
    right_label: str,
    left_only: int,
    shared: int,
    right_only: int,
):
    """
    Dessine un panneau Venn avec OC.
    """
    v = venn2(
        subsets=(left_only, right_only, shared),
        set_labels=(left_label, right_label),
        ax=ax,
    )

    style_venn(v)

    region_counts = {
        "10": left_only,
        "11": shared,
        "01": right_only,
    }

    total_union = left_only + shared + right_only

    if total_union <= 10:
        count_fontsize = 22
    elif total_union <= 100:
        count_fontsize = 18
    else:
        count_fontsize = 14

    for region_id, n in region_counts.items():
        label = v.get_label_by_id(region_id)

        if label is not None:
            label.set_text(fmt_int_space(n))
            label.set_fontsize(count_fontsize)
            label.set_fontweight("bold")
            label.set_color("black")

    for label in v.set_labels:
        if label is not None:
            label.set_fontsize(12)
            label.set_fontweight("bold")

    oc = compute_oc(left_only, shared, right_only)

    if np.isnan(oc):
        oc_text = "OC = NA"
    else:
        oc_text = f"OC = {fmt_decimal_comma(oc, 2)}"

    ax.set_title(
        f"{panel_label}. {title}",
        fontsize=14,
        fontweight="bold",
        pad=12,
    )

    ax.text(
        0.5,
        -0.10,
        oc_text,
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=13,
        fontweight="bold",
    )


def draw_single_venn(
    ax,
    *,
    title: str,
    left_label: str,
    right_label: str,
    left_only: int,
    shared: int,
    right_only: int,
):
    """
    Dessine une figure Venn simple.
    """
    v = venn2(
        subsets=(left_only, right_only, shared),
        set_labels=(left_label, right_label),
        ax=ax,
    )

    style_venn(v)

    for region_id, value in {
        "10": left_only,
        "11": shared,
        "01": right_only,
    }.items():
        label = v.get_label_by_id(region_id)

        if label is not None:
            label.set_text(fmt_int_space(value))
            label.set_fontsize(13)
            label.set_fontweight("bold")
            label.set_color("black")

    for label in v.set_labels:
        if label is not None:
            label.set_fontsize(13)
            label.set_fontweight("bold")

    ax.set_title(
        title,
        fontsize=15,
        fontweight="bold",
        pad=12,
    )


# ============================================================
# FIGURE 1 — 4 PANELS
# ============================================================

def run_ms_identification_overlap_4panel_figure(
    config: dict,
    outdir: Path,
) -> dict[str, Path]:
    """
    Génère la figure 4 panneaux des recouvrements MS.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    samples = get_config_samples(config)
    left_label = samples[0]
    right_label = samples[1]

    panel_data = MS_OVERLAP_COUNTS
    out_prefix = outdir / "figure_4_7_overlap_MS_identifications_4panels"

    print("\n[FIGURE] Recouvrement des identifications MS — 4 panneaux")
    print(f"[SAMPLES] {samples}")
    print(f"[OUTPUT PREFIX] {out_prefix}")

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    axes = axes.flatten()

    panel_order = [
        "canonical_genes",
        "canonical_variant_genes",
        "nc_wt_genes",
        "nc_variant_genes",
    ]

    for ax, key in zip(axes, panel_order):
        panel = panel_data[key]
        left_only, shared, right_only = panel["counts"]

        draw_overlap_panel(
            ax=ax,
            panel_label=panel["panel"],
            title=panel["title"],
            left_label=left_label,
            right_label=right_label,
            left_only=left_only,
            shared=shared,
            right_only=right_only,
        )

    plt.tight_layout()

    outputs = {}

    for ext in ["png", "svg", "pdf"]:
        out = out_prefix.with_suffix(f".{ext}")
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(f"[WRITE] {out}")
        outputs[f"ms_overlap_4panel_{ext}"] = out

    plt.close(fig)

    return outputs


# ============================================================
# FIGURE 2 — MUTATED SEQUENCES OVERLAP
# ============================================================

def run_mutated_sequences_overlap_figure(
    config: dict,
    outdir: Path,
) -> dict[str, Path]:
    """
    Génère la figure Venn des séquences protéiques mutées.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    samples = get_config_samples(config)
    left_label = samples[0]
    right_label = samples[1]

    mutated_data = MUTATED_SEQUENCE_OVERLAP

    title = mutated_data["title"]
    left_only, shared, right_only = mutated_data["counts"]

    out_prefix = outdir / "venn_mutated_sequences_overlap"

    print("\n[FIGURE] Recouvrement des séquences protéiques mutées")
    print(f"[SAMPLES] {samples}")
    print(f"[OUTPUT PREFIX] {out_prefix}")

    fig, ax = plt.subplots(figsize=(7, 6))

    draw_single_venn(
        ax=ax,
        title=title,
        left_label=left_label,
        right_label=right_label,
        left_only=left_only,
        shared=shared,
        right_only=right_only,
    )

    plt.tight_layout()

    outputs = {}

    for ext in ["png", "svg", "pdf"]:
        out = out_prefix.with_suffix(f".{ext}")
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(f"[WRITE] {out}")
        outputs[f"mutated_sequences_overlap_{ext}"] = out

    plt.close(fig)

    return outputs


# ============================================================
# MAIN RUNNER
# ============================================================

def run_global_ms_overlap_figures(
    config: dict,
    outdir: Path,
) -> dict[str, Path]:
    """
    Génère toutes les figures globales de recouvrement MS.
    """
    outputs = {}

    outputs.update(
        run_ms_identification_overlap_4panel_figure(
            config=config,
            outdir=outdir,
        )
    )

    outputs.update(
        run_mutated_sequences_overlap_figure(
            config=config,
            outdir=outdir,
        )
    )

    return outputs


if __name__ == "__main__":
    raise SystemExit(
        "Ce module est conçu pour être appelé depuis run_all_figures.py.\n"
        "Utilise plutôt : python main.py -c config.yaml"
    )