#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Figure globale — Séquences ajoutées aux bases personnalisées.

Cette figure est globale, car elle compare SH-SY5Y et SK-N-Be2.

Sorties :
    outdir/
        figure_4_5_sequences_ajoutees_bases_personnalisees.png
        figure_4_5_sequences_ajoutees_bases_personnalisees.svg
        figure_4_5_sequences_ajoutees_bases_personnalisees.pdf
"""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Patch
from matplotlib.lines import Line2D
import numpy as np


# ============================================================
# DONNÉES
# ============================================================

COLOR_CANONICAL_MISSENSE = "#1177CC"
COLOR_NC_MISSENSE = "#16A6A1"
COLOR_N_TERM = "#7B4BA0"
COLOR_C_TERM = "#D8A12A"
COLOR_NEW_ORF = "#58AD45"
COLOR_FUSION = "#E76F51"


fomonet_left = [
    ("uORF", 178, "#6B8DB8"),
    ("uoORF", 59, "#B7C9DD"),
    ("intORF", 83, "#7FB685"),
    ("doORF", 1002, "#F0B6BC"),
    ("dORF", 4224, "#E58C90"),
    ("Pseudo-ORF", 6680, "#8C6A91"),
    ("lncRNA-ORF", 45872, "#9FBEA8"),
    ("ncRNA-ORF", 21546, "#BDBDBD"),
]

fomonet_right = [
    ("CDS non annotées", 14643, "#E8D37A"),
    ("3' extension", 4, "#F28E54"),
    ("3' troncation", 109, "#D8A12A"),
    ("5' extension", 1509, "#F28E54"),
    ("5' troncation", 2947, "#D8A12A"),
    ("Autres", 316, "#E8C84A"),
]

cell_lines = ["SH-SY5Y", "SK-N-Be(2)"]
canonical_missense = [29672, 27523]
nc_missense = [4798, 3900]

rare_events = [
    ("Protéines de fusion", [25, 10], COLOR_FUSION),
    ("Nouveau ORF", [80, 57], COLOR_NEW_ORF),
    ("Altération C-terminale", [626, 480], COLOR_C_TERM),
    ("Altération N-terminale", [293, 242], COLOR_N_TERM),
]

venn_counts = {
    "SH_only": 23245,
    "shared": 12224,
    "SK_only": 19978,
}


# ============================================================
# HELPERS
# ============================================================

def fmt_int(n):
    return f"{int(n):,}".replace(",", " ")


def add_panel_label(ax, label, x=-0.08, y=1.08):
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=18,
        fontweight="bold",
        ha="left",
        va="top",
    )


def clean_spines(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def draw_track_bars(ax, data, max_val, x_offset=0, manual_labels=False):
    """
    Dessine des barres compactes sur fond gris.
    x_offset décale les barres vers la droite.
    manual_labels=True place les labels directement à gauche des barres.
    """
    labels = [f"{name} ({fmt_int(value)})" for name, value, _ in data]
    y = np.arange(len(data))[::-1]

    ax.barh(
        y,
        [max_val] * len(data),
        left=x_offset,
        color="#E6E6E6",
        edgecolor="white",
        height=0.74,
    )

    for i, (_name, value, color) in enumerate(data):
        yy = len(data) - 1 - i

        ax.barh(
            yy,
            value,
            left=x_offset,
            color=color,
            edgecolor="white",
            height=0.74,
        )

        if value / max_val > 0.18:
            x_text = x_offset + value * 0.55
            ha = "center"
        else:
            x_text = x_offset + value + max_val * 0.035
            ha = "left"

        ax.text(
            x_text,
            yy,
            fmt_int(value),
            va="center",
            ha=ha,
            fontsize=11,
            fontweight="bold",
            color="black",
        )

    if manual_labels:
        ax.set_yticks([])

        label_x = x_offset - max_val * 0.02

        for i, label in enumerate(labels):
            yy = len(data) - 1 - i
            ax.text(
                label_x,
                yy,
                label,
                va="center",
                ha="right",
                fontsize=12,
                color="black",
            )
    else:
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=12)
        ax.tick_params(axis="y", length=0)

    ax.set_xlim(0, x_offset + max_val * 1.05)
    ax.set_xticks([])

    for spine in ax.spines.values():
        spine.set_visible(False)


# ============================================================
# MAIN FIGURE
# ============================================================

def run_global_sequences_added_figure(
    config: dict,
    outdir: Path,
) -> dict[str, Path]:
    """
    Génère la figure globale des séquences ajoutées aux bases personnalisées.
    """

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    out_prefix = outdir / "figure_4_5_sequences_ajoutees_bases_personnalisees"

    print("\n[FIGURE] Séquences ajoutées aux bases personnalisées")
    print(f"[OUTPUT PREFIX] {out_prefix}")

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
    })

    fig = plt.figure(figsize=(13.5, 9), constrained_layout=False)

    gs = fig.add_gridspec(
        nrows=2,
        ncols=2,
        height_ratios=[1.30, 1.20],
        width_ratios=[1.35, 1.15],
        hspace=0.28,
        wspace=0.12,
    )

    # ============================================================
    # PANEL A — FOMOnet
    # ============================================================

    gsA = gs[0, :].subgridspec(1, 2, wspace=0.22)

    axA1 = fig.add_subplot(gsA[0, 0])
    axA2 = fig.add_subplot(gsA[0, 1])

    add_panel_label(axA1, "A", x=-0.12, y=1.18)

    axA1.text(
        0.00,
        1.16,
        "Répartition des ORF ajoutés par FOMOnet",
        transform=axA1.transAxes,
        fontsize=14,
        fontweight="bold",
        ha="left",
        va="center",
    )

    max_left = max(v for _, v, _ in fomonet_left)
    max_right = max(v for _, v, _ in fomonet_right)

    draw_track_bars(
        axA1,
        fomonet_left,
        max_left,
        x_offset=0,
        manual_labels=False,
    )

    draw_track_bars(
        axA2,
        fomonet_right,
        max_right,
        x_offset=max_right * 0.35,
        manual_labels=True,
    )

    # ============================================================
    # PANEL B — MUTATIONS
    # ============================================================

    gsB = gs[1, 0].subgridspec(
        1,
        2,
        width_ratios=[0.95, 1.15],
        wspace=0.65,
    )

    axB1 = fig.add_subplot(gsB[0, 0])
    axB2 = fig.add_subplot(gsB[0, 1])

    add_panel_label(axB1, "B", x=-0.23, y=1.18)

    axB1.text(
        -0.02,
        1.16,
        "Séquences protéiques associées aux mutations",
        transform=axB1.transAxes,
        fontsize=13,
        fontweight="bold",
        ha="left",
        va="center",
    )

    x = np.arange(len(cell_lines))
    bar_width = 0.62

    axB1.bar(
        x,
        canonical_missense,
        width=bar_width,
        color=COLOR_CANONICAL_MISSENSE,
        edgecolor="white",
        linewidth=0.8,
        label="ORF canonique missense",
    )

    axB1.bar(
        x,
        nc_missense,
        bottom=canonical_missense,
        width=bar_width,
        color=COLOR_NC_MISSENSE,
        edgecolor="white",
        linewidth=0.8,
        label="ORF non canonique missense",
    )

    for i in range(len(cell_lines)):
        axB1.text(
            x[i],
            canonical_missense[i] / 2,
            fmt_int(canonical_missense[i]),
            ha="center",
            va="center",
            color="white",
            fontsize=10,
        )

        axB1.text(
            x[i],
            canonical_missense[i] + nc_missense[i] / 2,
            fmt_int(nc_missense[i]),
            ha="center",
            va="center",
            color="white",
            fontsize=10,
            fontweight="bold",
        )

    axB1.set_title(
        "B1. Événements fréquents",
        loc="left",
        fontsize=11,
        fontweight="bold",
    )

    axB1.set_xticks(x)
    axB1.set_xticklabels(cell_lines, fontsize=12, fontweight="bold")
    axB1.tick_params(axis="x", labelsize=12)
    axB1.tick_params(axis="y", labelsize=11)
    axB1.set_ylabel("Nombre de séquences", fontsize=12)
    axB1.set_ylim(0, 40000)
    axB1.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.4)
    axB1.set_axisbelow(True)
    clean_spines(axB1)

    axB1.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.24),
        frameon=False,
        fontsize=12,
        ncol=1,
    )

    axB2.set_title(
        "B2. Événements rares",
        loc="left",
        fontsize=11,
        fontweight="bold",
    )

    ypos = np.arange(len(rare_events))[::-1]

    for idx, (event, values, color) in enumerate(rare_events):
        yy = len(rare_events) - 1 - idx

        sh_value, sk_value = values

        axB2.scatter(
            sh_value,
            yy + 0.08,
            marker="o",
            s=80,
            color=color,
            edgecolor="white",
            zorder=3,
        )

        axB2.scatter(
            sk_value,
            yy - 0.08,
            marker="s",
            s=80,
            color=color,
            edgecolor="white",
            zorder=3,
        )

        axB2.text(
            sh_value + 25,
            yy + 0.08,
            fmt_int(sh_value),
            va="center",
            ha="left",
            fontsize=9,
        )

        axB2.text(
            sk_value + 25,
            yy - 0.08,
            fmt_int(sk_value),
            va="center",
            ha="left",
            fontsize=9,
        )

    axB2.set_yticks(ypos)
    axB2.set_yticklabels([event for event, _, _ in rare_events], fontsize=11)
    axB2.tick_params(axis="y", pad=4, labelsize=11)
    axB2.tick_params(axis="x", labelsize=11)
    axB2.set_xlabel("Nombre de séquences", fontsize=12)
    axB2.set_xlim(0, 1000)
    axB2.set_ylim(-0.6, len(rare_events) - 0.4)
    axB2.grid(axis="x", linestyle="--", linewidth=0.6, alpha=0.4)
    axB2.set_axisbelow(True)
    clean_spines(axB2)

    shape_legend = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="black",
            linestyle="None",
            markersize=7,
            label="SH-SY5Y",
        ),
        Line2D(
            [0],
            [0],
            marker="s",
            color="black",
            linestyle="None",
            markersize=7,
            label="SK-N-Be(2)",
        ),
    ]

    event_legend = [
        Patch(facecolor=COLOR_N_TERM, label="Altération N-terminale"),
        Patch(facecolor=COLOR_C_TERM, label="Altération C-terminale"),
        Patch(facecolor=COLOR_NEW_ORF, label="Nouveau ORF"),
        Patch(facecolor=COLOR_FUSION, label="Protéines de fusion"),
    ]

    leg1 = axB2.legend(
        handles=shape_legend,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.22),
        frameon=False,
        fontsize=12,
        ncol=2,
    )

    axB2.add_artist(leg1)

    axB2.legend(
        handles=event_legend,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.42),
        frameon=False,
        fontsize=12,
        ncol=2,
    )

    # ============================================================
    # PANEL C — VENN MANUEL
    # ============================================================

    axC = fig.add_subplot(gs[1, 1])
    axC.set_aspect("equal", adjustable="box")
    axC.axis("off")

    axC.text(
        -0.08,
        1.10,
        "C",
        transform=axC.transAxes,
        fontsize=18,
        fontweight="bold",
        ha="left",
        va="top",
    )

    axC.text(
        0.00,
        1.05,
        "Recouvrement des séquences mutées entre lignées",
        transform=axC.transAxes,
        fontsize=13,
        fontweight="bold",
        ha="left",
        va="top",
    )

    left_center = (0.36, 0.58)
    right_center = (0.70, 0.58)
    radius = 0.40

    left_fill = "#F2B6B6"
    left_edge = "#D94F5C"
    right_fill = "#B7D9F2"
    right_edge = "#1177CC"
    shared_fill = "#D9D9D9"

    left_bg = Circle(
        left_center,
        radius,
        facecolor=left_fill,
        edgecolor="none",
        alpha=0.35,
        zorder=1,
    )

    right_bg = Circle(
        right_center,
        radius,
        facecolor=right_fill,
        edgecolor="none",
        alpha=0.35,
        zorder=1,
    )

    axC.add_patch(left_bg)
    axC.add_patch(right_bg)

    shared_patch = Circle(
        left_center,
        radius,
        facecolor=shared_fill,
        edgecolor="none",
        alpha=0.95,
        zorder=2,
    )

    shared_patch.set_clip_path(right_bg)
    axC.add_patch(shared_patch)

    left_outline = Circle(
        left_center,
        radius,
        facecolor="none",
        edgecolor=left_edge,
        linewidth=1.6,
        zorder=3,
    )

    right_outline = Circle(
        right_center,
        radius,
        facecolor="none",
        edgecolor=right_edge,
        linewidth=1.6,
        zorder=3,
    )

    axC.add_patch(left_outline)
    axC.add_patch(right_outline)

    axC.text(
        left_center[0] - 0.13,
        0.96,
        "SH-SY5Y",
        color=left_edge,
        fontsize=15,
        fontweight="bold",
        ha="center",
    )

    axC.text(
        right_center[0] + 0.13,
        0.96,
        "SK-N-Be(2)",
        color=right_edge,
        fontsize=15,
        fontweight="bold",
        ha="center",
    )

    axC.text(
        left_center[0] - 0.15,
        left_center[1],
        fmt_int(venn_counts["SH_only"]),
        ha="center",
        va="center",
        fontsize=15,
        fontweight="bold",
    )

    axC.text(
        0.53,
        left_center[1],
        fmt_int(venn_counts["shared"]),
        ha="center",
        va="center",
        fontsize=15,
        fontweight="bold",
    )

    axC.text(
        right_center[0] + 0.15,
        right_center[1],
        fmt_int(venn_counts["SK_only"]),
        ha="center",
        va="center",
        fontsize=15,
        fontweight="bold",
    )

    axC.text(
        0.53,
        0.08,
        "38 % des séquences protéiques uniques\nassociées aux mutations sont partagées",
        ha="center",
        va="center",
        fontsize=13,
        style="italic",
    )

    axC.set_xlim(-0.08, 1.14)
    axC.set_ylim(0.02, 1.08)

    # ============================================================
    # EXPORT
    # ============================================================

    outputs = {}

    for ext in ["png", "svg", "pdf"]:
        out = out_prefix.with_suffix(f".{ext}")
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(f"[WRITE] {out}")
        outputs[f"global_sequences_added_{ext}"] = out

    plt.close(fig)

    return outputs


if __name__ == "__main__":
    raise SystemExit(
        "Ce module est conçu pour être appelé depuis run_all_figures.py.\n"
        "Utilise plutôt : python main.py -c config.yaml"
    )