"""
Étape 8 — Génération des figures.

Ce module orchestre les figures finales du pipeline.

Sorties :
    analysis_outdir/{sample}/08_figures/
        figures propres à chaque lignée

    analysis_outdir/08_figures_global/
        figures globales comparant plusieurs lignées
"""

from pathlib import Path

from modules.figures.figure_global_sequences_added import (
    run_global_sequences_added_figure,
)

from modules.figures.figure_global_ms_overlap import (
    run_global_ms_overlap_figures,
)

from modules.figures.figure_global_canonical_wt_overlap import (
    run_global_canonical_wt_overlap_figures,
)

def get_analysis_outdir(config: dict) -> Path:
    path_cfg = config.get("path", {})

    if "analysis_outdir" not in path_cfg:
        raise KeyError(
            "Chemin analysis_outdir absent du config.yaml. Ajoute :\n"
            "path:\n"
            "  analysis_outdir: results"
        )

    return Path(path_cfg["analysis_outdir"])


def get_sample_figdir(config: dict, sample: str) -> Path:
    """
    Dossier des figures propres à une lignée.

    Exemple :
        results/SH-SY5Y/08_figures/
        results/SK-N-Be2/08_figures/
    """
    analysis_outdir = get_analysis_outdir(config)
    return analysis_outdir / sample / "08_figures"


def get_global_figdir(config: dict) -> Path:
    """
    Dossier des figures globales.

    Exemple :
        results/08_figures_global/
    """
    analysis_outdir = get_analysis_outdir(config)
    return analysis_outdir / "08_figures_global"


def run_figures(config: dict) -> dict[str, Path]:
    """
    Lance toutes les figures finales.

    Pour l'instant :
        - figure globale des séquences ajoutées aux bases personnalisées

    Plus tard :
        - figures par lignée
        - figures QC final
        - figures support peptidique
        - figures NDP
    """

    print("\n" + "#" * 100)
    print("[Étape 8] Génération des figures")
    print("#" * 100)

    outputs: dict[str, Path] = {}

    # ============================================================
    # Figures globales
    # ============================================================

    global_figdir = get_global_figdir(config)
    global_figdir.mkdir(parents=True, exist_ok=True)

    print(f"\n[GLOBAL FIGURES] {global_figdir}")

    # figure 1 - proteines 
    global_outputs = run_global_sequences_added_figure(
        config=config,
        outdir=global_figdir,
    )

    outputs.update(global_outputs)
    
    #figure 2 - venn
    ms_overlap_outputs = run_global_ms_overlap_figures(
        config=config,
        outdir=global_figdir,
    )

    outputs.update(ms_overlap_outputs)
    
    
    canonical_wt_overlap_outputs = run_global_canonical_wt_overlap_figures(
    config=config,
    outdir=global_figdir,
)

    outputs.update(canonical_wt_overlap_outputs)
    
    # ============================================================
    # Figures par lignée — à ajouter progressivement
    # ============================================================

    for sample in config["samples"]:
        sample_figdir = get_sample_figdir(config, sample)
        sample_figdir.mkdir(parents=True, exist_ok=True)

        print(f"\n[SAMPLE FIGURES] {sample} -> {sample_figdir}")

        # Exemple futur :
        # qc_outputs = run_sample_qc_figures(
        #     config=config,
        #     sample=sample,
        #     outdir=sample_figdir,
        # )
        # outputs.update(qc_outputs)

    print("\n[DONE] Étape 8 — Figures terminées.")

    return outputs
