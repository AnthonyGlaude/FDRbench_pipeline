#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Pipeline d'analyse protéomique modulaire.

Première partie : nettoyage
1. Cleanup basique
2. Vérification des I/L

Deuxième partie : filtrage des données
3. Filtrer les données

Troisième partie : génération des fichiers pour MS2PIP
4. Générer les fichiers input pour MS2PIP

Cette étape est un point d'arrêt :
les fichiers générés doivent ensuite être soumis au serveur MS2PIP.

Quatrième partie : analyses post-MS2PIP
5. Analyse des sorties MS2PIP
6. Analyse du support peptidique
7. QC final des peptides

Cinquième partie : figures
8. Génération des figures

Sixième partie optionnelle
9. Extraction des peaks d'intérêt
"""

from pathlib import Path
import argparse
import yaml

# MODE 1 — pré-MS2PIP
from modules.basic_cleanup import run_basic_cleanup
from modules.il_filtering import run_il_filtering
from modules.filter_data import run_filter_data
from modules.generate_ms2pip_input import run_generate_ms2pip_input

# MODE 2 — post-MS2PIP
from modules.run_ms2pip_postprocessing import run_ms2pip_experimental_scoring
from modules.peptide_support import run_peptide_support_analysis
from modules.final_qc import run_final_qc

# MODE 3 — figures
from modules.run_all_figures import run_figures


# ============================================================
# CONFIG
# ============================================================

def load_config(config_path):
    config_path = Path(config_path)

    if not config_path.exists():
        raise FileNotFoundError(f"Config introuvable : {config_path}")

    with open(config_path, "r") as f:
        return yaml.safe_load(f)


# ============================================================
# MS2PIP PLACEHOLDERS
# ============================================================

def prepare_ms2pip_intensity_folders(config):
    """
    Crée les dossiers où déposer manuellement les fichiers Spectronaut
    retournés par le serveur MS2PIP.

    Des fichiers vides sont créés comme placeholders.
    Ils doivent être remplacés ou remplis par les sorties MS2PIP.
    """

    path_cfg = config.get("path", {})

    if "ms2pip_output_root" in path_cfg:
        ms2pip_output_root = Path(path_cfg["ms2pip_output_root"])
    else:
        ms2pip_output_root = Path("../data/MS2PIP_intensities")

    ms2pip_output_root.mkdir(parents=True, exist_ok=True)

    expected_files = [
        "canonical_wt_ms2pip.spectronaut.tsv",
        "c_mut_ms2pip.spectronaut.tsv",
        "nc_wt_ms2pip.spectronaut.tsv",
        "nc_mut_ms2pip.spectronaut.tsv",
    ]

    for sample in config["samples"]:
        sample_dir = ms2pip_output_root / sample
        sample_dir.mkdir(parents=True, exist_ok=True)

        # Créer les fichiers vides attendus.
        # Le pipeline post-MS2PIP les ignore tant qu'ils sont vides.
        for name in expected_files:
            placeholder = sample_dir / name
            placeholder.touch(exist_ok=True)

        readme = sample_dir / "README_MS2PIP.txt"

        with open(readme, "w") as f:
            f.write(
                "Dossier de dépôt manuel des fichiers MS2PIP.\n\n"
                "Ce dossier contient des fichiers vides servant de placeholders.\n"
                "Après avoir lancé MS2PIP sur le serveur, remplace ou remplis ces fichiers "
                "avec les sorties Spectronaut correspondantes.\n\n"
                "Serveur MS2PIP :\n"
                "https://iomics.ugent.be/ms2pip/#run\n\n"
                "Paramètres à utiliser :\n"
                "- Fragmentation : HCD\n"
                "- Modèle : HCD (2021)\n"
                "- Format de sortie : Spectronaut CSV/TSV\n\n"
                "Fichiers attendus dans ce dossier :\n"
            )

            for name in expected_files:
                f.write(f"- {name}\n")

            f.write(
                "\nCorrespondance attendue :\n"
                "- canonical_wt_ms2pip.spectronaut.tsv : peptides canoniques WT\n"
                "- c_mut_ms2pip.spectronaut.tsv : peptides canoniques variants\n"
                "- nc_wt_ms2pip.spectronaut.tsv : peptides non canoniques WT\n"
                "- nc_mut_ms2pip.spectronaut.tsv : peptides non canoniques variants\n\n"
                "Une fois les fichiers remplis, relancer le pipeline et répondre 'n' à :\n"
                "Générer les fichiers input pour MS2PIP? (y/n)\n"
            )

    print(f"\nDossiers de dépôt MS2PIP créés dans : {ms2pip_output_root}")
    print("Des fichiers vides ont été créés comme placeholders.")
    print("Remplace ou remplis ces fichiers avec les sorties Spectronaut du serveur MS2PIP.")

    return ms2pip_output_root


# ============================================================
# MODES
# ============================================================

def run_mode_1_pre_ms2pip(config):
    """
    Mode 1 :
    Étapes 1 à 4.
    Génère les fichiers input pour MS2PIP, puis stop.
    """

    print("\n" + "#" * 100)
    print("MODE 1 — Préparation des fichiers MS2PIP")
    print("#" * 100)

    ##############################
    # Étape 1 - Cleanup basique
    ##############################
    print("\n[Étape 1] Cleanup basique")
    cleanup_outputs = run_basic_cleanup(config)

    ##############################
    # Étape 2 - Vérification I/L
    ##############################
    print("\n[Étape 2] Vérification des I/L")
    il_outputs = run_il_filtering(config)

    ##############################
    # Étape 3 - Filtrage / résumé des données
    ##############################
    print("\n[Étape 3] Filtrage des données")
    filter_outputs = run_filter_data(config)

    ##############################
    # Étape 4 - Fichiers input pour MS2PIP
    ##############################
    print("\n[Étape 4] Génération des fichiers input pour MS2PIP")
    peprec_outputs = run_generate_ms2pip_input(config)

    ##############################
    # Préparer les dossiers de dépôt MS2PIP
    ##############################
    ms2pip_output_root = prepare_ms2pip_intensity_folders(config)

    print("\nÉtapes 1 à 4 terminées.")
    print("Les fichiers input pour MS2PIP ont été générés.")
    print("Les dossiers de dépôt MS2PIP ont été préparés.")
    print(f"Dépôt attendu des sorties MS2PIP : {ms2pip_output_root}")
    print("Lance maintenant MS2PIP à l'extérieur du pipeline.")
    print("Ensuite, relance ce script et réponds 'n' à la première question.")

    return {
        "cleanup": cleanup_outputs,
        "il_filtering": il_outputs,
        "filter_data": filter_outputs,
        "ms2pip_input": peprec_outputs,
        "ms2pip_output_root": ms2pip_output_root,
    }


def run_mode_2_post_ms2pip(config):
    """
    Mode 2 :
    Étapes 5 à 7.
    Analyse les sorties MS2PIP, le support peptidique et le QC final.
    """

    print("\n" + "#" * 100)
    print("MODE 2 — Analyses post-MS2PIP")
    print("#" * 100)

    ##############################
    # Étape 5 - Analyse des sorties MS2PIP
    ##############################
    print("\n[Étape 5] Analyse des sorties MS2PIP")
    ms2pip_outputs = run_ms2pip_experimental_scoring(config)

    ##############################
    # Étape 6 - Analyse du support peptidique
    ##############################
    print("\n[Étape 6] Analyse du support peptidique")
    support_outputs = run_peptide_support_analysis(config)

    ##############################
    # Étape 7 - QC final des peptides
    ##############################
    print("\n[Étape 7] QC final des peptides")
    final_qc_outputs = run_final_qc(config)

    print("\nMode 2 terminé.")

    return {
        "ms2pip_scoring": ms2pip_outputs,
        "peptide_support": support_outputs,
        "final_qc": final_qc_outputs,
    }


def run_mode_3_figures(config):
    """
    Mode 3 :
    Étape 8.
    Génère les figures seulement.
    """

    print("\n" + "#" * 100)
    print("MODE 3 — Figures")
    print("#" * 100)

    ##############################
    # Étape 8 - Figures
    ##############################
    print("\n[Étape 8] Génération des figures")
    figure_outputs = run_figures(config)

    print("\nMode 3 terminé.")

    return {
        "figures": figure_outputs,
    }


# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Pipeline d'analyse protéomique modulaire"
    )

    parser.add_argument(
        "-c", "--config",
        default="../../config/config.yaml",
        help="Fichier config.yaml [défaut: ../../config/config.yaml]",
    )

    args = parser.parse_args()
    config = load_config(args.config)

    # ------------------------------------------------------------
    # Mode 1 : préparation MS2PIP
    # ------------------------------------------------------------
    run_before_ms2pip = input(
        "\nGénérer les fichiers input pour MS2PIP? (y/n) "
    ).strip().lower() == "y"

    if run_before_ms2pip:
        run_mode_1_pre_ms2pip(config)
        return

    # ------------------------------------------------------------
    # Mode 2 + 3 ou Mode 3 seulement
    # ------------------------------------------------------------
    run_post_ms2pip = input(
        "\nRouler les analyses post-MS2PIP avant les figures? (y/n) "
    ).strip().lower() == "y"

    if run_post_ms2pip:
        run_mode_2_post_ms2pip(config)

        run_figures_after = input(
            "\nGénérer les figures maintenant? (y/n) "
        ).strip().lower() == "y"

        if run_figures_after:
            run_mode_3_figures(config)
        else:
            print("\nFigures non générées.")
            print("Tu pourras relancer le script plus tard et choisir le mode figures seulement.")

        print("\nPipeline post-MS2PIP terminé.")
        return

    # ------------------------------------------------------------
    # Mode 3 seulement
    # ------------------------------------------------------------
    run_figures_only = input(
        "\nGénérer uniquement les figures? (y/n) "
    ).strip().lower() == "y"

    if run_figures_only:
        run_mode_3_figures(config)
        print("\nPipeline figures terminé.")
        return

    print("\nAucune étape lancée.")


if __name__ == "__main__":
    main()
