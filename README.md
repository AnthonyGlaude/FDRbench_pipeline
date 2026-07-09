# FDRBench pipeline

Ce dépôt contient un pipeline protéogénomique basé sur **Snakemake**, **MSFragger**, **PeptideProphet**, **FDRBench** et une étape de **post-analyse modulaire en Python**.

Le pipeline principal effectue l’analyse protéomique à partir de fichiers `.raw` jusqu’à la génération des fichiers finaux `canonical_final.tsv` et `nc_final.tsv`. Une seconde étape, située dans `workflow/post_analysis/`, permet ensuite de nettoyer, filtrer, annoter et résumer les identifications obtenues.

Ce pipeline a été adapté à partir du pipeline de Félix-Antoine Trifiro.

---

# 1. Structure générale du dépôt

La structure attendue du dépôt est la suivante :

```text
FDRbench_pipeline/
├── config/
│   ├── config.yaml
│   ├── closed_fragger.params
│   ├── closed_entrapment_fragger.params
│   ├── comet.params
│   └── comet_entrapment.params
│
├── envs/
│   ├── fdrbench.yml
│   └── mopepgen.yml
│
└── workflow/
    ├── Snakefile
    ├── data/
    │   ├── raw/
    │   ├── mzml/
    │   ├── fasta/
    │   └── contam/
    │
    ├── logs/
    ├── results/
    ├── scripts/
    └── post_analysis/
        ├── main.py
        ├── modules/
        └── results/

# 2. Prérequis

## 2.1 Snakemake

Le pipeline a été utilisé avec la version:

```text
7.32.4
```

Il est recommandé d’utiliser une version équivalente de Snakemake pour éviter des différences de comportement entre versions.

## 2.2 Singularity / Apptainer

Le pipeline utilise des conteneurs Singularity pour certaines étapes, notamment :

```text
docker://chambm/pwiz-skyline-i-agree-to-the-vendor-licenses
docker://spctools/tpp:version6.1
docker://spctools/tpp:version7.2.0
```



## 2.3 Conda

Le pipeline utilise aussi des environnements Conda définis dans :

```text
envs/
```

## 2.4 Java

MSFragger et FDRBench nécessitent Java.

Sur le serveur, Java est chargé dans certaines règles avec :

```bash
module load java
```

Il faut donc vérifier que le module Java est disponible dans l’environnement de calcul.

---

## 2.5 MSFragger

Le pipeline nécessite d’avoir **MSFragger** installé au préalable.

Le fichier `.jar` de MSFragger n’est pas inclus automatiquement dans le dépôt. Il faut le télécharger séparément, puis indiquer son chemin dans `config/config.yaml`.


## 2.6 FDRBench

Le pipeline nécessite aussi le fichier `.jar` de FDRBench.

Le chemin doit être indiqué dans `config/config.yaml`.

Exemple :

```yaml
fdrbench: "../../tools/FDRBench/FDRBench.jar"
```

---

## 2.7 Trans-Proteomic Pipeline

Le pipeline utilise TPP via des conteneurs Singularity pour :

* `decoyFASTA`;
* `xinteract`;
* PeptideProphet.

Les règles concernées utilisent notamment :

```text
docker://spctools/tpp:version6.1
docker://spctools/tpp:version7.2.0
```

---

## 2.8 ProteoWizard / msconvert

La conversion des fichiers `.raw` vers `.mzML` utilise `msconvert` via le conteneur :

```text
docker://chambm/pwiz-skyline-i-agree-to-the-vendor-licenses
```

---

# 3. Configuration

Le fichier principal de configuration est :

```text
config/config.yaml
```

Il contient notamment :

* les datasets;
* les samples;
* les chemins vers MSFragger;
* les chemins vers FDRBench;
* les fichiers de paramètres MSFragger;
* les seeds utilisées pour FDRBench.

---

# 4. Fichiers d’entrée

## 4.1 Fichiers RAW

Les fichiers `.raw` doivent être placés dans :

```text
workflow/data/raw/<dataset>/<sample>.raw
```

---

## 4.2 Bases FASTA

Les bases protéiques personnalisées doivent être placées dans :

```text
workflow/data/fasta/<dataset>/<dataset>.fasta
```


Ces bases peuvent contenir des protéines canoniques, non canoniques, variants, contaminants ajoutés plus tard, ou toute autre séquence utilisée pour l’analyse.

---

## 4.3 Contaminants

Le fichier de contaminants doit être placé dans :

```text
workflow/data/contam/crap_contam.fasta
```

Il sera ajouté automatiquement à la base protéique avec la règle :

```text
add_contams
```

---

# 5. Lancer le pipeline Snakemake

Le pipeline Snakemake doit être lancé depuis le dossier :

```text
workflow/
```

Commande utilisée :

```bash
cd workflow

snakemake --cores 12 --use-singularity --use-conda --singularity-args "-B /project/def-marie87/glaudea/FDRbench_pipeline"
```

L’argument :

```bash
--singularity-args "-B /project/def-marie87/glaudea/FDRbench_pipeline"
```

sert à monter le dossier du projet dans les conteneurs Singularity.

Si le dépôt est déplacé ailleurs, il faut modifier ce chemin.

---

# 6. Pipeline Snakemake principal

Le pipeline Snakemake effectue les étapes suivantes :

```text
RAW
 |
 v
mzML
 |
 v
FASTA + contaminants
 |
 v
base entrapment FDRBench
 |
 v
base decoy
 |
 v
MSFragger
 |
 v
PeptideProphet
 |
 v
conversion pepXML vers TSV
 |
 v
FDRBench FDP
 |
 v
canonical_final.tsv / nc_final.tsv
```

---

## 6.1 Règle `all`

La règle `all` définit les sorties finales attendues.

Sorties principales :

```text
workflow/results/<dataset>/ms_fragger/final/<sample>/canonical_final.tsv
workflow/results/<dataset>/ms_fragger/final/<sample>/nc_final.tsv
```

Ces fichiers sont les résultats finaux du pipeline Snakemake.

---

## 6.2 Ajout des contaminants

Règle :

```text
add_contams
```

Entrées :

```text
workflow/data/fasta/<dataset>/<dataset>.fasta
workflow/data/contam/crap_contam.fasta
```

Sortie :

```text
workflow/data/fasta/<dataset>/<dataset>_contam.fasta
```

Cette étape concatène la base protéique avec les contaminants.

Elle remplace aussi les caractères `*` par `X` dans la base finale.

---

## 6.3 Conversion RAW vers mzML

Règle :

```text
convert_raw
```

Entrée :

```text
workflow/data/raw/<dataset>/<sample>.raw
```

Sortie :

```text
workflow/data/mzml/<dataset>/<sample>.mzML
```

Cette étape utilise `msconvert`.

Options utilisées :

```text
--mzML
--zlib
peakPicking
zeroSamples removeExtra 1-
```

---

## 6.4 Génération de la base entrapment

Règle :

```text
generate_entrapment_db
```

Entrée :

```text
workflow/data/fasta/<dataset>/<dataset>_contam.fasta
```

Sorties :

```text
workflow/results/database/<dataset>/<dataset>_FDR_entrapment_<seed>.txt
workflow/results/database/<dataset>/<dataset>_FDR_entrapment_<seed>.fasta
```

Cette étape utilise FDRBench pour générer une base entrapment au niveau peptide.

Paramètres principaux :

```text
-enzyme 1
-fix_nc c
-level peptide
-uniprot
-miss_c 2
-minLength 9
-maxLength 35
-seed <seed>
```

---

## 6.5 Modification des headers entrapment

Règle :

```text
modify_headers
```

Script :

```text
workflow/scripts/modify_header.py
```

Sortie :

```text
workflow/results/database/<dataset>/<dataset>_<seed>_entrap_FDR_headers.fasta
```

Cette étape clarifie les headers de la base entrapment.

---

## 6.6 Génération des decoys pour la base entrapment

Règle :

```text
decoys_ent_db
```

Entrée :

```text
workflow/results/database/<dataset>/<dataset>_<seed>_entrap_FDR_headers.fasta
```

Sortie :

```text
workflow/results/database/<dataset>/<dataset>_<seed>_FDR_entrapment_decoys.fasta
```

Cette étape utilise `decoyFASTA` avec le préfixe :

```text
DECOY_
```

Le script suivant est ensuite appelé :

```text
workflow/scripts/reverse_decoy.py
```

---

## 6.7 Recherche MSFragger entrapment

Règle :

```text
msfragger_entrapment
```

Entrées :

```text
workflow/data/mzml/<dataset>/<sample>.mzML
workflow/results/database/<dataset>/<dataset>_<seed>_FDR_entrapment_decoys.fasta
```

Sortie :

```text
workflow/results/<dataset>/ms_fragger/msfragger_entrapment_search_<seed>/<sample>/<sample>_entrapment.pepXML
```

Cette étape lance MSFragger avec le fichier de paramètres défini dans :

```yaml
fragger_entrapment_params
```

Le fichier `.params` est copié dans le dossier de sortie, puis modifié automatiquement pour indiquer :

```text
database_name
num_threads
```

---

## 6.8 PeptideProphet entrapment

Règle :

```text
pep_prophet_entrapment
```

Entrée :

```text
workflow/results/<dataset>/ms_fragger/msfragger_entrapment_search_<seed>/<sample>/<sample>_entrapment.pepXML
```

Sortie :

```text
workflow/results/<dataset>/ms_fragger/msfragger_entrapment_search_<seed>/<sample>/<sample>_entrapment_prophet.pepXML
```

Cette étape utilise `xinteract` avec le préfixe decoy :

```text
-dDECOY_
```

---

## 6.9 Génération des decoys pour la recherche classique

Règle :

```text
decoys_db
```

Entrée :

```text
workflow/data/fasta/<dataset>/<dataset>_contam.fasta
```

Sortie :

```text
workflow/results/database/<dataset>/<dataset>_contam_decoys.fasta
```

Cette base est utilisée pour la recherche MSFragger classique.

---

## 6.10 Recherche MSFragger classique

Règle :

```text
msfragger_search
```

Entrées :

```text
workflow/data/mzml/<dataset>/<sample>.mzML
workflow/results/database/<dataset>/<dataset>_contam_decoys.fasta
```

Sortie :

```text
workflow/results/<dataset>/ms_fragger/search/<sample>/<sample>.pepXML
```

Cette étape lance MSFragger avec le fichier de paramètres défini dans :

```yaml
fragger_params
```

---

## 6.11 PeptideProphet classique

Règle :

```text
pep_prophet
```

Entrée :

```text
workflow/results/<dataset>/ms_fragger/search/<sample>/<sample>.pepXML
```

Sortie :

```text
workflow/results/<dataset>/ms_fragger/search/<sample>/<sample>_prophet.pepXML
```

Cette étape produit les résultats PeptideProphet utilisés pour la filtration finale.

---

## 6.12 Génération des fichiers pickle

Règle :

```text
pickles
```

Script :

```text
workflow/scripts/pickles.py
```

Sorties :

```text
workflow/results/database/<dataset>/<dataset>_<seed>_check_prot.pkl
workflow/results/database/<dataset>/<dataset>_<seed>_entrapment_pep.pkl
```

Ces fichiers servent à relier les peptides et protéines aux catégories canoniques et non canoniques.

---

## 6.13 Conversion des résultats entrapment en TSV

Règle :

```text
convert_xml_tsv
```

Script :

```text
workflow/scripts/tsvs_for_FDR.py
```

Entrée :

```text
workflow/results/<dataset>/ms_fragger/msfragger_entrapment_search_<seed>/<sample>/<sample>_entrapment_prophet.pepXML
```

Sorties :

```text
workflow/results/<dataset>/ms_fragger/msfragger_entrapment_search_<seed>/<sample>/canonicals.tsv
workflow/results/<dataset>/ms_fragger/msfragger_entrapment_search_<seed>/<sample>/non_canonicals.tsv
```

Cette étape sépare les identifications canoniques et non canoniques pour le calcul du FDP avec FDRBench.

---

## 6.14 Conversion de la recherche classique en TSV

Règle :

```text
convert_xml_tsv_search
```

Script :

```text
workflow/scripts/tsvs_for_FDR_search.py
```

Entrée :

```text
workflow/results/<dataset>/ms_fragger/search/<sample>/<sample>_prophet.pepXML
```

Sorties :

```text
workflow/results/<dataset>/ms_fragger/search/<sample>/canonicals_psm.tsv
workflow/results/<dataset>/ms_fragger/search/<sample>/non_canonicals_psm.tsv
```

Ces fichiers contiennent les PSMs de la recherche classique.

---

## 6.15 Calcul du FDP avec FDRBench

Règle :

```text
calculate_fdp
```

Entrées :

```text
canonicals.tsv
non_canonicals.tsv
<dataset>_FDR_entrapment_<seed>.txt
```

Sorties :

```text
workflow/results/<dataset>/ms_fragger/FDP/<sample>/can_FDR_fdrbench_<seed>.txt
workflow/results/<dataset>/ms_fragger/FDP/<sample>/nc_FDR_fdrbench_<seed>.txt
```

Cette étape utilise FDRBench pour estimer le FDP séparément pour les identifications canoniques et non canoniques.

---

## 6.16 Génération des fichiers finaux

Règle :

```text
final_output
```

Script :

```text
workflow/scripts/final_filtration.py
```

Sorties :

```text
workflow/results/<dataset>/ms_fragger/final/<sample>/canonical_final.tsv
workflow/results/<dataset>/ms_fragger/final/<sample>/nc_final.tsv
```

Ces deux fichiers sont les sorties principales de la portion Snakemake du pipeline.

---

# 7. Post-analyse modulaire

La post-analyse se trouve dans :

```text
workflow/post_analysis/
```

Elle est lancée avec :

```bash
cd workflow/post_analysis

python main.py -c ../../config/config.yaml
```

Cette étape prend comme entrée les fichiers finaux produits par Snakemake :

```text
workflow/results/<dataset>/ms_fragger/final/<sample>/canonical_final.tsv
workflow/results/<dataset>/ms_fragger/final/<sample>/nc_final.tsv
```

---

## 7.1 Objectif de la post-analyse

La post-analyse permet de :

* nettoyer les fichiers finaux;
* standardiser les colonnes;
* regrouper les résultats par dataset;
* distinguer les catégories canoniques, variants, ncORF WT et ncORF variants;
* filtrer les ambiguïtés I/L;
* préparer les entrées pour MS2PIP;
* intégrer les résultats MS2PIP;
* appliquer des règles de validation supplémentaires;
* exporter les spectres;
* générer les figures finales.

---

## 7.2 Point d’arrêt MS2PIP

Certaines étapes de la post-analyse dépendent de MS2PIP.

Le script peut donc demander au lancement :

```text
Lancer prédiction (y/n)?
```

Ce point d’arrêt est volontaire. Il permet de générer les fichiers d’entrée MS2PIP, d’exécuter MS2PIP séparément, puis de relancer la post-analyse pour intégrer les résultats.

Workflow recommandé :

```text
1. Lancer Snakemake
2. Obtenir canonical_final.tsv et nc_final.tsv
3. Lancer post_analysis/main.py
4. Générer les entrées MS2PIP
5. Exécuter MS2PIP
6. Relancer post_analysis/main.py
7. Générer les fichiers filtrés, les spectres et les figures
```

---

## 7.3 Sorties de post-analyse

Les sorties de post-analyse sont écrites dans des dossiers organisés par dataset.

Exemples :

```text
workflow/post_analysis/results/<dataset>/01_basic_cleanup/
workflow/post_analysis/results/<dataset>/02_IL_filtering/
workflow/post_analysis/results/<dataset>/04_ms2pip_input/
workflow/post_analysis/results/<dataset>/09_spectra/
workflow/post_analysis/results/<dataset>/figures/
```

Les noms exacts peuvent varier selon les modules activés dans `main.py`.

---

# 8. Sorties principales

## 8.1 Sorties Snakemake

```text
workflow/results/<dataset>/ms_fragger/final/<sample>/canonical_final.tsv
workflow/results/<dataset>/ms_fragger/final/<sample>/nc_final.tsv
```

## 8.2 Sorties post-analyse

```text
workflow/post_analysis/results/<dataset>/
```

Exemples possibles :

```text
01_basic_cleanup/
02_IL_filtering/
04_ms2pip_input/
09_spectra/
figures/
```

---

# 9. Logs

Les logs du pipeline Snakemake sont écrits dans :

```text
workflow/logs/
```

Exemples :

```text
workflow/logs/<dataset>/<sample>/convert_raw.log
workflow/logs/<dataset>/entrapment/
workflow/logs/<dataset>/database/
workflow/logs/<dataset>/classic/
workflow/logs/<dataset>/convert/
workflow/logs/<dataset>/<sample>/FINAL.log
```

Ces fichiers doivent être consultés en cas d’erreur.

---

# 10. Workflow résumé

```text
Fichiers RAW
    |
    v
convert_raw
    |
    v
mzML
    |
    v
add_contams
    |
    v
FASTA + contaminants
    |
    +-----------------------------+
    |                             |
    v                             v
generate_entrapment_db        decoys_db
    |                             |
    v                             v
modify_headers                MSFragger classique
    |                             |
    v                             v
decoys_ent_db                 PeptideProphet classique
    |                             |
    v                             v
MSFragger entrapment          canonicals_psm.tsv
    |                          non_canonicals_psm.tsv
    v
PeptideProphet entrapment
    |
    v
canonicals.tsv / non_canonicals.tsv
    |
    v
FDRBench FDP
    |
    v
final_output
    |
    v
canonical_final.tsv / nc_final.tsv
    |
    v
post_analysis/main.py
    |
    v
fichiers filtrés, spectres et figures
```

---

---

## Lancer le pipeline complet Snakemake

```bash
cd workflow

snakemake --cores 12 --use-singularity --use-conda --singularity-args "-B /project/def-marie87/glaudea/FDRbench_pipeline"
```

---

## Lancer une sortie spécifique

Exemple :

```bash
cd workflow

snakemake results/SH-SY5Y/ms_fragger/final/SHSY5Y_rep1/nc_final.tsv \
  --cores 12 \
  --use-singularity \
  --use-conda \
  --singularity-args "-B /project/def-marie87/glaudea/FDRbench_pipeline"
```

---

## Lancer la post-analyse

```bash
cd workflow/post_analysis

python main.py -c ../../config/config.yaml
```

---


ne correspondent pas aux datasets. Ce sont des catégories biologiques ou analytiques utilisées dans la post-analyse.


