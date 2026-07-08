#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Étape 6 — Analyse du support peptidique.

Cette étape est lancée après le scoring expérimental MS2PIP.

Entrées attendues :

    analysis_outdir/{sample}/02_IL_filtering/
        canonical_variant.tsv
        nc_variant.tsv
        nc_wt.tsv

    ../data/fasta/{sample}/{sample}_contam.fasta

    results_root/{sample}/ms_fragger/search/*/
        canonicals_psm.tsv
        non_canonicals_psm.tsv

    path.seq_to_ids_pkl dans config.yaml :
        seq_to_ids_pkl:
          SH-SY5Y: /chemin/vers/SRR43_seq_to_ids_cut.pkl
          SK-N-Be2: /chemin/vers/SRR46_seq_to_ids_cut.pkl

Sorties :

    analysis_outdir/{sample}/06_peptide_support/
        {sample}_{category}_support.tsv
        ALL_categories_support.tsv
        summary_counts.tsv

Conçu pour être appelé depuis main.py :

    from modules.peptide_support import run_peptide_support_analysis
    support_outputs = run_peptide_support_analysis(config)
"""

from pathlib import Path
from collections import Counter, defaultdict
import ast
import pickle
import re

import pandas as pd
import numpy as np


# ============================================================
# CATEGORIES
# ============================================================

CATEGORIES = {
    "canonical_variant": {
        "input_name": "canonical_variant.tsv",
        "psm_name": "canonicals_psm.tsv",
        "reject_reason": "singleton_variant_without_protein_support",
    },
    "nc_variant": {
        "input_name": "nc_variant.tsv",
        "psm_name": "non_canonicals_psm.tsv",
        "reject_reason": "singleton_nc_variant_without_orf_support",
    },
    "nc_wt": {
        "input_name": "nc_wt.tsv",
        "psm_name": "non_canonicals_psm.tsv",
        "reject_reason": "singleton_nc_wt_without_orf_support",
    },
}


# ============================================================
# REGEX
# ============================================================

AA_ONLY = re.compile(r"[A-Z]")
ENSP_RE = re.compile(r"(ENSP\d+)")
ENST_RE = re.compile(r"(ENST\d+)")
ENS_ANY_RE = re.compile(r"(ENS[TPG]\d+)")
DUP_RE = re.compile(r"_dup(\d+)\b")


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


def get_results_root(config: dict) -> Path:
    """
    Dossier contenant les résultats MSFragger/FDRBench originaux.

    Exemple config :
        path:
          results_root: ../results

    Si absent, fallback :
        ../results
    """
    path_cfg = config.get("path", {})
    return Path(path_cfg.get("results_root", "../results"))


def get_custom_fasta(config: dict, sample: str) -> Path:
    """
    FASTA custom attendue pour chaque sample/cell line :

        ../data/fasta/{sample}/{sample}_contam.fasta

    Exemple :
        ../data/fasta/SH-SY5Y/SH-SY5Y_contam.fasta
        ../data/fasta/SK-N-Be2/SK-N-Be2_contam.fasta
    """
    path_cfg = config.get("path", {})

    fasta_root = Path(path_cfg.get("fasta_root", "../data/fasta"))

    return fasta_root / sample / f"{sample}_contam.fasta"


def get_seq_to_ids_pkl(config: dict, sample: str) -> Path:
    """
    Chemin du pickle seq_to_ids pour chaque sample.

    Exemple config :
        path:
          seq_to_ids_pkl:
            SH-SY5Y: /chemin/vers/SRR43_seq_to_ids_cut.pkl
            SK-N-Be2: /chemin/vers/SRR46_seq_to_ids_cut.pkl
    """
    path_cfg = config.get("path", {})
    pkl_cfg = path_cfg.get("seq_to_ids_pkl", {})

    if not isinstance(pkl_cfg, dict):
        raise TypeError(
            "path.seq_to_ids_pkl doit être un dictionnaire dans config.yaml."
        )

    if sample not in pkl_cfg:
        raise KeyError(
            f"seq_to_ids_pkl absent pour {sample} dans config.yaml.\n"
            "Ajoute par exemple :\n"
            "path:\n"
            "  seq_to_ids_pkl:\n"
            "    SH-SY5Y: /chemin/vers/SRR43_seq_to_ids_cut.pkl\n"
            "    SK-N-Be2: /chemin/vers/SRR46_seq_to_ids_cut.pkl"
        )

    return Path(pkl_cfg[sample])


def get_input_tsv(config: dict, sample: str, category: str) -> Path:
    """
    Entrée issue de l'étape I/L filtering :

        analysis_outdir/{sample}/02_IL_filtering/{input_name}

    Exemple :
        results/SH-SY5Y/02_IL_filtering/canonical_variant.tsv
    """
    analysis_outdir = get_analysis_outdir(config)

    return (
        analysis_outdir
        / sample
        / "02_IL_filtering"
        / CATEGORIES[category]["input_name"]
    )


def get_psm_glob(config: dict, sample: str, category: str) -> str:
    """
    PSM MSFragger/FDRBench attendus :

        results_root/{sample}/ms_fragger/search/*/{psm_name}

    Exemple :
        ../results/SH-SY5Y/ms_fragger/search/*/canonicals_psm.tsv
        ../results/SH-SY5Y/ms_fragger/search/*/non_canonicals_psm.tsv
    """
    results_root = get_results_root(config)
    psm_name = CATEGORIES[category]["psm_name"]

    return str(
        results_root
        / sample
        / "ms_fragger"
        / "search"
        / "*"
        / psm_name
    )


def get_support_outdir(config: dict, sample: str) -> Path:
    """
    Sortie par sample :

        analysis_outdir/{sample}/06_peptide_support/
    """
    analysis_outdir = get_analysis_outdir(config)

    return analysis_outdir / sample / "06_peptide_support"


# ============================================================
# GENERAL HELPERS
# ============================================================

def to_list(x):
    if isinstance(x, (list, tuple, set)):
        return [str(i).strip() for i in x if str(i).strip()]

    if isinstance(x, np.ndarray):
        return [str(i).strip() for i in x.tolist() if str(i).strip()]

    if x is None:
        return []

    try:
        if pd.isna(x):
            return []
    except Exception:
        pass

    s = str(x).strip()

    if not s or s.lower() in {"nan", "none"}:
        return []

    if s.startswith("[") and s.endswith("]"):
        try:
            v = ast.literal_eval(s)
            if isinstance(v, (list, tuple, set)):
                return [str(i).strip() for i in v if str(i).strip()]
        except Exception:
            pass

    return [p.strip() for p in re.split(r"[;,]\s*", s) if p.strip()]


def pep_key(x) -> str:
    """
    Normalisation peptide :
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

    s = str(x).upper()
    seq = "".join(AA_ONLY.findall(s))

    return seq.replace("I", "J").replace("L", "J")


def find_peptide_col(df: pd.DataFrame) -> str:
    candidates = [
        "peptide",
        "Peptide",
        "modified_peptide",
        "Modified Peptide",
        "ModifiedPeptide",
        "Peptide Sequence",
        "Peptide_Sequence",
        "ASSIGNED",
    ]

    for c in candidates:
        if c in df.columns:
            return c

    raise ValueError(
        f"Colonne peptide introuvable. Colonnes disponibles: {list(df.columns)}"
    )


def find_protein_col(df: pd.DataFrame) -> str:
    candidates = [
        "protein",
        "Protein",
        "Proteins",
        "protein_raw",
        "Protein Accession",
        "Protein_Accession",
    ]

    for c in candidates:
        if c in df.columns:
            return c

    raise ValueError(
        f"Colonne protein introuvable. Colonnes disponibles: {list(df.columns)}"
    )


# ============================================================
# PROTEIN / ORF ID NORMALIZATION
# ============================================================

def strip_version_from_ens_id(x: str) -> str:
    """
    ENSP00000379535.3 -> ENSP00000379535
    ENST00000422379.7 -> ENST00000422379
    """
    return re.sub(r"^(ENS[TPG]\d+)\.\d+$", r"\1", str(x))


def protein_aliases(x):
    """
    Crée des aliases pour matcher les IDs entre :
    - TSV d'entrée
    - PSM
    - FASTA custom

    Exemples :
        ENSP00000379535.3 -> ENSP00000379535
        ENSP00000415081_mut_C(G)_pos=... -> ENSP00000415081
        ENST00000422379_mut_C(G)_pos=... -> ENST00000422379
    """
    if x is None:
        return set()

    s = str(x).strip()

    if not s or s.lower() in {"nan", "none"}:
        return set()

    out = set()

    out.add(s)

    first = s.split()[0]
    out.add(first)
    out.add(strip_version_from_ens_id(first))

    pipe_first = first.split("|")[0]
    out.add(pipe_first)
    out.add(strip_version_from_ens_id(pipe_first))

    for m in ENS_ANY_RE.finditer(s):
        out.add(m.group(1))

    for m in re.finditer(r"(ENS[TPG]\d+\.\d+)", s):
        out.add(m.group(1))
        out.add(strip_version_from_ens_id(m.group(1)))

    return {v for v in out if v and v.lower() not in {"nan", "none"}}


def normalize_alias_set(items):
    out = set()

    for x in to_list(items):
        out.update(protein_aliases(x))

    return out


# ============================================================
# CLEAN PROTEIN LIST
# ============================================================

def clean_protein_list(protein_list):
    """
    Nettoie les accessions :
    - garde les entrées custom avec ORF
    - pour les ENSP, préfère la version sans _dupN
    - déduplique
    """
    protein_list = to_list(protein_list)

    if not protein_list:
        return []

    best_ensp = {}
    others = []

    for tok in protein_list:
        t = str(tok).strip()

        if not t:
            continue

        m = ENSP_RE.search(t)

        if not m:
            others.append(t)
            continue

        ensp = m.group(1)

        mdup = DUP_RE.search(t)
        rank = int(mdup.group(1)) if mdup else 0

        if ensp not in best_ensp:
            best_ensp[ensp] = (rank, t)
        else:
            prev_rank, _ = best_ensp[ensp]
            if rank < prev_rank:
                best_ensp[ensp] = (rank, t)

    cleaned = others + [
        tok for _, tok in sorted(best_ensp.values(), key=lambda x: x[0])
    ]

    seen = set()
    out = []

    for x in cleaned:
        if x not in seen:
            seen.add(x)
            out.append(x)

    return out


# ============================================================
# FASTA / PKL
# ============================================================

def load_fasta_id2seq(fasta_path: Path):
    """
    Charge une FASTA en dict :
        ID/alias -> séquence
    """
    id2seq = {}
    cur_header = None
    cur_seq = []

    def flush_record(header, seq_chunks):
        if header is None:
            return

        seq = "".join(seq_chunks)
        first = header.split()[0]

        keys = set()
        keys.add(first)
        keys.update(protein_aliases(first))
        keys.update(protein_aliases(header))

        for k in keys:
            if k and k not in id2seq:
                id2seq[k] = seq

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                flush_record(cur_header, cur_seq)
                cur_header = line[1:].strip()
                cur_seq = []
            else:
                cur_seq.append(line)

        flush_record(cur_header, cur_seq)

    return id2seq


def load_seq_to_ids_pkl(pkl_path: Path):
    with open(pkl_path, "rb") as f:
        return pickle.load(f)


# ============================================================
# TA / GROUP CLASSIFICATION
# ============================================================

def add_TA_column(df: pd.DataFrame, id2seq: dict, seq_to_ids: dict) -> pd.DataFrame:
    """
    Ajoute TA :
        protein accession -> sequence FASTA -> IDs partageant cette même séquence.
    """
    df = df.copy()

    def row_to_ta(row):
        proteins = to_list(row.get("protein", []))
        out = set()

        for pid in proteins:
            seq = None

            for alias in protein_aliases(pid):
                seq = id2seq.get(alias)
                if seq:
                    break

            if not seq:
                continue

            ids = seq_to_ids.get(seq) or seq_to_ids.get(seq.rstrip("*"))

            if not ids:
                continue

            protein_alias_set = normalize_alias_set(proteins)

            for x in ids:
                sx = str(x).strip()

                if not sx:
                    continue

                if sx in proteins:
                    continue

                if protein_aliases(sx) & protein_alias_set:
                    continue

                out.add(sx)

        return sorted(out)

    df["TA"] = df.apply(row_to_ta, axis=1)

    return df


def classify_unique_group(df: pd.DataFrame) -> pd.DataFrame:
    """
    Classe selon protein + TA :
    - unique : 1 ENSP et aucun ENST
    - group : >=2 ENSP ou >=1 ENST
    - unclassified : rien d'exploitable
    """
    df = df.copy()

    def classify_row(row):
        all_ids = to_list(row.get("protein", [])) + to_list(row.get("TA", []))

        ensps = set()
        ensts = set()

        for x in all_ids:
            for alias in protein_aliases(x):
                m_ensp = ENSP_RE.search(alias)
                m_enst = ENST_RE.search(alias)

                if m_ensp:
                    ensps.add(m_ensp.group(1))

                if m_enst:
                    ensts.add(m_enst.group(1))

        if len(ensps) == 1 and len(ensts) == 0:
            return "unique"

        if len(ensps) >= 2 or len(ensts) >= 1:
            return "group"

        return "unclassified"

    df["cat"] = df.apply(classify_row, axis=1)

    return df


def group_set_full_ids(row):
    proteins = to_list(row.get("protein", []))
    ta = to_list(row.get("TA", []))

    return normalize_alias_set(proteins + ta)


# ============================================================
# PSM INDEX
# ============================================================

def load_psm_files(psm_glob: str):
    if psm_glob.startswith("/"):
        psm_files = sorted(Path("/").glob(psm_glob[1:]))
    else:
        psm_files = sorted(Path().glob(psm_glob))

    if not psm_files:
        raise FileNotFoundError(f"Aucun fichier PSM trouvé avec glob: {psm_glob}")

    return psm_files


def read_psm_file(psm_path: Path):
    df = pd.read_csv(psm_path, sep="\t")

    pep_col = find_peptide_col(df)
    prot_col = find_protein_col(df)

    run = psm_path.parent.name

    out = pd.DataFrame({
        "run": run,
        "peptide_raw": df[pep_col],
        "pep_key": df[pep_col].map(pep_key),
        "protein_raw": df[prot_col],
    })

    out = out[out["pep_key"] != ""].copy()

    return out


def build_psm_index(psm_files):
    pep2psm = Counter()
    pep2prot = defaultdict(set)
    rows = []

    for p in psm_files:
        d = read_psm_file(p)

        for _, r in d.iterrows():
            k = r["pep_key"]

            if not k:
                continue

            pep2psm[k] += 1

            proteins = to_list(r["protein_raw"])

            for pid in proteins:
                aliases = protein_aliases(pid)

                for alias in aliases:
                    pep2prot[k].add(alias)

                base_aliases = sorted({
                    a for a in aliases
                    if ENSP_RE.fullmatch(a) or ENST_RE.fullmatch(a)
                })

                if base_aliases:
                    for alias in base_aliases:
                        rows.append({
                            "run": r["run"],
                            "pep_key": k,
                            "protein": alias,
                            "protein_raw": pid,
                        })
                else:
                    rows.append({
                        "run": r["run"],
                        "pep_key": k,
                        "protein": str(pid),
                        "protein_raw": pid,
                    })

    psm_long = pd.DataFrame(rows)

    return {
        "pep2psm": pep2psm,
        "pep2prot": pep2prot,
        "psm_long": psm_long,
    }


def add_peptide_psm_count(df: pd.DataFrame, psm_index: dict) -> pd.DataFrame:
    df = df.copy()
    pep2psm = psm_index["pep2psm"]

    df["pep_key"] = df["peptide"].map(pep_key)
    df["psm_count"] = df["pep_key"].map(lambda k: int(pep2psm.get(k, 0)))

    return df


# ============================================================
# SUPPORT EXCLUDING CURRENT PEPTIDE
# ============================================================

def add_intra_inter_support(
    df: pd.DataFrame,
    psm_index: dict,
    top_n_list: int = 20,
) -> pd.DataFrame:
    df = df.copy()

    pep2prot = psm_index["pep2prot"]
    pep2psm = psm_index["pep2psm"]

    def row_calc(row):
        group_set = group_set_full_ids(row)

        if not group_set:
            return pd.Series(
                [0, 0, [], [], []],
                index=[
                    "psm_intra",
                    "psm_inter",
                    "pep_intra_list",
                    "pep_inter_list",
                    "group_aliases_used",
                ],
            )

        current_k = row.get("pep_key", pep_key(row.get("peptide", "")))

        intra = []
        inter = []

        for k, prots in pep2prot.items():
            if k == current_k:
                continue

            if not (prots & group_set):
                continue

            n_psm = int(pep2psm.get(k, 0))

            prots_base = {
                a for a in prots
                if ENSP_RE.fullmatch(a) or ENST_RE.fullmatch(a)
            }

            group_base = {
                a for a in group_set
                if ENSP_RE.fullmatch(a) or ENST_RE.fullmatch(a)
            }

            if prots_base and group_base:
                outside = prots_base - group_base
            else:
                outside = prots - group_set

            if len(outside) == 0:
                intra.append((k, n_psm))
            else:
                inter.append((k, n_psm, len(outside)))

        psm_intra = sum(x[1] for x in intra)
        psm_inter = sum(x[1] for x in inter)

        intra_sorted = sorted(intra, key=lambda t: t[1], reverse=True)[:top_n_list]
        inter_sorted = sorted(inter, key=lambda t: t[1], reverse=True)[:top_n_list]

        return pd.Series(
            [
                psm_intra,
                psm_inter,
                intra_sorted,
                inter_sorted,
                sorted(group_set),
            ],
            index=[
                "psm_intra",
                "psm_inter",
                "pep_intra_list",
                "pep_inter_list",
                "group_aliases_used",
            ],
        )

    df[
        [
            "psm_intra",
            "psm_inter",
            "pep_intra_list",
            "pep_inter_list",
            "group_aliases_used",
        ]
    ] = df.apply(row_calc, axis=1)

    df["total_support_excl_current"] = df["psm_intra"] + df["psm_inter"]

    return df


# ============================================================
# FILTER
# ============================================================

def apply_singleton_filter(df: pd.DataFrame, reject_reason: str) -> pd.DataFrame:
    df = df.copy()

    df["support_decision"] = "KEEP"
    df["support_reject_reason"] = ""

    psm_count = (
        pd.to_numeric(df["psm_count"], errors="coerce")
        .fillna(0)
        .astype(int)
    )

    support = (
        pd.to_numeric(df["total_support_excl_current"], errors="coerce")
        .fillna(0)
        .astype(int)
    )

    mask = (psm_count == 1) & (support == 0)

    df.loc[mask, "support_decision"] = "REJECT"
    df.loc[mask, "support_reject_reason"] = reject_reason

    return df


# ============================================================
# ANALYSIS
# ============================================================

def analyse_one(config: dict, sample: str, category: str) -> pd.DataFrame:
    print(f"\n=== Processing {sample} | {category} ===")

    input_tsv = get_input_tsv(config, sample, category)
    custom_fasta = get_custom_fasta(config, sample)
    pkl_path = get_seq_to_ids_pkl(config, sample)
    psm_glob = get_psm_glob(config, sample, category)
    reject_reason = CATEGORIES[category]["reject_reason"]

    print(f"[PATH] input_tsv    : {input_tsv}")
    print(f"[PATH] custom_fasta : {custom_fasta}")
    print(f"[PATH] seq_to_ids   : {pkl_path}")
    print(f"[PATH] psm_glob     : {psm_glob}")

    for p in [input_tsv, custom_fasta, pkl_path]:
        if not p.exists():
            raise FileNotFoundError(f"Fichier introuvable: {p}")

    print(f"[1] Loading input: {input_tsv}")
    df = pd.read_csv(input_tsv, sep="\t")

    if "peptide" not in df.columns:
        pep_col = find_peptide_col(df)
        df = df.rename(columns={pep_col: "peptide"})

    if "protein" not in df.columns:
        prot_col = find_protein_col(df)
        df = df.rename(columns={prot_col: "protein"})

    df["cell_line"] = sample
    df["sample"] = sample
    df["category"] = category

    if "row_id" not in df.columns:
        df.insert(0, "row_id", np.arange(len(df)))

    print(f"    rows before dedup: {len(df)}")

    df = df.drop_duplicates(
        subset=["sample", "category", "peptide"],
        keep="first",
    ).copy()

    print(f"    rows after dedup: {len(df)}")

    print("[2] Cleaning protein column")
    df["protein"] = df["protein"].apply(to_list)
    df["protein"] = df["protein"].apply(clean_protein_list)

    print(f"[3] Loading FASTA: {custom_fasta}")
    id2seq = load_fasta_id2seq(custom_fasta)

    print(f"[4] Loading pickle: {pkl_path}")
    seq_to_ids = load_seq_to_ids_pkl(pkl_path)

    print("[5] Adding TA")
    df = add_TA_column(df, id2seq=id2seq, seq_to_ids=seq_to_ids)

    print("[6] Classifying unique/group")
    df = classify_unique_group(df)

    print("[7] Loading PSM files")
    psm_files = load_psm_files(psm_glob)

    print(f"    n PSM files: {len(psm_files)}")
    for p in psm_files:
        print(f"      - {p}")

    print("[8] Building PSM index")
    psm_index = build_psm_index(psm_files)

    print("[9] Counting current peptide PSMs")
    df = add_peptide_psm_count(df, psm_index)

    print("[10] Computing support excluding current peptide")
    df = add_intra_inter_support(df, psm_index)

    print("[11] Applying singleton unsupported filter")
    df = apply_singleton_filter(df, reject_reason=reject_reason)

    print("[SUMMARY]")
    print(df["support_decision"].value_counts(dropna=False))
    print(df["cat"].value_counts(dropna=False))

    return df


def run_peptide_support_analysis(config: dict) -> dict[str, Path]:
    """
    Étape 6 — Analyse du support peptidique.

    Boucle sur les samples dans config["samples"] et sur les catégories :
        - canonical_variant
        - nc_variant
        - nc_wt

    Sorties par sample :
        analysis_outdir/{sample}/06_peptide_support/
            {sample}_{category}_support.tsv
            ALL_categories_support.tsv
            summary_counts.tsv
    """

    outputs: dict[str, Path] = {}

    print("\n" + "#" * 100)
    print("[Étape 6] Analyse du support peptidique")
    print("#" * 100)

    for sample in config["samples"]:
        print("\n" + "=" * 100)
        print(f"[SAMPLE] {sample}")
        print("=" * 100)

        outdir = get_support_outdir(config, sample)
        outdir.mkdir(parents=True, exist_ok=True)

        print(f"[OUTPUT] {outdir}")

        sample_dfs = []

        for category in CATEGORIES:
            df = analyse_one(config, sample, category)

            out_file = outdir / f"{sample}_{category}_support.tsv"
            df.to_csv(out_file, sep="\t", index=False)

            print(f"[WRITE] {out_file}")

            outputs[f"{sample}_{category}"] = out_file
            sample_dfs.append(df)

        if not sample_dfs:
            print(f"[WARN] Aucun résultat produit pour {sample}")
            continue

        sample_all_df = pd.concat(sample_dfs, ignore_index=True)

        out_all = outdir / "ALL_categories_support.tsv"
        sample_all_df.to_csv(out_all, sep="\t", index=False)

        print(f"[WRITE] {out_all}")

        summary = (
            sample_all_df
            .groupby(
                ["sample", "category", "cat", "support_decision"],
                dropna=False,
            )
            .size()
            .reset_index(name="n")
        )

        out_summary = outdir / "summary_counts.tsv"
        summary.to_csv(out_summary, sep="\t", index=False)

        print(f"[WRITE] {out_summary}")

        print("\n=== FINAL SUMMARY ===")
        print(summary)

        outputs[f"{sample}_all_categories"] = out_all
        outputs[f"{sample}_summary_counts"] = out_summary

    if not outputs:
        raise RuntimeError("Aucun fichier produit par l'analyse du support peptidique.")

    print("\n[DONE] Étape 6 — Analyse du support peptidique terminée.")

    return outputs


if __name__ == "__main__":
    ()