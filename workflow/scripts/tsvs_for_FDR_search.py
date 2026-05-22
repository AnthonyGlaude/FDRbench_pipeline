import pandas as pd
from pyteomics import pepxml, auxiliary
import pickle as pkl
import sys

sys.stdout = open(snakemake.log[0], "w")
sys.stderr = sys.stdout

XML = snakemake.input["pepxml"]
CHECK = snakemake.input["check_prot"]
OUTPATH = snakemake.params["outpath"]

df = pepxml.DataFrame(XML)
print("=== DEBUG pepXML -> DataFrame ===")
print("XML =", XML)
print("Rows =", len(df))
print("Columns =", list(df.columns))

for idx in range(min(3, len(df))):
    print(f"\n--- Row {idx} ---")
    print("peptide =", df.iloc[idx].get("peptide"))
    print("prob   =", df.iloc[idx].get("peptideprophet_probability"))
    print("protein=", df.iloc[idx].get("protein"))

df["pep"] = 1 - df["peptideprophet_probability"]

df = auxiliary.target_decoy.qvalues(
    df, pep="pep", q_label="q_value", full_output=True
)

def double_check(peptide, check_list): # Make sur non canonical isn't anywhere in a canonical prot or contam
    for prot in check_list:
        if peptide in prot:
            return 'NOT UNIQUE'
    return 'NONCANONICAL'

def get_type(prots, pep, check_prots):
    if any(p.startswith('CONTAM') for p in prots):
        return 'CONTAM'
    if any(p.startswith('ENSP') or p.startswith('CANONICAL') for p in prots):
        return 'CANONICAL'
    if any(p.startswith('ENST') or p.startswith('NONCANONICAL') for p in prots):
        return double_check(pep, check_prots)
    return None

# --- PATCH: helpers for "real search" mode ---
def is_decoy_protein(p):
    return p.startswith('DECOY_') or p.startswith('rev_') or p.startswith('DECOY')

# detect whether this pepXML looks like entrapment-style IDs
has_target_token = False
has_entrapment_token = False
for prots in df['protein'].head(200):
    if any('TARGET' in p for p in prots):
        has_target_token = True
    if any('ENTRAPMENT' in p for p in prots):
        has_entrapment_token = True
    if has_target_token or has_entrapment_token:
        break

print("MODE:", "ENTRAPMENT" if (has_target_token or has_entrapment_token) else "SEARCH_NORMAL")
# --- END PATCH ---

check_prot = pkl.load(open(CHECK, 'rb'))
print('CHECK START')

pep_assignations = {}
df['TYPE'] = None
df['ASSIGNED'] = None

for i, peptide in df.iterrows():
    prots = peptide['protein']
    pep = peptide['peptide']

    # --- PATCH: if no TARGET/ENTRAPMENT tokens, handle as real search ---
    if not (has_target_token or has_entrapment_token):
        if any(is_decoy_protein(p) for p in prots):
            df.at[i,'TYPE'] = 'DECOY'
            df.at[i,'ASSIGNED'] = pep
            continue

        pep_type = get_type(prots, pep, check_prot)
        df.at[i,'TYPE'] = pep_type
        df.at[i,'ASSIGNED'] = pep  # keep it simple for search-normal
        continue
    # --- END PATCH ---

    if any('TARGET' in prot for prot in prots):
        valid_prots = [prot for prot in prots if 'TARGET' in prot]  # PATCH tiny bug: add 'in prot'
        pep_type = get_type(valid_prots, pep, check_prot)
        if pep_type not in ['CANONICAL', 'NONCANONICAL']:
            df.at[i,'TYPE'] = pep_type
            df.at[i,'ASSIGNED'] = 'NA'
            continue

        all_seq = [prot.split('_')[3] for prot in valid_prots if pep_type in prot]
        shortest_seq = min(all_seq, key=len)
        pep_assignations[pep] = shortest_seq

        df.at[i,'TYPE'] = pep_type
        df.at[i,'ASSIGNED'] = shortest_seq

    elif any('ENTRAPMENT' in prot for prot in prots):
        valid_prots = [prot for prot in prots if 'ENTRAPMENT' in prot]  # PATCH tiny bug: add 'in prot'
        pep_type = get_type(valid_prots, pep, check_prot)
        if pep_type not in ['CANONICAL', 'NONCANONICAL']:
            df.at[i,'TYPE'] = pep_type
            df.at[i,'ASSIGNED'] = 'NA'
            continue

        all_seq = [prot.split('_')[3] for prot in valid_prots if pep_type in prot]
        shortest_seq = min(all_seq, key=len)
        pep_assignations[pep] = shortest_seq

        df.at[i,'TYPE'] = pep_type
        df.at[i,'ASSIGNED'] = shortest_seq

    else:
        df.at[i,'TYPE'] = 'DECOY'
        df.at[i,'ASSIGNED'] = pep

pkl.dump(pep_assignations, open(f"{OUTPATH}assigned_peptides.pkl", "wb"))

print("\n=== DEBUG TYPE counts ===")
print(df["TYPE"].value_counts(dropna=False))

print('DIVIDING FOR TSV')

def clean_df(df_in, outfile, out_psm=None):
    # colonnes de base qu'on veut toujours
    base = ["peptide", "ASSIGNED", "peptideprophet_probability", "q_value", "TYPE", "protein"]

    # colonnes PSM utiles (présentes selon pepXML)
    extras = [
        "spectrum", "spectrumNativeID",
        "start_scan", "end_scan", "scan",
        "assumed_charge", "charge",
        "retention_time_sec",
        "precursor_neutral_mass", "precursor_mz"
    ]

    keep = [c for c in base if c in df_in.columns] + [c for c in extras if c in df_in.columns]

    copied = df_in[keep].copy()

    # uniformiser / renommer
    copied.rename(columns={
        "peptideprophet_probability": "score",
        "start_scan": "scan",   # si start_scan existe, on l'appelle scan
        "charge": "assumed_charge"  # si jamais c'est "charge" chez vous
    }, inplace=True)

    # Option: extraire un "run" depuis spectrum (souvent: RUN.1234.1234.2)
    if "spectrum" in copied.columns and "run" not in copied.columns:
        copied["run"] = copied["spectrum"].astype(str).str.replace(r"\..*$", "", regex=True)

    # 1) écrire le fichier PSM complet (toutes les lignes)
    if out_psm is not None:
        copied.to_csv(out_psm, sep="\t", index=False)

    # 2) écrire le fichier "unique peptides": meilleur score par peptide
    if "peptide" in copied.columns:
        sorted_df = copied.sort_values(by="score", ascending=False)
        df_no_duplicates = sorted_df.drop_duplicates(subset=["peptide"], keep="first")
        df_no_duplicates.to_csv(outfile, sep="\t", index=False)
    else:
        # fallback si jamais "peptide" n'existe pas (rare)
        copied.to_csv(outfile, sep="\t", index=False)


canonicals = df[df['TYPE'] == 'CANONICAL']
non_canonicals = df[df['TYPE'] == 'NONCANONICAL']

clean_df(canonicals,
         f"{OUTPATH}canonicals.tsv",
         out_psm=f"{OUTPATH}canonicals_psm.tsv")

clean_df(non_canonicals,
         f"{OUTPATH}non_canonicals.tsv",
         out_psm=f"{OUTPATH}non_canonicals_psm.tsv")
