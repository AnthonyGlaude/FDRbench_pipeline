import pandas as pd
from pyteomics import pepxml, auxiliary
import pickle as pkl
import math
import sys

sys.stdout = open(snakemake.log[0], "w")
sys.stderr = sys.stdout

PROPHET = snakemake.input["prophet"]
FDR_C = snakemake.input["c_FDR"]
FDR_NC = snakemake.input["nc_FDR"]
CHECK = snakemake.input["check_prot"]
ENT = snakemake.input["entrapment_pep"]
OUTPATH = snakemake.params["outpath"]
SEEDS = snakemake.params["seeds"]

print("FDR empirical threshold")


# To remove entrapment peptide for the ranking
def remove_entrapment(df, entrapment_pep):
    df['is_entrapment'] = None
    for i, line in df.iterrows():
        peptide=line['peptide']
        if peptide in entrapment_pep:
            df.at[i,'is_entrapment'] = True
        else:
            df.at[i,'is_entrapment'] = False
    return df

def get_max_rank(FDR_df):
    sorted_df = FDR_df.sort_values(by='score', ascending=False)
    sorted_df['Rank'] = sorted_df['score'].rank(method='dense',ascending=False)
    first_index = sorted_df[sorted_df["paired_fdp"] > 0.01].index[0]
    fdp_df = sorted_df[:first_index]
    if len(fdp_df) == 0:
        return 0
    max_rank = max(fdp_df['Rank'])
    return max_rank

ranks = {'CAN':[], 'NC':[]}
for seed in SEEDS:
    print(f"Starting seed:{seed}")
    index = seed-1

    print(FDR_C[index])
    print(FDR_NC[index])
    print(ENT[index])

    can_FDR_df = pd.read_csv(FDR_C[index], sep=',')
    nc_FDR_df = pd.read_csv(FDR_NC[index], sep=',')
    entrapment_pep = pkl.load(open(ENT[index], 'rb'))

    ent_can_FDR = remove_entrapment(can_FDR_df, entrapment_pep)
    ent_nc_FDR = remove_entrapment(nc_FDR_df, entrapment_pep)

    clean_can_FDR = ent_can_FDR[ent_can_FDR['is_entrapment'] == False]
    clean_nc_FDR = ent_nc_FDR[ent_nc_FDR['is_entrapment'] == False]

    can_max_rank = get_max_rank(clean_can_FDR)
    nc_max_rank = get_max_rank(clean_nc_FDR)
    ranks['CAN'].append(can_max_rank)
    ranks['NC'].append(nc_max_rank)

print('CALCULATING THE MEAN')

CAN_RANK = int(math.floor(sum(ranks['CAN'])/len(ranks['CAN'])))
NC_RANK = int(math.floor(sum(ranks['NC'])/len(ranks['NC'])))

print(CAN_RANK, NC_RANK)

print("No Entrapment Filtration")

prophet_df = pepxml.DataFrame(PROPHET)
check_prot = pkl.load(open(CHECK, 'rb'))

def check_prot_seq(can_list, peptide):
    for prot in can_list:
        if peptide in prot:
            return True
    return False


prophet_df['class'] = None
for i, line in prophet_df.iterrows():
    if any(item.startswith('CONTAM') for item in line['protein']):
        prophet_df.at[i,'class'] = 'contam'
    elif any(item.startswith('ENSP') for item in line['protein']):
        prophet_df.at[i,'class'] = 'canonical'
    elif any(item.startswith('ENST') for item in line['protein']):
        peptide = line['peptide']
        check = check_prot_seq(check_prot, peptide)
        if check == True:
            prophet_df.at[i,'class'] = 'canonical'
        elif check == False:
            prophet_df.at[i,'class'] = 'non-canonical'
    elif any(item.startswith('DECOY_') for item in line['protein']):
            prophet_df.at[i,'class'] = 'decoy'
    else:
        print('WARNING NO ASSIGNED PROTEIN')

prophet_decoy_df = prophet_df[prophet_df["class"] == 'decoy']
prophet_contam_df = prophet_df[prophet_df["class"] == 'contam']
prophet_can_df = prophet_df[prophet_df["class"] == 'canonical']
prophet_no_can_df = prophet_df[prophet_df["class"] == 'non-canonical']

missing_pep = len(prophet_df) - (len(prophet_decoy_df) + len(prophet_contam_df) + len(prophet_can_df) + len(prophet_no_can_df))

print(f"{missing_pep} missing peptides | {PROPHET}")
print(f"{len(prophet_decoy_df)} decoy peptides | {PROPHET}")
print(f"{len(prophet_contam_df)} contam peptides | {PROPHET}")
if missing_pep > 0:
    print('WARNING WARNING WARNING')


def get_final_df(PROPHET_df, max_rank):
    sorted_df= PROPHET_df.sort_values(by='peptideprophet_probability', ascending=False)
    df_no_duplicates = sorted_df.drop_duplicates(subset=['peptide'], keep='first')
    df_no_duplicates['Rank'] = df_no_duplicates['peptideprophet_probability'].rank(method='dense',ascending=False)
    final_df = df_no_duplicates[df_no_duplicates['Rank'] <= max_rank]
    return final_df


get_final_df(prophet_can_df, CAN_RANK).to_csv(f"{OUTPATH}canonical_final.tsv", sep='\t', index=False)
get_final_df(prophet_no_can_df, NC_RANK).to_csv(f"{OUTPATH}nc_final.tsv", sep='\t', index=False)
