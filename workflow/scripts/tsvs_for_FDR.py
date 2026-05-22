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


df["pep"] = 1 - df["peptideprophet_probability"]

df = auxiliary.target_decoy.qvalues(
    df, pep="pep", q_label="q_value", full_output=True
)

def double_check(peptide, check_list): # Make sur non canonical isn't anywhere in a canonical prot or contam
    for prot in check_list:
        if peptide in prot:
            return 'NOT UNIQUE' # We don't want thos anywhere because they are not semi-tryptic for a canonical but the sequence is in a canonical prot so can't include in non canonical
    return 'NONCANONICAL'

def get_type(prots, pep, check_prots):
    if any(prot.startswith('CONTAM') for prot in prots):
        return 'CONTAM'
    elif any(prot.startswith('CANONICAL') for prot in prots):
        return 'CANONICAL'
    elif any(prot.startswith('NONCANONICAL') for prot in prots):
        return double_check(pep, check_prots)


check_prot = pkl.load(open(CHECK, 'rb')) # check file containing canonical and contam prots

print('CHECK START')

pep_assignations = {} # Dictionnary with(pep identified : pep assigned to it) The assigned represent the complete peptide if the identified one is semi-tryptic

df['TYPE'] = None
df['ASSIGNED'] = None

for i, peptide in df.iterrows():
    prots = peptide['protein']
    pep = peptide['peptide']
    if any('TARGET' in prot for prot in prots): # Look at all the possible peptides if there is a target
        valid_prots = [prot for prot in prots if 'TARGET'] # Get all the targets of the correct type from which the peptide is coming 
        pep_type = get_type(valid_prots, pep, check_prot) # Look for the type of this target
        if pep_type not in ['CANONICAL', 'NONCANONICAL']: # If not in those categories, no need to assigned a peptide
            df.at[i,'TYPE'] = pep_type
            df.at[i,'ASSIGNED'] = 'NA'
            continue
        
        all_seq = [prot.split('_')[3] for prot in valid_prots if pep_type in prot]
        shortest_seq = min(all_seq, key=len) #Select the shortest sequence (closest to the identified peptide) for identification of the right category by FDR bench
        pep_assignations[pep] = shortest_seq # Save the assigned peptide in a pkl

        df.at[i,'TYPE'] = pep_type
        df.at[i,'ASSIGNED'] = shortest_seq

    
    elif any('ENTRAPMENT' in prot for prot in prots): # Same with all entrapment peptides
        valid_prots = [prot for prot in prots if 'ENTRAPMENT']
        pep_type = get_type(valid_prots, pep, check_prot)
        if pep_type not in ['CANONICAL', 'NONCANONICAL']: # If not in those categories, no need to assigned a peptide
            df.at[i,'TYPE'] = pep_type
            df.at[i,'ASSIGNED'] = 'NA'
            continue

        all_seq = [prot.split('_')[3] for prot in valid_prots if pep_type in prot]
        shortest_seq = min(all_seq, key=len)
        pep_assignations[pep] = shortest_seq

        df.at[i,'TYPE'] = pep_type
        df.at[i,'ASSIGNED'] = shortest_seq

    else: # The rest are decoys exclude from both TSV
        df.at[i,'TYPE'] = 'DECOY'
        df.at[i,'ASSIGNED'] = pep

pkl.dump(pep_assignations, open(f"{OUTPATH}assigned_peptides.pkl", "wb"))

all_type = df['TYPE'].unique()
if None in all_type:
    print('WARNING some peptide types not assigned')

all_assigned = df['ASSIGNED'].unique()
if None in all_type:
    print('WARNING some peptide no sequence assigned')

print('DIVIDING FOR TSV')

def clean_df(df, outfile,out_psm=None):
    copied = df[['ASSIGNED','peptideprophet_probability','q_value']].copy()
    copied.rename(columns={'peptideprophet_probability': 'score'}, inplace=True)
    copied.rename(columns={'ASSIGNED': 'peptide'}, inplace=True)

    dup = copied[copied['peptide'].duplicated()==True]

    if out_psm is not None:
        copied.to_csv(out_psm, sep="\t", index=False)

    sorted_df_for_tsv = copied.sort_values(by='score', ascending=False)
    df_no_duplicates = sorted_df_for_tsv.drop_duplicates(subset=['peptide'], keep='first')

    if len(df_no_duplicates) + len(dup) == len(copied):
        df_no_duplicates.to_csv(outfile, sep='\t', index=False)
    else:
        print('WARNING element missing')

canonicals = df[df['TYPE'] == 'CANONICAL']
non_canonicals = df[df['TYPE'] == 'NONCANONICAL']


clean_df(canonicals,
         f"{OUTPATH}canonicals.tsv",
         out_psm=f"{OUTPATH}canonicals_psm.tsv")

clean_df(non_canonicals,
         f"{OUTPATH}non_canonicals.tsv",
         out_psm=f"{OUTPATH}non_canonicals_psm.tsv")
