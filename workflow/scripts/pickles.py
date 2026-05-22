import sys

sys.stdout = open(snakemake.log[0], "w")
sys.stderr = sys.stdout

import pickle as pkl
from Bio import SeqIO

SEARCH = snakemake.input["search"]
ENTRAPMENT = snakemake.input["entrapment"]
OUTPATH = snakemake.params["outpath"]

entrapment = []

print("ENTRAPMENT DATABASE")
fasta = SeqIO.to_dict(SeqIO.parse(ENTRAPMENT, "fasta"))
for prot in fasta:
    if 'ENTRAPMENT' in prot:
        entrapment.append(str(fasta[prot].seq))

print(len(fasta), len(entrapment))

pkl.dump(entrapment, open(f"{OUTPATH}entrapment_pep.pkl", "wb"))

print("PARSING FASTA")

fasta = SeqIO.to_dict(SeqIO.parse(SEARCH, "fasta"))

check_prot = []
for prot in fasta:
    if prot.startswith('ENSP'):
        check_prot.append(str(fasta[prot].seq))
    elif prot.startswith('CONTAM'):
        check_prot.append(str(fasta[prot].seq))

print(len(check_prot))

pkl.dump(check_prot, open(f"{OUTPATH}check_prot.pkl", "wb"))
