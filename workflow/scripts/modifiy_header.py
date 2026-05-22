import sys

sys.stdout = open(snakemake.log[0], "w")
sys.stderr = sys.stdout

TXT = snakemake.input["txt"]
FINAL = snakemake.output["fasta_entrapment_headers"]

new_headers = {}

with open(TXT, 'r') as f:
    for line in f:
        if line.startswith('sequence'):
            continue
        line = line.strip().split('\t')
        pep = line[0]
        prots = line[2].split(';')
        pep_type = line[3]
        index = line[4]

        if any(prot.startswith('CONTAM') for prot in prots):
            cat = 'CONTAM'
        elif any(prot.startswith('ENSP')  for prot in prots):
            cat = 'CANONICAL'
        elif any(prot.startswith('ENST')  for prot in prots):
            cat = 'NONCANONICAL'
        else:
            print('WARNING NO TYPE FOUND FOR ', line)

        if str(pep_type) == 'p_target':
            entrapment = "ENTRAPMENT"
        elif str(pep_type) == 'target':
            entrapment = "TARGET"
        else:
            print('WARNING ENTRAPMENT NOT FOUND FOR ', line, str(pep_type))

        new_name = '_'.join([cat, entrapment, index, pep])
        new_headers[pep] = new_name

print('QC')
print(len(list(new_headers.values())), len(set(new_headers.values())))

with open(FINAL, 'w') as f:
    for pep, header in new_headers.items():
        f.writelines(['>', header, '\n', pep, '\n'])
