import argparse 

parser = argparse.ArgumentParser()
parser.add_argument("--output", required=True)
parser.add_argument("--log", required=True)
args = parser.parse_args()
FASTA = args.output

fasta_dict = {}
order = []

with open(FASTA, "r") as f:
    for line in f:
        if line.startswith('>'):
            id = line.strip()
            order.append(id)
        else:
            fasta_dict[id] = line.strip()

with open(FASTA, "w") as f:
    for id in order:
        seq = fasta_dict[id]
        if id.startswith(">DECOY_") and seq.startswith(("K", "R")):
            new_seq = seq[1:] + seq[0]
            fasta_dict[id] = new_seq
            
        f.writelines([id, '\n', fasta_dict[id], '\n'])
