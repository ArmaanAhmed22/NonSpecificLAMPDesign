import pandas as pd
from Bio import SeqIO
import json
from tqdm import tqdm

params = json.load(open(snakemake.input[2],"r"))
ambiguous_to_nt = json.load(open(snakemake.input[3],"r"))
df_primers = pd.read_csv(snakemake.input[0])

targets = [0 for _ in range(df_primers.shape[0])]

num_to_region = {
	0:"F3",
	1:"F2",
	2:"F1c",
	3:"B1c",
	4:"B2",
	5:"B3"
}


def target_sequence(sequence,primers):

	for primer_num,primer in enumerate(primers):
		mismatch = 0
		for nuc_relative_pos,nuc in enumerate(primer["sequence"]):
			nuc_pos = nuc_relative_pos + primer["position"]
			#SHOULD CHANGE: MISMATCH TOLERANCE DEPENDS ON REGION + primer
			if mismatch > params["mismatch_tolerance"]["tolerance"]:
				return False
			if nuc in ambiguous_to_nt.keys():
				if not (sequence[nuc_pos] in ambiguous_to_nt[nuc]):
					mismatch+=1
			elif nuc != sequence[nuc_pos]:
				mismatch += 1
	return True

total = 0

untargetable = [[] for _ in range(df_primers.shape[0])]
for record in tqdm(SeqIO.parse(open(snakemake.input[1],"r"),"fasta")):
	sequence = str(record.seq)
	total+=1
	targeted = False
	for i,row in df_primers.iterrows():
		if target_sequence(sequence,[{"sequence" : row[f"{i}_sequence"], "position" : row[f"{i}_position"]} for i in range(6)]):
			targets[i]+=1
			targeted = True
	if not targeted:
		untargetable[i].append(record)

for target in targets:
	print(target/total)

SeqIO.write(untargetable[0], open(snakemake.output[1],"w"), "fasta")
open(snakemake.output[0],"w")

