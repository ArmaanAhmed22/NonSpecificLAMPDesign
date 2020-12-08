import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

#Open reference to get sequence length
with open(snakemake.input[0],"r") as handle:
	reference = next(SeqIO.parse(handle,"fasta")).seq


#Open quasispecies file
sequences = list(SeqIO.parse(open(snakemake.input[1],"r"),"fasta"))

#dict storing nucleotide prevalences
prevalences = [{"A":0,"G":0,"C":0,"T":0} for _ in range(len(reference))]

#Create nucleotide prevalence. (Note to User: Make sure there aren't any insertions/deletions/ambiguous nucleotides!)
for seq_ref in tqdm(sequences):

	for nuc_pos in range(len(prevalences)):
		prevalences[nuc_pos][seq_ref.seq[nuc_pos]]+=1

#Output to file
pd.DataFrame(prevalences).to_csv(snakemake.output[0])
