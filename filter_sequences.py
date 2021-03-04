import argparse
from Bio import SeqIO
from tqdm import tqdm

def parser():
	args = argparse.ArgumentParser()
	args.add_argument("--fasta","-f",help="Input fasta")
	return args.parse_args()

input_file_name = parser().fasta

qs = list(SeqIO.parse(open(input_file_name,"r"),"fasta"))

def checkFrag(sequence):
	return sequence[0] == "X" or sequence[-1] == "X"

normal = "ACGT"
good_seqs = []
for record in tqdm(qs):
	if checkFrag(record.seq):
		continue
	bad = False
	for nuc in record.seq:
		if not nuc in normal:
			bad = True
			break
	if bad:
		continue
	good_seqs.append(record)

print(f"OLD LENGTH: {len(qs)}\nNEW LENGTH: {len(good_seqs)}\nPERCENT PERSERVED: {len(good_seqs)/len(qs)}")
SeqIO.write(good_seqs, open(input_file_name.split(".")[0]+"_filtered.fasta","w"), "fasta")
