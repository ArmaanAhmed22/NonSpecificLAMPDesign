from Bio import SeqIO
from Bio.Alphabet import generic_dna

with open(snakemake.input[0]) as handle:
	reference = SeqIO.parse(handle,"fasta",generic_dna)[0].seq
	
with open(snakemake.input[1]) as handle:
	sequences = SeqIO.parse(handle,"fasta",generic_dna)