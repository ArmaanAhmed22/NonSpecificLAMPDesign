from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from tqdm import tqdm

non_amb_nucs = ["A","G","C","T"]

reference = next(SeqIO.parse(open(snakemake.input[0],"r"),"fasta")).seq
sequences = list(SeqIO.parse(open(snakemake.input[1],"r"),"fasta"))

aligner = PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 1
aligner.mismatch_score = -2
aligner.gap_score = -2

sequence_output = []

for seq_ref in tqdm(sequences):
	seq = seq_ref.seq
	seq_str = str(seq_ref.seq).upper()

	not_good = False
	for nuc in seq_str:
		if nuc not in non_amb_nucs:
			not_good = True
			break
	if not_good:
		continue

	seq_ref.seq = Seq(seq_str)

	if len(seq) != len(reference):
		continue
	alignment = aligner.align(reference, seq)[0]

	if "-" in alignment.target or "-" in alignment.query:
		continue

	sequence_output.append(seq_ref)

SeqIO.write(sequence_output,open(snakemake.output[0],"w"),"fasta")

