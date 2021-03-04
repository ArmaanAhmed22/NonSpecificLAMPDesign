from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def reverse_sequences():
	reference_rev = [SeqRecord(reference.seq.reverse_complement(),id=reference.id,name=reference.name,description=reference.description)]
	quasispecies_rev = [SeqRecord(record.seq.reverse_complement(),id=record.id,name=record.name,description=record.description) for record in quasispecies]
	return (reference_rev,quasispecies_rev)

reference = next(SeqIO.parse(open(snakemake.input[0],"r"),"fasta"))
quasispecies = SeqIO.parse(open(snakemake.input[1],"r"),"fasta")

reference_rev,quasispecies_rev = reverse_sequences()
SeqIO.write(reference_rev,open(snakemake.output[0],"w"),"fasta")
SeqIO.write(quasispecies_rev,open(snakemake.output[1],"w"),"fasta")