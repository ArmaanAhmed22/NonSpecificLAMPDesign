import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import json
from tqdm import tqdm
import primer3
from collections import Counter

thermo_imp = __import__("thermo")

params = json.load(open(snakemake.input[2],"r"))
ambiguous_to_nt = json.load(open(snakemake.input[3],"r"))
df_primers = pd.read_csv(snakemake.input[0])
assay_params = json.load(open(snakemake.input[4],"r"))
thermo_params = json.load(open(snakemake.input[5],"r"))

targets = [0 for _ in range(df_primers.shape[0])]

num_to_region = {
	0:"F3",
	1:"F2",
	2:"F1c",
	3:"B1c",
	4:"B2",
	5:"B3"
}

def inEnd(nuc_pos,end_len,rest_len,region_name):
	if region_name == "F3" or region_name == "F2":
		return not (nuc_pos < rest_len)
	elif region_name == "B1c":
		return (nuc_pos < end_len)
	elif region_name == "F1c":
		return not (nuc_pos < rest_len)
	elif region_name == "B2" or region_name == "B3":
		return (nuc_pos < end_len)
	else:
		raise Exception("Not legal!!")

def target_sequence(sequence,primers):

	for primer_num,primer in enumerate(primers):
		mismatch = 0
		region_name = num_to_region[primer_num]
		end_region_len = params["length"][region_name]["minimum"] // 3
		rest_len = len(primer["sequence"]) - end_region_len
		for nuc_relative_pos,nuc in enumerate(primer["sequence"]):
			nuc_pos = nuc_relative_pos + primer["position"]
			#SHOULD CHANGE: MISMATCH TOLERANCE DEPENDS ON REGION + primer
			
			if mismatch > params["mismatch_tolerance"]["tolerance"]:
				return False
			if nuc in ambiguous_to_nt.keys():
				if not (sequence[nuc_pos] in ambiguous_to_nt[nuc]):
					mismatch+=1
					if inEnd(nuc_relative_pos,end_region_len,rest_len,region_name):
						return False
			elif nuc != sequence[nuc_pos]:
				mismatch += 1
				if inEnd(nuc_relative_pos,end_region_len,rest_len,region_name):
					return False
	return True

def good_thermo(sequence,primers):
	reg_seq = {num_to_region[i] : get_primers_spread(primer["sequence"]) for i,primer in enumerate(primers)}
	for region_name,sequence_var in reg_seq.items():
		primer3_params_end = {
		"dv_conc":assay_params["DV_CONC"],
		"dntp_conc":assay_params["DNTP_CONC"],
		"dna_conc":assay_params[region_name],
		"temp_c":assay_params["TEMP_C"]
		}
		primer3_params_tm = {
		"dv_conc":assay_params["DV_CONC"],
		"dntp_conc":assay_params["DNTP_CONC"],
		"dna_conc":assay_params[region_name],
		}
		good = False
		for sequence in sequence_var:
			res = thermo_imp.calc_median_gc(region_name,Counter([Seq(sequence)]))
			if thermo_params["GC"][region_name]["minimum"] < res and res < thermo_params["TM"][region_name]["maximum"]:
				good = True
				break
		if not good:
			return False
	return True



def get_primers_spread(primer):
	spread = [""]
	for nuc in primer:
		#l_nucs = ambiguous_to_nt[nuc]
		if not (nuc in ["A","G","C","T"]):
			l_nucs = ambiguous_to_nt[nuc]
			spread = spread*len(l_nucs)
			for i,l_nuc in enumerate(l_nucs):
				for j in range(i*len(spread)//len(l_nucs),(i+1)*len(spread)//len(l_nucs)):
					spread[j]+=l_nuc
		else:
			for i in range(len(spread)):
				spread[i]+=nuc
	return spread

total = 0

untargetable = [[] for _ in range(df_primers.shape[0])]
for record in tqdm(SeqIO.parse(open(snakemake.input[1],"r"),"fasta")):
	sequence = str(record.seq)
	total+=1
	targeted = False
	for i,row in df_primers.iterrows():
		primer_data = [{"sequence" : row[f"{i}_sequence"], "position" : row[f"{i}_position"]} for i in range(6)]
		if target_sequence(sequence,primer_data):
			targets[i]+=1
			targeted = True
	if not targeted:
		untargetable[i].append(record)

with open(snakemake.output[0],"w") as h:
	for i,target in enumerate(targets):
		h.write(f"{target/total}")
		if i != len(targets)-1:
			h.write("\n")


SeqIO.write(untargetable[0], open(snakemake.output[1],"w"), "fasta")

