from Bio import SeqIO
from Bio.Seq import Seq
import json
from tqdm import tqdm
import pandas as pd
import numpy as np
import primer3
from collections import Counter
import math
thermo_imp = __import__("thermo")




def calc_base_pair_shannon(nuc_prevalence: pd.DataFrame):
	nuc_prevalence["total"] = nuc_prevalence["A"] + nuc_prevalence["G"] + nuc_prevalence["C"] + nuc_prevalence["T"]
	shannon = None
	for nuc in ["A","G","C","T"]:
		percent = nuc_prevalence[nuc] / nuc_prevalence["total"]
		if shannon is None:
			shannon = np.where(percent == 0,percent,-percent * np.log2(percent))
		else:
			shannon += np.where(percent == 0,percent,-percent * np.log2(percent))
	nuc_prevalence["shannon"] = shannon

def get_partial_sequences(region_pos,region_len) -> Counter:
	holding = []
	for record in quasispecies:
		region = record.seq[region_pos:region_pos+region_len]
		holding.append(region)
	return Counter(holding)

def score_length_regions(pos1: int,pos2: int,region_len: int,type_region: str):
	
	primer3_params_end = {
			"dv_conc":assay_params["DV_CONC"],
			"dntp_conc":assay_params["DNTP_CONC"],
			"dna_conc":assay_params[type_region],
			"temp_c":assay_params["TEMP_C"]
	}
	primer3_params_tm = {
		"dv_conc":assay_params["DV_CONC"],
		"dntp_conc":assay_params["DNTP_CONC"],
		"dna_conc":assay_params[type_region],
	}

	primer_end = params["length"][type_region]["minimum"] // 3
	primer_rest = region_len - primer_end
	out = []
	for i in tqdm(range(pos1,pos2 - (region_len -1))):
		part_seqs = get_partial_sequences(i,region_len)
		if type_region == "LoopF":
			entropy_end = nt_prev["shannon"][i:i+primer_end].sum()
			entropy = nt_prev["shannon"][i+primer_end:i+region_len].sum()
		elif type_region == "LoopB":
			entropy_end = nt_prev["shannon"][i+primer_rest:i+region_len].sum()
			entropy = nt_prev["shannon"][i:i+primer_rest].sum()
		else:
			raise TypeError("Incorrect")
		tm = thermo_imp.calc_median_tm(type_region,part_seqs,primer3_params_tm)
		#Calc End stability
		end = thermo_imp.calc_median_end_stability(type_region,part_seqs,primer3_params_end)
		#Calc GC
		gc = thermo_imp.calc_median_gc(type_region,part_seqs)
		#Calc hairpin ()
		hairpin = thermo_imp.calc_hairpin(type_region,part_seqs,primer3_params_end,0.1)

		out.append({"region":type_region,"pos":i,"length":region_len,"entropy":entropy,"entropy_end":entropy_end,"tm":tm,"end":end,"gc":gc,"hairpin":hairpin})
	return out

def score_region(region_name:str,pos1:int,pos2:int):
	minimum = params["length"][region_name]["minimum"]
	maximum = params["length"][region_name]["maximum"]
	data = []
	for length in tqdm(range(minimum,maximum+1),desc=region_name):
		data.extend(score_length_regions(pos1,pos2,length,region_name))
	return data


def score_with_ideal_min_max(minimum,maximum,ideal,elem,weight=1) -> float:
	if elem < minimum:
		elem_m = minimum - elem
		return 1/(weight*elem_m+1) - 1
	elif elem > maximum:
		elem_m = elem - maximum
		return 1/(weight*elem_m+1) - 1
	if elem < ideal:
		return 1/(ideal-minimum) * (elem - minimum)
	else:
		return 1/(ideal-maximum) * (elem - maximum)

def get_ambiguous(A,G,C,T,cutoff=0.1):
	total = A+G+C+T
	possible = []
	for nuc,prevalence in {"A":A,"G":G,"C":C,"T":T}.items():
		if prevalence / total >= cutoff:
			possible.append(nuc)
	if len(possible) == 0:
		#print("yes")
		return None
	if len(possible) == 1:
		return possible[0]
	for amb,nucs in ambiguous_to_nt.items():
		should_use_this_amb = True
		for cur_nuc in possible:
			if not (cur_nuc in nucs):
				should_use_this_amb = False
				break
		if not should_use_this_amb:
			continue
		return amb
	return None

def main() -> None:
	FLoop_start = primers.iloc[0]["1_position"] + len(primers.iloc[0]["1_sequence"])
	FLoop_end = primers.iloc[0]["2_position"]
	BLoop_start = primers.iloc[0]["3_position"] + len(primers.iloc[0]["3_sequence"])
	BLoop_end = primers.iloc[0]["4_position"]

	data = []
	data.extend(score_region("LoopF",FLoop_start,FLoop_end))
	data.extend(score_region("LoopB",BLoop_start,BLoop_end))

	df = pd.DataFrame(data)
	
	cur_scored_primers = df.copy()[df.hairpin == False]
	scored_tm = []
	scored_gc = []
	scored_end = []
	scored_entropy = []
	scored_entropy_end = []
	for region in ["LoopF","LoopB"]:
		scored_tm+=list(cur_scored_primers[cur_scored_primers.region == region].tm.apply(lambda elem : score_with_ideal_min_max(thermo_params["TM"][region]["minimum"],thermo_params["TM"][region]["maximum"],thermo_params["TM"][region]["ideal"],elem,2.3333)).to_numpy()) #weight results in -0.7 score for 1 degree off from minimum or maximum
		scored_gc+=list(cur_scored_primers[cur_scored_primers.region == region].gc.apply(lambda elem : score_with_ideal_min_max(thermo_params["GC"][region]["minimum"],thermo_params["GC"][region]["maximum"],thermo_params["GC"][region]["ideal"],elem,23.333)).to_numpy()) #weight results in -0.7 score for 10% off from minimum or maximum
	scored_end+=list(cur_scored_primers.end.apply(lambda elem: 1 if elem <= thermo_params["END_STABILITY"]["maximum"] else 0).to_numpy())
	scored_entropy+=list(cur_scored_primers.entropy.apply(lambda elem : 1-(elem / cur_scored_primers.entropy.max())).to_numpy())
	#scored_entropy[region] = np.true_divide(scored_entropy[region],np.amax(scored_entropy[region]))
	scored_entropy_end+=list(cur_scored_primers.entropy_end.apply(lambda elem : 1-(elem / cur_scored_primers.entropy_end.max())).to_numpy())
	total_score = np.array(scored_tm) + np.array(scored_gc) + np.array(scored_end) #+ np.array(scored_entropy)*1/max(scored_entropy) + np.array(scored_entropy_end) * 1/max(scored_entropy_end)

	cur_scored_primers["score"] = total_score
	print(cur_scored_primers.sort_values(by="score",ascending=False).head())

	best_row_f = cur_scored_primers[cur_scored_primers.region == "LoopF"].sort_values(by="score",ascending=False).iloc[0]
	best_row_b = cur_scored_primers[cur_scored_primers.region == "LoopB"].sort_values(by="score",ascending=False).iloc[0]
	out_df = pd.DataFrame()
	for loop,row in (("LoopF",best_row_f),("LoopB",best_row_b)):
		seq = ""
		for offset in range(row["length"]):
			count = {}
			for l in get_partial_sequences(row["pos"]+offset,1).most_common():
				count[str(l[0])] = l[1]
			seq+=get_ambiguous(count.get("A",0),count.get("G",0),count.get("C",0),count.get("T",0))
		if loop == "LoopF":
			row["Seq"] = str(Seq(seq).reverse_complement())
		elif loop == "LoopB":
			row["Seq"] = str(Seq(seq))
		else:
			raise Exception("NO")
		out_df = out_df.append(row)
	out_df.head().to_csv(snakemake.output[0],index=False)
	




	

primers = pd.read_csv(snakemake.input[0])

open_fasta = lambda file: list(SeqIO.parse(open(file),"fasta"))

ref = open_fasta(snakemake.input[1])
quasispecies = open_fasta(snakemake.input[2])
nt_prev = pd.read_csv(snakemake.input[3])
params = json.load(open(snakemake.input[4]))
assay_params = json.load(open(snakemake.input[5]))
thermo_params = json.load(open(snakemake.input[6]))
ambiguous_to_nt = json.load(open(snakemake.input[7]))

calc_base_pair_shannon(nt_prev)

main()


