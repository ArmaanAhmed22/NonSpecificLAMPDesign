from Bio import SeqIO
import json
from tqdm import tqdm
import pandas as pd
import numpy as np
import primer3
from collections import Counter
import math


def calc_base_pair_shannon(nuc_prevalence):
	nuc_prevalence["total"] = nuc_prevalence["A"] + nuc_prevalence["G"] + nuc_prevalence["C"] + nuc_prevalence["T"]
	shannon = None
	for nuc in ["A","G","C","T"]:
		percent = nuc_prevalence[nuc] / nuc_prevalence["total"]
		if shannon is None:
			shannon = np.where(percent == 0,percent,-percent * np.log2(percent))
		else:
			shannon += np.where(percent == 0,percent,-percent * np.log2(percent))
	nuc_prevalence["shannon"] = shannon

def calc_median_end_stability(region,partial_sequences,primer3_params):
	#Maybe adjust?
	end_stabilities = []
	for partial_seq,number in partial_sequences.items():
		if region == "F3" or region == "F2" or region == "B1c": #✓
			target = partial_seq.reverse_complement() #Already 3'
			primer = partial_seq
			if region == "B1c": #Change to 5'
				target = target[::-1]
				primer = primer[::-1]
		elif region == "F1c" or region == "B2" or region == "B3": #✓
			target = partial_seq
			primer = partial_seq.reverse_complement() #F1c, 5' already
			if region != "F1c": #Change to 3'
				target = target[::-1]
				primer = primer[::-1]

		target = str(target)
		primer = str(primer)
		end_stability = primer3.bindings.calcEndStability(primer,target,**primer3_params).dg / 1000
		end_stabilities.extend([end_stability for _ in range(number)])
	
	end_stabilities = np.array(end_stabilities)
	median = np.median(end_stabilities)
	return median

def calc_median_tm(region,partial_sequences,primer3_params):
	tms = []
	for partial_seq,number in partial_sequences.items():
		if region == "F3" or region == "F2" or region == "B1c": #✓
			primer = str(partial_seq)
		elif region == "F1c" or region == "B2" or region == "B3": #✓
			primer = str(partial_seq.reverse_complement()[::-1]) #Keep everything 5' -> 3'
		tm = primer3.bindings.calcTm(primer,**primer3_params)
		tms.extend([tm for _ in range(number)])

	tms = np.array(tms)
	median = np.median(tms)
	return median

def calc_median_gc(region,partial_sequences):
	gcs = []
	for partial_seq,number in partial_sequences.items():
		if region == "F3" or region == "F2" or region == "B1c": #✓
			primer = str(partial_seq)
		elif region == "F1c" or region == "B2" or region == "B3": #✓
			primer = str(partial_seq.reverse_complement()[::-1])
		gc_amt = 0
		for nuc in primer:
			if nuc == "G" or nuc == "C":
				gc_amt+=1
		gcs.extend([gc_amt / len(primer) for _ in range(number)])
	gcs = np.array(gcs)
	median = np.median(gcs)
	return median

def calc_hairpin(region,partial_sequences,primer3_params):
	hairpins = []
	for partial_seq,number in partial_sequences.items():
		if region == "F3" or region == "F2" or region == "B1c":
			primer = str(partial_seq)
		elif region == "F1c" or region == "B2" or region == "B3":
			primer = str(partial_seq.reverse_complement()[::-1])
		struc_found = 1 if primer3.bindings.calcHairpin(primer,**primer3_params).structure_found else 0
		hairpins.extend([struc_found for _ in range(number)])
	hairpins = np.array(hairpins)
	median = np.median(hairpins)
	return math.floor(median) == 1



def get_primer_locs(start,primer_lengths,primer_spacing):
	primer_locs = {}

	primer_locs["F3"] = start
	primer_locs["F2"] = primer_locs["F3"] + primer_lengths["F3"] + primer_spacing["F3|3',F2|5'"]
	primer_locs["F1c"] = primer_locs["F2"] + primer_lengths["F2"] + (primer_spacing["F2|5',F1c|3'"] - primer_lengths["F2"])

	primer_locs["B2"] = primer_locs["F2"] + primer_lengths["F2"] + (primer_spacing["F2|5',B2|5'"] - primer_lengths["B2"] - primer_lengths["F2"])
	primer_locs["B1c"] = primer_locs["B2"] + primer_lengths["B2"] - (primer_spacing["B1c|3',B2|5'"] - primer_lengths["B1c"])
	primer_locs["B3"] = primer_locs["B2"]+primer_lengths["B2"]+primer_spacing["B2|5',B3|3'"]

	return primer_locs

def get_partial_sequences(region_pos,region_len):
	holding = []
	for record in quasispecies:
		region = record.seq[region_pos:region_pos+region_len]
		holding.append(region)
	return Counter(holding)



def individual_primer_generate(reference,primer_lengths,primer_spacing,nuc_prevalence):
	#df = pd.DataFrame(columns=["region","pos","length","entropy"])
	def generate_region(region,downstream_requirements,upstream_requirements,depend_on_length=True):
		primer3_params_end = {
			"dv_conc":assay_params["DV_CONC"],
			"dntp_conc":assay_params["DNTP_CONC"],
			"dna_conc":assay_params[region],
			"temp_c":assay_params["TEMP_C"]
		}
		primer3_params_tm = {
			"dv_conc":assay_params["DV_CONC"],
			"dntp_conc":assay_params["DNTP_CONC"],
			"dna_conc":assay_params[region],
		}
		data = []
		for region_len in tqdm(range(primer_lengths[region]["minimum"],primer_lengths[region]["maximum"]+1),desc=region):
			primer_end = params["length"][region]["minimum"] // 3
			primer_rest = region_len - primer_end
			region_pos = downstream_requirements
			while (region_pos + ((region_len if depend_on_length else 0)+ upstream_requirements) <= len(reference)):
				partial_sequences = get_partial_sequences(region_pos,region_len)
				entropy = nuc_prevalence["shannon"][region_pos:region_pos+region_len].sum()

				if region == "F3" or region == "F2":
					entropy = nuc_prevalence["shannon"][region_pos:region_pos+primer_rest].sum()
					entropy_end = nuc_prevalence["shannon"][region_pos+primer_rest:region_pos+region_len].sum()
				elif region == "B1c":
					entropy_end = nuc_prevalence["shannon"][region_pos:region_pos+primer_end].sum()
					entropy = nuc_prevalence["shannon"][region_pos+primer_end:region_pos+region_len].sum()
				elif region == "F1c":
					entropy = nuc_prevalence["shannon"][region_pos:region_pos+primer_rest].sum()
					entropy_end = nuc_prevalence["shannon"][region_pos+primer_rest:region_pos+region_len].sum() 
				elif region == "B2" or region == "B3":
					entropy_end = nuc_prevalence["shannon"][region_pos:region_pos+primer_end].sum()
					entropy = nuc_prevalence["shannon"][region_pos+primer_end:region_pos+region_len].sum()
				#Calc Tm
				tm = calc_median_tm(region,partial_sequences,primer3_params_tm)
				#Calc End stability
				end = calc_median_end_stability(region,partial_sequences,primer3_params_end)
				#Calc GC
				gc = calc_median_gc(region,partial_sequences)
				#Calc hairpin
				hairpin = calc_hairpin(region,partial_sequences,primer3_params_end)
				data.append({"region":region,"pos":region_pos,"length":region_len,"entropy":entropy,"entropy_end":entropy_end,"tm":tm,"end":end,"gc":gc,"hairpin":hairpin})
				region_pos+=1
		return data

	data_hold = []

	#Generate F3 primer data
	data_hold.extend(
	generate_region("F3",
		0,
		primer_spacing["F3|3',F2|5'"]["minimum"] + primer_spacing["F2|5',B2|5'"]["minimum"] + primer_spacing["B2|5',B3|3'"]["minimum"] + primer_lengths["B3"]["minimum"])
	)
	
	#Generate F2 primer data
	data_hold.extend(
	generate_region("F2",
		primer_lengths["F3"]["minimum"] + primer_spacing["F3|3',F2|5'"]["minimum"],
		primer_spacing["F2|5',B2|5'"]["minimum"] + primer_spacing["B2|5',B3|3'"]["minimum"] + primer_lengths["B3"]["minimum"],
		depend_on_length=False)
	)
	
	#Generate F1c primer data
	data_hold.extend(
	generate_region("F1c",
		primer_lengths["F3"]["minimum"]+primer_spacing["F3|3',F2|5'"]["minimum"]+primer_spacing["F2|5',F1c|3'"]["minimum"],
		-primer_spacing["F2|5',F1c|3'"]["maximum"] + primer_spacing["F2|5',B2|5'"]["minimum"] + primer_spacing["B2|5',B3|3'"]["minimum"] + primer_lengths["B3"]["minimum"],
		depend_on_length=False
		)
	)

	#Generate B1c primer data
	data_hold.extend(
	generate_region("B1c",
		primer_lengths["F3"]["minimum"] + primer_spacing["F3|3',F2|5'"]["minimum"] + primer_spacing["F2|5',B2|5'"]["minimum"] - primer_spacing["B1c|3',B2|5'"]["maximum"] - primer_lengths["B1c"]["maximum"],
		primer_spacing["B1c|3',B2|5'"]["minimum"]+primer_spacing["B2|5',B3|3'"]["minimum"] + primer_lengths["B3"]["minimum"]
		)
	)

	#Generate B2 primer data
	data_hold.extend(
	generate_region("B2",
		primer_lengths["F3"]["minimum"] + primer_spacing["F3|3',F2|5'"]["minimum"] + primer_spacing["F2|5',B2|5'"]["minimum"] - primer_lengths["B2"]["maximum"],
		primer_spacing["B2|5',B3|3'"]["minimum"] + primer_lengths["B3"]["minimum"]
		)
	)
	# Generate B3 primer data
	data_hold.extend(
	generate_region("B3",
		primer_lengths["F3"]["minimum"] + primer_spacing["F3|3',F2|5'"]["minimum"] + primer_spacing["F2|5',B2|5'"]["minimum"] + primer_spacing["B2|5',B3|3'"]["minimum"],
		0
		)
	)
	
	result = pd.DataFrame(data_hold)
	return result


#open reference
reference = next(SeqIO.parse(open(snakemake.input[0],"r"),"fasta")).seq

#Quasispecies sequences
quasispecies = list(SeqIO.parse(open(snakemake.input[1],"r"),"fasta"))

#open nuc prevalence to calculate shannon entropy
nuc_prevalence = pd.read_csv(snakemake.input[2])

#Data from https://primerexplorer.jp/e/v5_manual/
params = json.load(open(snakemake.input[3],"r"))

assay_params = json.load(open(snakemake.input[4],"r"))

thermo_params = json.load(open(snakemake.input[5],"r"))

calc_base_pair_shannon(nuc_prevalence)

df = individual_primer_generate(reference,params["length"],params["spacing"],nuc_prevalence)
df.to_csv(snakemake.output[0],index=False)