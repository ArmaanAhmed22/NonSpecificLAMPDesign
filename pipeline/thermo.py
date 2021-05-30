from Bio import SeqIO
import json
from tqdm import tqdm
import pandas as pd
import numpy as np
import primer3
from collections import Counter
import math

def calc_median_end_stability(region,partial_sequences,primer3_params):
	#Maybe adjust?
	end_stabilities = []
	for partial_seq,number in partial_sequences.items():
		if region == "F3" or region == "F2" or region == "B1c" or region == "LoopB": #✓
			target = partial_seq.reverse_complement() #Already 3'
			primer = partial_seq
			if region == "B1c": #Change to 5'
				target = target[::-1]
				primer = primer[::-1]
		elif region == "F1c" or region == "B2" or region == "B3" or region == "LoopF": #✓
			target = partial_seq
			primer = partial_seq.reverse_complement() #3' already
			if region == "F1c": #Change to 5'
				target = target[::-1]
				primer = primer[::-1]
		else:
			raise Exception(f"{region} region does not exist")

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
		if region == "F3" or region == "F2" or region == "B1c" or region == "LoopB": #✓
			primer = str(partial_seq)
		elif region == "F1c" or region == "B2" or region == "B3" or region == "LoopF": #✓
			primer = str(partial_seq.reverse_complement()) #Keep everything 5' -> 3'
		else:
			raise Exception(f"{region} region does not exist")
		tm = primer3.bindings.calcTm(primer,**primer3_params)
		tms.extend([tm for _ in range(number)])

	tms = np.array(tms)
	median = np.median(tms)
	return median

def calc_median_gc(region,partial_sequences):
	gcs = []
	for partial_seq,number in partial_sequences.items():
		if region == "F3" or region == "F2" or region == "B1c" or region == "LoopB": #✓
			primer = str(partial_seq)
		elif region == "F1c" or region == "B2" or region == "B3" or region == "LoopF": #✓
			primer = str(partial_seq.reverse_complement())
		else:
			raise Exception(f"{region} region does not exist")
		gc_amt = 0
		for nuc in primer:
			if nuc == "G" or nuc == "C":
				gc_amt+=1
		gcs.extend([gc_amt / len(primer) for _ in range(number)])
	gcs = np.array(gcs)
	median = np.median(gcs)
	return median

def calc_hairpin(region,partial_sequences,primer3_params,threshold):
	hairpins = []
	for partial_seq,number in partial_sequences.items():
		if region == "F3" or region == "F2" or region == "B1c" or region == "LoopB":
			primer = str(partial_seq)
		elif region == "F1c" or region == "B2" or region == "B3" or region == "LoopF":
			primer = str(partial_seq.reverse_complement())
		else:
			raise Exception(f"{region} region does not exist")
		struc_found = 1 if primer3.bindings.calcHairpin(primer,**primer3_params).structure_found else 0
		hairpins.extend([struc_found for _ in range(number)])
	hairpins = np.array(hairpins)
	mean = np.mean(hairpins)
	return mean >= threshold
