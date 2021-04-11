#!/usr/bin/env python3
import os
import argparse
import re

global MATCHING
MATCHING = "_\d+$"

global ROOT_CUR_NUMBERS
ROOT_CUR_NUMBERS = []



def parse():
	parser = argparse.ArgumentParser()
	parser.add_argument("-n", default=1, help="Number of iterations", type=int)
	parser.add_argument("--inputs",nargs="+",type=str)
	args = parser.parse_args()

	roots_temp = []
	for inp in args.inputs:
		match = re.search(MATCHING, inp)
		if match is None:
			raise Exception(f"{inp} needs to end with {MATCHING}")
		else:
			roots_temp.append(inp[:match.span()[0]])
			ROOT_CUR_NUMBERS.append(int(inp[match.span()[0]+1:]))
	args.inputs = tuple(roots_temp)
	return args

arguments = parse()
print(arguments)
def run(names):
	exit_code_rev = os.system(f"snakemake -j1 --use-conda ../Input/{names[0]}_{ROOT_CUR_NUMBERS[0]}_rev/quasispecies.fasta") >> 8
	if exit_code_rev == 1:
		raise Exception("Error with snakemake execution")
	exit_code_run = os.system(f"snakemake -j2 --use-conda ../Output/"+"{"+f"{names[0]}_{ROOT_CUR_NUMBERS[0]}"+f",{names[0]}_{ROOT_CUR_NUMBERS[0]}_rev"+"}"+"/untargetable.fasta")
	if exit_code_run == 1:
		raise Exception("Error with snakemake execution")

def next_iteration(names):
	value_norm = 0
	with open(f"../Output/{names[0]}_{ROOT_CUR_NUMBERS[0]}/primer_sets_coverage.txt") as h:
		value_norm = float(h.readline())
	value_rev = 0
	with open(f"../Output/{names[0]}_{ROOT_CUR_NUMBERS[0]}_rev/primer_sets_coverage.txt") as h:
		value_rev = float(h.readline())

	if value_norm >= value_rev:
		os.mkdir(f"../Input/{names[0]}_{ROOT_CUR_NUMBERS[0]+1}")
		with open(f"../Input/{names[0]}_{ROOT_CUR_NUMBERS[0]+1}/quasispecies.fasta",mode="w") as h_write,open(f"../Output/{names[0]}_{ROOT_CUR_NUMBERS[0]}/untargetable.fasta",mode="w") as h_read:
			for line in h_read:
				h_write.write(line)
	else:
		os.mkdir(f"../Input/{names[0]}_{ROOT_CUR_NUMBERS[0]+1}_rev")
		with open(f"../Input/{names[0]}_{ROOT_CUR_NUMBERS[0]+1}_rev/quasispecies.fasta",mode="w") as h_write,open(f"../Output/{names[0]}_{ROOT_CUR_NUMBERS[0]}_rev/untargetable.fasta",mode="w") as h_read:
			for line in h_read:
				h_write.write(line)

	for i in range(len(ROOT_CUR_NUMBERS)):
		ROOT_CUR_NUMBERS[i]+=1
	
				

	

run(arguments.inputs)
