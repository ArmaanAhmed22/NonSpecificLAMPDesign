#!/usr/bin/env python3
import os
import argparse
import re
from Bio import SeqIO

def parse() -> tuple:
	parser = argparse.ArgumentParser()
	parser.add_argument("-n", default=1, nargs="+",help="Number of iterations", type=int)
	parser.add_argument("--inputs",nargs="+",type=str)
	args = parser.parse_args()

	if args.inputs == None:
		raise Exception("Invalid inputs")

	files_len = len(args.input)
	generations = []
	for n in args.n:
		if n <= 0:
			raise Exception("Invalid inputs")
		else:
			generations.append(n)
	if len(generations) == 1:
		generations = generations * len(files)
	elif len(generations) != files_len:
		raise Exception("Mismatch 'n' and 'inputs' sizes")

	files = [File.construct_file(file_name,max_gen) for file_name,max_gen in zip(args.inputs,generations)]

	return files

class File:
	PARITY_IDENTIFIER = "_rev"
	def __init__(self,parity: bool,base_name: str,number: int, generation: int, maximum_generation: int):
		"""
		"True" means forward, "False" means reverse complement
		"""
		self.parity = parity
		self.base_name = base_name
		self.number = number
		self.generation = generation
		self._maximum_generation = maximum_generation

	def get_full_name(self):
		return self.base_name + f"_{self.number}" if self.parity else self.base_name + f"_rev" + f"_{self.number}"

	def get_double_negative_full_name(self):
		if not self.parity:
			raise Exception("parity needs to be positive for this method!")
		return self.base_name + "_rev_rev" + f"_{self.number}"

	def get_expanded_full_name(self):
		return get_full_name() if self.parity else get_double_negative_full_name()

	@staticmethod
	def construct_file(full_name:str,maximum_generation:int):
		match_object = re.search("_\d+$",full_name)
		if match_object is not None:
			result = match_object.span()
			is_rev = len(full_name) - (result[1] - result[0]) >= len(File.PARITY_IDENTIFIER) and full_name[result[0] - len(File.PARITY_IDENTIFIER):result[0]] == File.PARITY_IDENTIFIER
			return File(not is_rev, full_name[:result[0] - len(File.PARITY_IDENTIFIER)] if is_rev else full_name[:result[0]],int(full_name[result[0]+1:result[1]]),1,maximum_generation)
		else:
			raise Exception("Incorrect file name")

	def create_opposite(self):
		return File(not self.parity,self.base_name,self.number,self.generation,self._maximum_generation)

	@staticmethod
	def snakemake_construct_opposites(files,maximum_cores=4):
		cores_used = min(maximum_cores,len(files))
		new_files = [file.create_opposite() for file in files]
		command = "snakemake -j" + str(cores_used) + " --use-conda ../Input/{" + ",".join([file.get_expanded_full_name() for file in new_files]) + "}/quasispecies.fasta"
		os.system(command)
		for file in new_files:
			if file.parity:
				os.system(f"mv ../Input/{file.get_double_negative_full_name()} ../Input/{file.get_full_name()}")
		return new_files

	@staticmethod
	def snakemake_run(files_total,maximum_cores=4):
		cores_used = min(maximum_cores,len(files))
		command = "snakemake -j" + str(cores_used) + " --use-conda ../Output/{" + ",".join([file.get_full_name() for file in files_total]) + "}/ready_primer_sets.json"
		os.system(command)

	@staticmethod
	def evaluate_and_create_next_files_if_possible(files,files_opposite):
		new_files = []
		for file_pair in zip(files,files_opposite):
			lengths = [len(list(SeqIO.parse(open(f"../Output/{file.get_full_name()}/untargetable.fasta"),"fasta"))) for file in file_pair]
			if lengths[0] > lengths[1]:
				chosen = file_pair[1]
				unchosen = file_pair[0]
			else:
				chosen = file_pair[0]
				unchosen = file_pair[1]

			put_away_unnecessary = "mv ../Output/"+unchosen.get_full_name()+" ../Output/"+unchosen.get_full_name()+"_UNSELECTED"
			if file_pair[0]._maximum_generation == file_pair[0].generation:
				os.system(put_away_unnecessary)
				continue
			new_file = File(chosen.parity,chosen.base_name,chosen.number+1,chosen.generation+1,chosen._maximum_generation)
			command = "mkdir ../Input/" + new_file.get_full_name() + "; cp ../Output/"+chosen.get_full_name()+"/untargetable.fasta ../Input/"+new_file.get_full_name()+"/quasispecies.fasta; cp ../Input/" + chosen.get_full_name() + "/reference.fasta ../Input/"+new_file.get_full_name()+"/reference.fasta"
			os.system(command)
			os.system(put_away_unnecessary)
			new_files.append(new_file)
		return new_files



	def __str__(self):
		return str(self.__dict__)



files = parse()


while len(files) >= 1:
	opposites = File.snakemake_construct_opposites(files)
	File.snakemake_run(files+opposites)
	files = File.evaluate_and_create_next_files_if_possible(files,opposites)


