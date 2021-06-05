#!/usr/bin/env python3
import os
import argparse
import re
from Bio import SeqIO

def parse() -> tuple:
	parser = argparse.ArgumentParser()
	parser.add_argument("-n", default=1, nargs="+",help="Number of iterations", type=int)
	parser.add_argument("--inputs",nargs="+",type=str)
	parser.add_argument("--continuation",type=bool,default=False)
	args = parser.parse_args()


	if args.inputs == None:
		raise Exception("Invalid inputs")

	files_len = len(args.inputs)
	generations = []
	for n in args.n:
		if n <= 0:
			raise Exception("Invalid inputs")
		else:
			generations.append(n)
	if len(generations) == 1:
		generations = generations * files_len
	elif len(generations) != files_len:
		raise Exception("Mismatch 'n' and 'inputs' sizes")

	files = [File.construct_file(file_name,max_gen) for file_name,max_gen in zip(args.inputs,generations)]

	return (files,args.continuation)

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
		return self.base_name + f"_{self.number}" if self.parity else self.base_name  + f"_{self.number}" + File.PARITY_IDENTIFIER

	def get_double_negative_full_name(self):
		if not self.parity:
			raise Exception("parity needs to be True for this method!")
		return self.base_name + f"_{self.number}" + File.PARITY_IDENTIFIER+File.PARITY_IDENTIFIER

	def get_expanded_full_name(self):
		return self.get_full_name() if not self.parity else self.get_double_negative_full_name()

	def file_increment(self,change:int):
		return File(self.parity,self.base_name,self.number+change,self.generation,self.maximum_generation)
	
	@staticmethod
	def construct_file(full_name:str,maximum_generation:int):
		match_object_fr = re.search("_\d+$",full_name)
		match_object_rev = re.search("_\d+"+File.PARITY_IDENTIFIER+"$",full_name)
		if match_object_fr is not None:
			result = match_object_fr.span()
			return File(True, full_name[:result[0]],int(full_name[result[0]+1:result[1]]),1,maximum_generation)
		elif match_object_rev:
			result = match_object_rev.span()
			return File(False,full_name[:result[0]],int(full_name[result[0]+1:result[1] - len(File.PARITY_IDENTIFIER)]),1,maximum_generation)
		else:
			raise Exception("Incorrect file name")

	def create_opposite(self):
		return File(not self.parity,self.base_name,self.number,self.generation,self._maximum_generation)

	@staticmethod
	def snakemake_construct_opposites(files,maximum_cores=4):
		cores_used = min(maximum_cores,len(files))
		new_files = [file.create_opposite() for file in files]
		brace_begin = "{" if len(files) > 1 else ""
		brace_end = "}" if len(files) > 1 else ""
		command = "snakemake -j" + str(cores_used) + " --use-conda ../Input/"+brace_begin + ",".join([file.get_expanded_full_name() for file in new_files]) + brace_end+"/quasispecies.fasta"
		assert os.system(command) == 0
		for file in new_files:
			if file.parity:
				assert os.system(f"mv ../Input/{file.get_double_negative_full_name()} ../Input/{file.get_full_name()}") == 0
		return new_files

	@staticmethod
	def snakemake_run(files_total,maximum_cores=4):
		cores_used = min(maximum_cores,len(files_total))
		brace_begin = "{" if len(files_total) > 1 else ""
		brace_end = "}" if len(files_total) > 1 else ""
		command = "snakemake -j" + str(cores_used) + " --use-conda ../Output/"+brace_begin + ",".join([file.get_full_name() for file in files_total]) + brace_end+"/untargetable.fasta"
		print(command)
		assert os.system(command) == 0

	@staticmethod
	def create_next_files(files):
		new_files = []
		for file in files:
			new_file = File(file.parity,file.base_name,file.number+1,file.generation+1,file._maximum_generation)
			command = "mkdir ../Input/" + new_file.get_full_name() + "; cp ../Output/"+file.get_full_name()+"/untargetable.fasta ../Input/"+new_file.get_full_name()+"/quasispecies.fasta; cp ../Input/" + file.get_full_name() + "/reference.fasta ../Input/"+new_file.get_full_name()+"/reference.fasta"
			assert os.system(command) == 0
			new_files.append(new_file)
		return new_files

	@staticmethod
	def evaluate_and_create_next_files_if_possible(files,files_opposite):
		new_files = []
		for file_pair in zip(files,files_opposite):
			lengths = [len(list(SeqIO.parse(open(f"../Output/{file.get_full_name()}/untargetable.fasta"),"fasta"))) for file in file_pair]
			if lengths[0] > lengths[1]:
				chosen = file_pair[1]
				unchosen = file_pair[0]
			else:
				if lengths[0] == lengths[1] and lengths[1] == 0:
					continue
				chosen = file_pair[0]
				unchosen = file_pair[1]

			put_away_unnecessary = "mv ../Output/"+unchosen.get_full_name()+" ../Output/"+unchosen.get_full_name()+"_UNSELECTED"
			if file_pair[0]._maximum_generation == file_pair[0].generation:
				assert os.system(put_away_unnecessary) == 0
				continue
			new_file = File(chosen.parity,chosen.base_name,chosen.number+1,chosen.generation+1,chosen._maximum_generation)
			command = "mkdir ../Input/" + new_file.get_full_name() + "; cp ../Output/"+chosen.get_full_name()+"/untargetable.fasta ../Input/"+new_file.get_full_name()+"/quasispecies.fasta; cp ../Input/" + chosen.get_full_name() + "/reference.fasta ../Input/"+new_file.get_full_name()+"/reference.fasta"
			assert os.system(command) == 0
			assert os.system(put_away_unnecessary) == 0
			new_files.append(new_file)
		return new_files



	def __str__(self):
		return str(self.__dict__)



files,cont = parse()

if cont:
	files = File.create_next_files(files)

while len(files) >= 1:
	for f in files:
		print(f)
	opposites = File.snakemake_construct_opposites(files)
	File.snakemake_run(files+opposites)
	files = File.evaluate_and_create_next_files_if_possible(files,opposites)


