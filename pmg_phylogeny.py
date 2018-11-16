#!/home/frankaylward/anaconda_ete/bin/python
import argparse, os, sys, subprocess, re, shlex, pandas, operator
from collections import defaultdict
from Bio import SeqIO

# specify directory where the script is located
homedir = "/home/frankaylward/Desktop/marker_gene_benchmarking/pangenomics/Multilocus_Phylogenetics"

# This script requires that HMMER3 and Prodigal are installed in the PATH. 

# create and wipe the LOGFILE before beginning
o = open("LOGFILE.txt", "w")
o.close()

#########################################
### define prodigal launcher function ###
#########################################
suffix = ".fa"
def predict_genes(folder):
	for filenames in os.listdir(folder):
		if filenames.endswith(".fna") or filenames.endswith(".fa"):
			s = filenames.split(".")
			suffix = "."+s[-1]

			#print "Predicting genes..."
			input_file = os.path.join(folder, filenames)
			protein_file = re.sub(suffix, ".faa", input_file)
			gff_file = re.sub(suffix, ".gff", input_file)
			acc = re.sub(suffix, "", filenames)

			cmd = "prodigal -i "+ input_file +" -f gff -o "+ gff_file +" -a "+ protein_file 
			#print cmd
			cmd2 = shlex.split(cmd)
			subprocess.call(cmd2, stdout=open("LOGFILE.txt", 'w'), stderr=open("LOGFILE.txt", 'a'))


def make_nr(folder):
	for filenames in os.listdir(folder):
		if filenames.endswith(".faa"):

			protein_file = os.path.join(folder, filenames)
			acc1 = re.sub(".faa", "", filenames)
			acc = re.sub("_protein", "", acc1)

			# make protein names non-redundant
			output_seqs = []
			handle = open(protein_file, "r")
			for record in SeqIO.parse(handle, "fasta"):
				name = record.id
				if ".." in name:
					pass
				else:
					record.id = acc +".."+ name
					#print record.id, acc, filenames, name

				output_seqs.append(record)

			handle.close()
			output = open(protein_file, "w")
			SeqIO.write(output_seqs, output, "fasta")

# end

#######################################
### define HMMER3 launcher function ###
#######################################
def run_hmmer(folder, db_path, cutoff, cpus):
	for filenames in os.listdir(folder):
		if filenames.endswith(".faa"):
			input_file = os.path.join(folder, filenames)
			hmm_out = re.sub(".faa", ".hmmout", input_file)

			cmd = cmd = "hmmsearch --cpu "+ cpus +" --tblout " + hmm_out + " "+ cutoff + " " + db_path + " " + input_file
			cmd2 = shlex.split(cmd)
			subprocess.call(cmd2, stdout=open("LOGFILE.txt", 'w'), stderr=open("LOGFILE.txt", 'a'))
# end

######################################
######## get protein dictionary ######
######################################
def parse_faa(folder, dictionary):
	seq_dict = {}
	for filenames in os.listdir(folder):
		if filenames.endswith(".faa"):
			input_file = os.path.join(folder, filenames)
			for j in SeqIO.parse(input_file, "fasta"):
				#print j.id
				if j.id in dictionary:
					#print j.id
					seq_dict[j.id] = j
	return seq_dict
# end

#############################################
######## get protein length dictionary ######
#############################################
def get_length(input_file):
	seq_dict = {}
	for j in SeqIO.parse(input_file, "fasta"):
		seq_dict[j.id] = len(j.seq)
	return seq_dict
# end


####################################
#### define HMM output parser ######
####################################
def hmm_parser(folder, output):
	full_hit_dict = {}
	cog_dict = defaultdict(list)
	score_list = {}
	prot_list = []
	combined_output = open(output, "w")
	#combined_output.write("protein\tacc\thit\tscore\tlength\tcategory\tspecies\tdomain\tphylum\torder\tp2\tp3\n")
	hits = []
	bit_dict = {}
	df = pandas.DataFrame()

	for filenames in os.listdir(folder):
		if filenames.endswith(".hmmout"):

			protein_file = re.sub(".hmmout", ".faa", filenames)
			length_dict = get_length(os.path.join(folder, protein_file))

			acc = re.sub(".hmmout", "", filenames)
			f = open(os.path.join(folder, filenames), 'r')
			hit_dict = {}
			bit_dict = defaultdict(int)
			hit_type = {}
			marker_dict = {}

			for line in f.readlines():
				if line.startswith("#"):
					pass
				else:
					newline = re.sub( '\s+', '\t', line)
					list1 = newline.split('\t')
					ids = list1[0]
					hit = list1[2]
					acc2 = list1[3]
					if "COG" in hit:
						pass
					else:
						hit = acc2

					bit_score = list1[5]
					score = float(bit_score)

					if score > bit_dict[ids]:
						hit_dict[ids] = hit
						bit_dict[ids] = score

			bit_sorted = sorted(bit_dict.items(), key=operator.itemgetter(1), reverse=True)
			output_list = []
			for item in bit_sorted:
				entry = item[0]
				if entry in hit_dict:
					#if entry in bit_dict and entry in hit_dict:
					output_list.append(entry +"\t"+ str(hit_dict[entry]) +"\t"+ str(bit_dict[entry]))
					full_hit_dict[entry] = entry
					#print entry
				#else:
				#	print entry

			#parsed = open(folder+"/"+filenames+".parsed", 'r')
			hit_profile = defaultdict(int)
			done = []
			for line in output_list:
				line1 = line.rstrip()
				tabs = line1.split("\t")
				ids = tabs[0]
				hits.append(ids)
				cog = tabs[1]
				score = float(tabs[2])
				nr = acc +"_"+ cog

				if nr in done:
					#combined_output.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ score +"\tNH\t"+ "\n")
					pass
				else:
					#if score >= 0:
					#print score
					combined_output.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ str(score) +"\t"+ str(length_dict[ids]) +"\tBH\t"+ "\n")
					done.append(nr)

			#s1 = pandas.DataFrame(pandas.Series(hit_profile, name = acc))
			#df = pandas.concat([df, s1], axis=1)
	return full_hit_dict
# end


def run_program(input_dir, cogs_file, project, faa_out, db, parse, nonredundant, genes, cpus, verbose):

	if os.path.isdir(project):
		pass
	else:
		os.mkdir(project)

	project_name = project+".all_hmm_out.txt"
	######## Check if genes have already been predicted.
	if not (genes):
		if verbose:
			print "Predicting genes...\n"

		predict_genes(input_dir)
	else:
		if verbose:
			print "Genes already predicted.\n"

	####### Check if predicted proteins should be made non-redundant
	if (nonredundant):
		if verbose:
			print "Making protein entries non-redundant...\n"

		make_nr(input_dir)
	else:
		if verbose:
			print "Genes already non-redundant.\n"

	###### Get a list of the HMM entries that should be used
	cogs = []
	list_file = open(os.path.join(homedir, "hmm/", db+".list"), "r")
	for i in list_file.readlines():
		line = i.rstrip()
		cogs.append(line)

	###########################################
	################ run functions ############
	###########################################
	cutoff_dict = {}
	if "checkm" in db:
		cutoff = "--cut_nc"
	else:
		cutoff = " "
		cutoff_file = open(os.path.join(homedir, "hmm/embl_cutoffs.txt"), "r")
		for i in cutoff_file.readlines():
			line = i.rstrip()
			tabs = line.split("\t")
			name = tabs[0]
			value = tabs[1]
			cutoff_dict[name] = float(value)

	db_path = os.path.join("hmm", db+".hmm")
	if not (parse):
		if verbose:
			print "Running HMMER3...\n"

		run_hmmer(input_dir, db_path, cutoff, cpus)

		if verbose:
			print "Parsing HMMER3 outputs...\n"

		all_hits = hmm_parser(input_dir, os.path.join(project, project_name))
	else:
		if verbose:
			print "Only parsing all.hmmout.txt file\n"

	###########################################################################
	######### Final compilation of files for input to ete3 ####################
	###########################################################################
	if verbose:
		print "Compilation final files...\n"

	hmm_out = open(os.path.join(project, project_name), "r")
	acc_dict = {}
	cog_dict = defaultdict(list)
	all_hits_dict = {}
	for j in hmm_out.readlines():
		line = j.rstrip()
		tabs = line.split("\t")
		hit_type = tabs[5]
		if hit_type == "BH":
			protein = tabs[0]
			all_hits_dict[protein] = protein

			acc = tabs[1]
			cog = tabs[2]
			score = tabs[3]
			if db == "embl":
				if score >= cutoff_dict[cog]:
					acc_dict[protein] = acc
					cog_dict[cog].append(protein)
			elif db == "RNAP":
				if score >= 100:
					acc_dict[protein] = acc
					cog_dict[cog].append(protein)			
			else:
				acc_dict[protein] = acc
				cog_dict[cog].append(protein)

	# get a dictionary that links protein IDs to their SeqIO records
	if verbose:
		print "Getting protein dictionary...\n"

	seq_dict = parse_faa(input_dir, all_hits_dict)

	cog_file = hmm_out = open(os.path.join(project, cogs_file), "w")
	out_proteins = open(os.path.join(project, faa_out), "w")
	output_seqs = []
	tally = []
	tally2 = []
	already_done = {}
	for n in cogs:
		protein_list = cog_dict[n]
		#print protein_list
		for i in protein_list:
			acc = acc_dict[i]
			acc2 = re.sub("_", ".", acc)
			new_name = acc2 +"_"+ n
		
			if new_name in already_done:
				print n, i, acc, acc2

			else:
				already_done[new_name] = new_name

				cog_file.write(new_name +"\t")
				tally.append(new_name)

				record = seq_dict[i]
				record.description = record.id
				#record.description = ""
				record.id = new_name
				output_seqs.append(record)
				tally2.append(record.id)

		cog_file.write("\n")

	SeqIO.write(output_seqs, out_proteins, "fasta")
	#print len(tally), len(tally2), len(set(tally)), len(set(tally2))

########################################################################
##### use argparse to run through the command line options given #######
########################################################################
def main(argv=None):
	args_parser = argparse.ArgumentParser(description="Program for identifying phylogenetic marker genes (PMGs) in sequencing data\nFrank O. Aylward, Assistant Professor, Virginia Tech Department of Biological Sciences <faylward at vt dot edu>", epilog="Virginia Tech Department of Biological Sciences")
	args_parser.add_argument('-i', '--input_folder', required=True, help='Input folder of FNA sequence files')
	args_parser.add_argument('-c', '--cogs', required=True, help='output file for the COGs file, used in ete3')
	args_parser.add_argument('-p', '--project', required=False, default="projectx", help='provide a project name for this analysis')
	args_parser.add_argument('-o', '--output_faa', required=True, help='amino acid fasta output file of the PMGs identified')
	args_parser.add_argument('-t', '--cpus', required=False, default=str(1), help='Number of CPUs to use for HMMER3 search (default=1)')
	args_parser.add_argument('-v', '--verbose', required=False, default=True, help='Verbosity')
	args_parser.add_argument('-f', '--genes', type=bool, default=False, const=True, nargs='?', help='Use this genes have already been predicted, and .faa files containing predicted proteins are present in the input folder)')
	args_parser.add_argument('-pa', '--parse_only', type=bool, default=False, const=True, nargs='?', help='Implies that the all_hmm.out file is already present and should be parsed. Skips HMM search, gene prediction, etc)')
	args_parser.add_argument('-n', '--nr', type=bool, default=False, const=True, nargs='?', help='Use this if you would like pmg_phylogeny to go through and re-name all of the proteins to ensure they are non-redundant')
	args_parser.add_argument('-d', '--db', required=True, help='Database for matching. Accepted sets are "embl", "checkm_bact", "checkm_arch", "rpl1_arch", "rpl1_bact", or "rpl2"')
	args_parser = args_parser.parse_args()

	# set up object names for input/output/database folders
	input_dir = args_parser.input_folder
	cogs_file = args_parser.cogs
	project = args_parser.project
	faa_out = args_parser.output_faa
	db = args_parser.db
	parse = args_parser.parse_only
	nonredundant = args_parser.nr
	genes = args_parser.genes
	cpus = args_parser.cpus
	verbose = args_parser.verbose

	run_program(input_dir, cogs_file, project, faa_out, db, parse, nonredundant, genes, cpus, verbose)

	return 0

if __name__ == '__main__':
	status = main()
	sys.exit(status)

# end

	# and then run through ete3 using the command "ete3 build -a proteins_for_phylogeny.faa --cogs cog_file.txt -w standard_trimmed_fasttree -m sptree_fasttree_90 --clearall -C 8 -o test_tree"






