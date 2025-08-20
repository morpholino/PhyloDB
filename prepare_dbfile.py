#!/usr/bin/python
import argparse
import hashlib
import os,shutil
import pandas as pd
from Bio import SeqIO
from ete3 import NCBITaxa

# Initialize NCBITaxa
ncbi = NCBITaxa()

def parse_arguments():
	# Create the argument parser
	parser = argparse.ArgumentParser(description="Reformat FASTA headers and log new taxonomies")
	# Add the filename argument
	parser.add_argument("-f", "--filename", type=str, help="Input FASTA file")
	# Add the species argument with nargs='+' to accept one or more arguments
	parser.add_argument("-s", "--species", nargs='+', type=str, help="Species name")
	# Parse the command-line arguments
	args = parser.parse_args()
	
	return args


def _def_variables():
	# Where the data is
	global datadir
	datadir = "/tscc/projects/ps-allenlab/projdata/common/db/PhyloDB_2.0/"
	
	# List of taxonomic ranks to keep
	# Do not keep "subphylum" for viruses
	global valid_ranks
	valid_ranks = [
		'superkingdom', 'kingdom', 'division', 'clade', 'phylum', 'subphylum',
		'class', 'order', 'family', 'genus', 'species'
	]
	
	global debug
	debug = False
	if debug:
		print("debug/verbose mode on!")


def _clean_name(string):
	badchars = "'[]/#()*+,-:="
	string = string.replace(" ", "_")
	return "".join([x for x in string if x not in badchars])


def _calculate_compatibility_score(query_annotations, pr2_annotations):
	query_ranks = set(query_annotations)
	pr2_ranks = set(pr2_annotations.values())
	common_ranks = query_ranks.intersection(pr2_ranks)
	score = len(common_ranks)
	return score, common_ranks


def _find_best_match(query_annotations, query_genus, pr2_database):
	best_match = None
	best_score = 0
	if debug:
		print("queried taxon:", query_genus)
	for index, row in pr2_database.iterrows():
		pr2_annotations = row.to_dict()
		score, common_ranks = _calculate_compatibility_score(query_annotations, pr2_annotations)
		if score > best_score:
			if debug:
				print(pr2_annotations)
				print(score, common_ranks)
			pr2_genus = pr2_annotations["genus"]
			best_match = row
			best_score = score
	if query_genus != pr2_genus:
		return best_match, 0
	else:
		return best_match, best_score


def _evaluate(identifier, best_match, best_score, ncbilineage):
	# Evaluate best score and report failed searches
	with open("lineages_check.log", "at") as errfile:
		if best_score > 0:
			# Decent chance we have a PR2 lineage match
			lineage = ";".join(list(best_match))
		else:
			# We need to create a new PR2 lineage item;
			# best_match is originally a named Series, need to convert to list
			matching_sp = list(best_match)[-1] # PR2 species matching the query
			best_match = list(best_match)[:-2] # PR2 lineage matching the query
			best_match.append(identifier.split()[0]) #query genus
			best_match.append(identifier.replace(" ", "_")) #query species_
			lineage = ";".join(best_match) # suggested lineage for manual curation
			errfile.write(f"{identifier}\t{_clean_name(identifier)}\t<{matching_sp}\n{lineage}\n{ncbilineage}\n")
	
	# Inform about the evaluation
	if best_match is not None:
		print("Best Match:", identifier, lineage)
		#print(best_match)
		print(f"Compatibility Score: {best_score}")
	else:
		print("No match found in PR2 database.")
		
	return lineage
	

def _get_known_taxonomies():
	#these are the available semi-manually parsed taxonomies
	taxonomyfile = datadir + "phylodb_2.0.species_lineages.tsv"
	taxonomy = {}
	with open(taxonomyfile) as f:
		for l in f:
			d = l.strip().split("\t")
			if len(d) < 7:
				if len(l) > 1:
					print("error parsing", l)
				continue
			#species	species_name	species_NCBI	lineage_8	lineage_9	lineage_NCBI	taxid
			#[0]		[1]				[2]				[3]			[4]			[5]				[6]
			taxonomy[d[0]] = d
	return taxonomy


def _get_NCBI_taxonomy_info(query):
	def fetch_taxid(query):
		try:
			name_to_taxid = ncbi.get_name_translator([query])
			if not name_to_taxid:
				return None
			return name_to_taxid[query][0]
		except Exception:
			return None

	def unique_tuples(tuples_list):
		seen = set()
		unique_list = [t for t in tuples_list if not (t in seen or seen.add(t))]
		return unique_list

	# Try to get taxid for the full query
	taxid = fetch_taxid(query)
	if not taxid:
		# Try to get taxid for species (first two words)
		species_query = ' '.join(query.split()[:2])
		taxid = fetch_taxid(species_query)
		if not taxid:
			# Try to get taxid for genus (first word)
			genus_query = query.split()[0]
			taxid = fetch_taxid(genus_query)
			if not taxid:
				return None, None, f"Could not find taxid for query: {query}"

	try:
		lineage = ncbi.get_lineage(taxid)
		lineage_ranks = ncbi.get_rank(lineage)
		lineage_names = ncbi.get_taxid_translator(lineage)
		
		#this is now ordered
		filtered_lineage = [
			(lineage_names[taxid], lineage_ranks[taxid])
			for taxid in lineage
			if lineage_ranks[taxid] in valid_ranks
		]
		
		# Include the query itself if it belongs to a valid rank
		# ChatGPT suggestion - probably redundant
		#if lineage_ranks[taxid] in valid_ranks:
		#	filtered_lineage.append((lineage_names[taxid], lineage_ranks[taxid]))
		
		if debug:
			print("NCBI filtered lineage (names, ranks):", filtered_lineage)
		
		return taxid, unique_tuples(filtered_lineage), None
	except Exception as e:
		return None, None, str(e)


def _get_pr2_lineages(identifier, filtered_lineage):
	#read PR2 taxonomies
	pr2_8ranks = pd.read_csv(datadir + "PR2_4.14.0_taxonomies_8-level-up.tsv", sep='\t')
	pr2_9ranks = pd.read_csv(datadir + "PR2_5.0.0_taxonomies_9-level-up.tsv", sep='\t')
	
	#convert NCBI lineages to a dictionary and a set
	query_annotations = {}
	for x in filtered_lineage:
		if x[1] not in query_annotations:
			query_annotations[x[1]] = x[0].replace(" ", "_") 
	query_lineage = ";".join([x[0].replace(" ", "_") for x in filtered_lineage])
	query_set = {x[0].replace(" ", "_") for x in filtered_lineage}
	
	#handle non-existent taxid:
	if not query_annotations:
		quit(f"taxid non-existent: {identifier}")
	if "species" in query_annotations:
		if query_annotations["species"] != identifier.replace(" ", "_"):
			# Replaced spaces to align with above query_annotations definition;
			# test if matches submitted identifier
			new_taxid = query_annotations["species"]
			print(f"taxid species disagreement: {identifier} => {new_taxid}")
		else:
			print(f"taxid species OK: {identifier}")
			pass # this is okay
		# if species in query_annotations, there will be a genus too
		# use genus to identify PR2 lineages
		query_clade = query_annotations["genus"] 
	else:
		if debug:
			print("NCBI annotations dict (only first of redundant ranks!)", query_annotations)
		if "genus" in query_annotations:
			query_clade = query_annotations["genus"]
			print(f"taxid missing species: {identifier}; genus => {query_clade}")
		else:
			if "family" in query_annotations:
				query_clade = query_annotations["family"]
				print(f"taxid missing genus: {identifier}; family => {query_clade}")
			else:
				quit(f"taxid missing species, genus and family: {identifier}; quitting") #need genus for _find_best_match
	
	best_8rank_match, best_8rank_score = _find_best_match(query_set, query_clade, pr2_8ranks)
	best_9rank_match, best_9rank_score = _find_best_match(query_set, query_clade, pr2_9ranks)
	
	lineage_8 = _evaluate(identifier, best_8rank_match, best_8rank_score, query_lineage)
	lineage_9 = _evaluate(identifier, best_9rank_match, best_9rank_score, query_lineage)
	
	return lineage_8, lineage_9


def _move_to_temp(file_path, directory_path):
	# Check if the directory exists
	if os.path.exists(directory_path):
		# Move the file into the directory
		shutil.move(file_path, directory_path)
		print(f"File moved to {directory_path}")
	else:
		print(f"Directory {directory_path} does not exist")


def _reformat_fasta(filename, output_prefix):
	output_filename = f"{output_prefix}_phylodb.faa"
	with open(output_filename, 'w') as output_file:
		# Iterate through the sequences in the input file
		for seq in SeqIO.parse(filename, "fasta"):
			# Extract the sequence ID from the header
			seq_id = seq.name
			# Read the sequence data
			sequence = str(seq.seq)
			# Calculate the MD5 hash of the sequence
			md5seq = hashlib.sha256(sequence.encode()).hexdigest()
			# Write the reformatted header+sequence
			output_file.write(f">{seq_id}\t{md5seq}\t{output_prefix}\n{sequence}\n")


if __name__ == "__main__":
	args = parse_arguments()
	_def_variables()
	
	# Join the species arguments into a single string with spaces
	sciname = ' '.join(args.species)
	handle = _clean_name(sciname)
	if sciname.startswith("Candidatus"):
		genus = args.species[1] #sciname.split(" ")[1]
		species = " ".join(args.species[1:3])
		strain = " ".join(args.species[3:])
	else:
		genus = args.species[0] #sciname.split(" ")[0]
		species = " ".join(args.species[:2])
		strain = " ".join(args.species[2:])
	
	# Call the reformat_fasta function with the provided arguments
	_reformat_fasta(args.filename, handle)

	# Make a new item in the taxonomy file
	taxonomy = _get_known_taxonomies()
	with open("new_handles-taxonomy.tsv", "at") as result:
		if species in taxonomy:
			d = taxonomy[species]
			#handle	species	species_name	species_NCBI	taxid	strain	lineage_8	lineage_9	lineage_NCBI	lineage_phylodb
			result.write(f"{handle}\t{species}\t{d[1]}\t{d[2]}\t{d[6]}\t{strain}\t{d[3]}\t{d[4]}\t{d[5]}\tNA\n")
		elif species.endswith("sp.") and genus in taxonomy:
			d = taxonomy[genus]
			#handle	species	species_name	species_NCBI	taxid	strain	lineage_8	lineage_9	lineage_NCBI	lineage_phylodb
			result.write(f"{handle}\t{species}\t{d[1]}\t{d[2]}\t{d[6]}\t{strain}\t{d[3]}\t{d[4]}\t{d[5]}\tNA\n")
		else:
			print("taxonomy not found in reference file")
			#write a function to obtain data from sciname/species/genus:
			#taxid
			#genus;species NCBI
			#lineage_NCBI
			#PR2 lineages
			outfile = open("new_species_lineages.tsv", "at")
			taxid, filtered_lineage, error = _get_NCBI_taxonomy_info(sciname)

			if error:
				outfile.write(f"{query}\tError:{error}\n\n")
				quit(f"{query}\tError:{error}\n\n")
				
			if filtered_lineage:
				# No need to sort further; already properly sorted by the function				
				species_ = species.replace(" ", "_")
				lineage_NCBI = []
				for name, rank in filtered_lineage:
					lineage_NCBI.append(f"{rank.capitalize()}:{name}")
				species_NCBI = ";".join([genus, species])
				lineage_NCBI = ";".join(lineage_NCBI)
				lineage_8, lineage_9 = _get_pr2_lineages(sciname, filtered_lineage)
				
				# Write the results to the output file
				#species	species_name	species_NCBI	lineage_8	lineage_9	lineage_NCBI	taxid
				#[0]		[1]				[2]				[3]			[4]			[5]				[6]
				outfile.write(f"{species}\t{species_}\t{species_NCBI}\t{lineage_8}\t{lineage_9}\t{lineage_NCBI}\t{taxid}\n")
				
				# Write lineages to result file
				result.write(f"{handle}\t{species}\t{species_}\t{species_NCBI}\t{taxid}\t{strain}\t{lineage_8}\t{lineage_9}\t{lineage_NCBI}\tNA\n")
				
				#print(f"lineage length {query}: {len(filtered_lineage)}")
				
			else:
				outfile.write(f"{query}\tNo taxonomy information found.\n")
			
			outfile.close()
				
			#result.write(f"{species}\t{handle}\t\n")
	if not debug:
		_move_to_temp(args.filename, "tmp_original")
