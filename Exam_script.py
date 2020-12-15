#!/usr/bin/python3

import os, sys, subprocess, shutil
import re
import string
from Bio.Seq import Seq
import pandas as pd
from Bio import SeqIO



## user input will be put in a library 
details={}
details["organism"]   = input(" Which organism do you want to make a BLAST database for? ")
details["gene"]   = input("\n Which gene do you want to search for? (leave empty if no specific gene is wanted) ")
details["db"]    = input("\n Do you want to make a BLAST database for protein or nucleotide?\n Please answer p or n ")
details["complete"] = input ("\n Do you want to include partial/incomplete result?\n Please answer yes or no ")


## yes and no function which return True or False to use for conditions 
def yes_no(answer):
	yes = set(['yes','y'])		# yes and y are put into a set; both yes or y will return True
	no = set(['no','n']) 		# no and n are put into a set; both no or n will return False
	choice = answer.lower()		# set all types of answer to lower case
	# loop forever until something is returned
	while True:		
		if choice in yes:		# both yes or y will return True
			return True
		elif choice in no:		# both no or n will return False
			return False
		else:					# error trap, all other input will cause the function to ask "yes or no" over and over until a desired answer is received
			print ("\n Please respond with 'yes' or 'no' \n")
			choice = input(" yes or no?").lower()


## function which return True or False to use for conditions 
def protein_nt(answer):
	prot = set(['protein','p'])		# protein and p are put into a set; both protein or p will return True
	nucl = set(['nucleotide','n']) 		# nucleotide and n are put into a set; both nucleotide or n will return False
	choice = answer.lower()		# set all types of answer to lower case
	# loop forever until something is returned
	while True:		
		if choice in prot:		# both protein or p will return True
			return True
		elif choice in nucl:		# both nucleotide or n will return False
			return False
		else:					# error trap, all other input will cause the function to ask "p or n" over and over until a desired answer is received
			print ("\n Please respond with 'p' or 'n' \n")
			choice = input(" protein or nucleotide?").lower()
			

# pass the library values variables
organism,gene,db,complete = list(details.values())

# check user inputs are valid and define the search term
if len(organism)==0 or str.isspace(organism):		# check organism is not empty or spaces 
	print("\n Sorry, no organism was chosen, please try again..\n") 	# error message; error trap
elif len(gene)==0 or str.isspace(gene):				# check if gene is empty 
	print("\n Ok, no specific gene will be searched for. \n")
	if protein_nt(db): 		# if user choose protein database
		if yes_no(complete):		# if user choose to include partial search
			print("\n Search will be done on following information:\n\tOrganism:",organism,"\n\tDatabase: protein","\n\tInclude partial/incomplete: Yes")
			es = "esearch -db protein -query \" "+ organism +"[organism] \" "		# prepare for esearch
			ef = "|efetch -db protein -format fasta"		# prepare for efetch
			ff = ".prot.fa"		# prepare for output file name 
		else:		# if user choose to not include partial search
			print("\n Search will be done on following information:\n\tOrganism:",organism,"\n\tDatabase: protein","\n\tInclude partial/incomplete: No")
			es = "esearch -db protein -query \" "+ organism +"[organism] Not partial \" "
			ef = "|efetch -db protein -format fasta"
			ff = ".prot.fa"
	else: 		# if user choose nucleotide database
		if yes_no(complete):		# if user choose to include partial search
			print("\n Search will be done on following information:\n\tOrganism:",organism,"\n\tDatabase: nucleotide","\n\tInclude partial/incomplete: Yes")
			es = "esearch -db nucleotide -query \" "+ organism +"[organism] \" "
			ef = "|efetch -db nucleotide -format fasta"
			ff = ".nuc.fa"
		else:		# if user choose to include partial search
			print("\n Search will be done on following information:\n\tOrganism:",organism,"\n\tDatabase: nucleotide","\n\tInclude partial/incomplete: No")
			es = "esearch -db nucleotide -query \" "+ organism +"[organism] Not partial \" "
			ef = "|efetch -db nucleotide -format fasta"
			ff = ".nuc.fa"
else:		# if gene is not empty 
	if protein_nt(db):  		# if user choose protein database
		if yes_no(complete):		# if user choose to include partial search
			print("\n Search will be done on following information:\n\tOrganism:",organism,"\n\tGene:",gene,"\n\tDatabase: protein","\n\tInclude partial/incomplete: Yes")
			es = "esearch -db protein -query \" "+ organism +"[organism] AND " + gene +"[Gene name] \" "
			ef = "|efetch -db protein -format fasta"
			ff = ".prot.fa"
		else:		# if user choose to not include partial search
			print("\n Search will be done on following information:\n\tOrganism:",organism,"\n\tGene:",gene,"\n\tDatabase: protein","\n\tInclude partial/incomplete: No")
			es = "esearch -db protein -query \" "+ organism +"[organism] AND " + gene +"[Gene name] Not partial \" "
			ef = "|efetch -db protein -format fasta"
			ff = ".prot.fa"
	else: 		# if user choose nucleotide database
		if yes_no(complete):		# if user choose to include partial search
			print("\n Search will be done on following information:\n\tOrganism:",organism,"\n\tGene:",gene,"\n\tDatabase: nucleotide","\n\tInclude partial/incomplete: Yes")
			es = "esearch -db nucleotide -query \" "+ organism +"[Organism] AND " + gene +"[Gene name] Not partial \" "
			ef = "|efetch -db nucleotide -format fasta"
			ff = ".nuc.fa"
		else:		# if user choose to include partial search
			print("\n Search will be done on following information:\n\tOrganism:",organism,"\n\tGene:",gene,"\n\tDatabase: nucleotide","\n\tInclude partial/incomplete: No")
			es = "esearch -db nucleotide -query \" "+ organism +"[Organism] AND " + gene +"[Gene name] Not partial \" "
			ef = "|efetch -db nucleotide -format fasta"
			ff = ".nuc.fa"

# search function to search and download sequences
def search(organism,gene,db,complete):
	if es:
		es_number = es + "|grep -i \"count\"|awk \'{split($0,a,\"<|>\");print a[3];}\'" 		# grep and awk to get esearch result number
		print("\n This is what is running in shell: \n\n " + es + "\n\n Please wait... \n")
		seq_number = subprocess.check_output(es_number,shell=True)		# get the search result number from shell 
		if int(seq_number) > 1000:		# if search result is over 1000, print warning message, end the script
			print("\n ** Warning: Over 1000 sequences found, continue is not recommended, please narrow down your search,"+ \
				"\n otherwise very slow processing speed and probably taking too much space! Thank you! \n")
			quit()
		if int(seq_number) == 0:		# if search result is 0, print warning message, end the script
			print("\n Sorry, no sequence was found! Likely spelling mistakes. Please try again. Thank you! \n")
			quit()
		else:
			# print amount of the sequence found, nubmer converted from byte string to normal string
			print("\n------\n "+ str(seq_number.decode('ascii').rstrip()) +" sequences was found! Nice choice! \n\n")
			
			# provide choices for user by ask if they want to download the sequences or not
			dow = input(" Do you want to download the sequences on your server? (Please note: organism name will be used as output file name.) \n\n BLAST database will only be available if sequences are downloaded. \n\n If your changed your mind about your search, you can reply no and start over again. \n\n Otherwise, please respond yes. \n\n Please respond yes or no. ")
			# yes_no function: only start downloading if response returns True
			if yes_no(dow):
				print("\n------\n Downloading sequences...\n\n Please wait... \n")
				# output file name is based on user's taxon input
				dow_file_name = '_'.join(organism.split())
				# download sequence with efetch and save sequances in file_name based on user's taxon input
				efet = es + ef + " >"+ dow_file_name + ff
				# print processing message
				print("\n This is what is running in shell: \n\n " + efet + "\n\n Please wait... \n")
				# call download in shell
				subprocess.call(efet,shell=True)
				# open the downloaded file and and confirm protein sequence number by counting ">" 
				print ("\n------\n Sequence downloaded! Checking " + dow_file_name + ff + " for content... \n")  
				file_contents = open(dow_file_name + ff).read()
				count = file_contents.count('>')
				# print confirmation message
				print ("\n------\n Check completed. " + str(count) + " sequences were successfully retrieved! Sequences are saved in " \
				+ dow_file_name + ff + " \n" )
			# if user do not want to continue downloading sequences, quit the script
			else:
				print("\n Thank you for searching! Bye! \n")
				quit()				


## function to check how many species there are if protein sequences were downloaded
def species_number(organism):
	# file name 
	dow_file_name = '_'.join(organism.split())
	# only carry on if the file exists
	if os.path.isfile(dow_file_name + ".prot.fa"):
		# read downloaded sequences
		file_contents = open(dow_file_name + ".prot.fa").read()
		# sequence number count 
		count = file_contents.count('>')
		# count unique species number, using regex to find all the text with [] around
		spe = list(set(re.findall('\[.*?\]',file_contents)))
		genus = list(set(re.findall('\[\w*',file_contents)))
		# if more than one species in the downloaded file do following
		if len(spe) >1:
			#print species count and ask if user want to continue or not
			print("\n Sequences are from " + str(len(genus))+ " different genera and there are " + str(len(spe)) + " different species in total. Do you wish to continue? \n")
			answer = input(" Please respond yes or no. ")
			if yes_no(answer)==False :
				print("\n Thank you for searching! Bye! \n")
		else:
			print("\n Sequences are from the same species " + str(spe[0].strip("[]")) + ". ")
			

##function to make blast database
def blastdb(organism,gene,db,complete):
	# downloaded file name
	dow_file_name = '_'.join(organism.split())
	# check the file exists
	if os.path.isfile(dow_file_name + ff):
		# ask user for answer
		blast_q = input("\n------\n Do you wish to make a BLAST database based on downloaded sequences? \n\n Please respond yes or no. ")
		# evaluate answer by yew_no function, carry on if yes
		if yes_no(blast_q):
			# if nucleotide database
			if ff == ".nuc.fa":
				print("\n------\n Preparing a nucleotide database... ") 
				# make database in shell
				database= "makeblastdb -in " + dow_file_name + ff +" -dbtype nucl -out " + dow_file_name
				subprocess.call(database,shell=True)
				# check output and print success message
				if os.path.isfile(dow_file_name + ".nhr"):
					print("\n------\n Nucleotide database has been successfully made! \n\n Nucleotide sequence headers are in " + dow_file_name + ".nhr file.\n Nucleotide indexes are in " + dow_file_name + ".nin file.\n Compressed nucleotide sequences are in " + dow_file_name + ".nsq file.\n")
			else:	# if protein database
				print("\n------\n Preparing a protein database... ") 
				# make database in shell
				database= "makeblastdb -in " + dow_file_name + ff +" -dbtype prot -out " + dow_file_name
				subprocess.call(database,shell=True)
				# check output and print success message
				if os.path.isfile(dow_file_name + ".phr"):
					print("\n------\n Protein database has been successfully made! \n\n Protein sequence headers are in " + dow_file_name + ".phr file.\n Protein indexes are in " + dow_file_name + ".pin file.\n Compressed protein sequences are in " + dow_file_name + ".psq file.\n")
		#  user does not want a database
		else:
			print("\n Sorry to know that you do not want a database. Bye! \n") 


## function to do blast analysis			
def blast(organism,gene,db,complete):
	# downloaded file name
	dow_file_name = '_'.join(organism.split())
	# only carry on if a database is made
	if os.path.isfile(dow_file_name + ".nhr") or os.path.isfile(dow_file_name + ".phr"):
		# ask user for input
		blast_in = input("\n------\n Do you wish to do BLAST for downloaded search results to find the most similar sequences?\n Please answer yes or no ")
		# yes_no function , carry on if yes
		if yes_no(blast_in):
			# if nucleotide database made
			if ff == ".nuc.fa":
				# read the downloaded sequence using biopython 
				records = SeqIO.parse("./"+ dow_file_name + ".nuc.fa","fasta")
				# ask user for restriction for output
				hit_number = input("\n How many top hits would you like for each sequence? ")
				# check if user input is invalid e.g 0 or white space or not a digit
				if str(hit_number) == '0' or str.isspace(str(hit_number)) or str(hit_number).isdigit()== False:
					print("\n The input is invalid. Resetting value to 10... \n")
					print("\n Processing blastn analysis... Please wait....\n ") 
					# reset the value to 10
					hit_number = 10 
					# loop over the downloaded file
					for i in records:
						acc = i.id	# accessions are ids of records
						single_seq = open('single_seq.fasta', 'w')	# open a new file and write; next loop will overwrite the file 
						single_seq.write(str(i.seq))	# write a single sequence into the file
						single_seq.close()	# close the file connection
						# nucleotide sequence against nucleotide database use blastn
						blastn = "blastn -db " + dow_file_name + " -query single_seq.fasta -outfmt 7 >> blast.out"	
						# call shell to do the blast 
						subprocess.call(blastn,shell=True)
						# trim all lines with '#' or self hit, then take the top 10 hit and write into a tsv file 
						grep = "grep -v \"#\\|" + str(acc) + '" ' + " blast.out | head -n" + str(hit_number) + ">> blast.tsv"
						subprocess.call(grep,shell=True)
						# replace the query with accession of the query sequence
						change = "sed -i \'s/Query_1/" + str(acc) + "/g\' blast.tsv"
						subprocess.call(change,shell=True)
					# set df1 to none type for later to check for change
					df1 = None
					# read tsv to df1 with no header and tab-delimiter  
					df1 = pd.read_csv('./blast.tsv',sep="\t",header=None)
					# set header
					df1.columns=['query acc.', 'subject acc.', '% identity', 'alignment length', 'mismatches', 'gap opens', 'q. start', 'q. end', 's. start', 's. end', 'e_value', 'bit score']
					# reorder the data frame by e-value and bit score
					df1.sort_values(['e_value','bit score'], ascending=True, inplace=True)
					df1.to_csv(r'./hitresult.csv',sep='\t')
					# check df1 is changed and print success messages
					if df1 is not None:
						print("\n Dataframe containing all hit results is ready in hitresult.csv file. The similarity of sequences is the highest at the top. \n")
						print("\n Showing first 10 lines of the hitresult.csv file: \n ")
						# show first 10 lines of the blast output
						show = " head -n10 hitresult.csv "
						subprocess.call(show,shell=True)				
				# user input for top hit number is valid
				else:
					print("\n Processing blastn analysis... Please wait....\n ") 
					# loop over the downloaded file					
					for i in records:
						acc = i.id	# accessions are ids of records
						single_seq = open('single_seq.fasta', 'w')	# open a new file and write; next loop will overwrite the file 
						single_seq.write(str(i.seq))	# write a single sequence into the file
						single_seq.close()	# close the file connection
						# nucleotide sequence against nucleotide database use blastn
						blastn = "blastn -db " + dow_file_name + " -query single_seq.fasta -outfmt 7 >> blast.out"	
						# call shell to do the blast 
						subprocess.call(blastn,shell=True)
						# trim all lines with '#' or self hit, then take the top 10 hit and write into a tsv file
						grep = "grep -v \"#\\|" + str(acc) + '" ' + " blast.out | head -n" + str(int(hit_number)) + ">> blast.tsv"
						subprocess.call(grep,shell=True)
						# replace the query with accession of the query sequence
						change = "sed -i \'s/Query_1/" + str(acc) + "/g\' blast.tsv"
						subprocess.call(change,shell=True)
					# set df1 to none type for later to check for change
					df1 = None
					# read tsv to df1 with no header and tab-delimiter  
					df1 = pd.read_csv('./blast.tsv',sep="\t",header=None)
					# set header
					df1.columns=['query acc.', 'subject acc.', '% identity', 'alignment length', 'mismatches', 'gap opens', 'q. start', 'q. end', 's. start', 's. end', 'e_value', 'bit score']
					# reorder the data frame by e-value and bit score
					df1.sort_values(['e_value','bit score'], ascending=True, inplace=True)
					df1.to_csv(r'./hitresult.csv',sep='\t')
					# check df1 is changed and print success messages
					if df1 is not None:
						print("\n Dataframe containing all hit results is ready in hitresult.csv file. The similarity is the highest for comparison sequences at the top. \n")
						print("\n Showing first 10 lines of the hitresult.csv file: \n ")
						# show first 10 lines of the blast output
						show = " head -n10 hitresult.csv "
						subprocess.call(show,shell=True)
			# if protein database made
			else:
				print("\n Protein BLAST analysis can be done similarly as nucleotide but using blastp instead of blastn \n")
				records = SeqIO.parse("./"+ dow_file_name + ".prot.fa","fasta")
				hit_number = input("\n How many top hits would you like for each sequence? ")
				# check if user input is invalid e.g 0 or white space or not a digit
				if str(hit_number) == '0' or str.isspace(str(hit_number)) or str(hit_number).isdigit()== False:
					print("\n The input is invalid. Resetting value to 10... \n")
					print("\n Processing blastp analysis... Please wait....\n ") 
					# reset the value to 10
					hit_number = 10 
					# loop over the downloaded file
					for i in records:
						acc = i.id	# accessions are ids of records
						single_seq = open('single_seq.fasta', 'w')	# open a new file and write; next loop will overwrite the file 
						single_seq.write(str(i.seq))	# write a single sequence into the file
						single_seq.close()	# close the file connection
						# protein sequence against protein database use blastp
						blastn = "blastp -db " + dow_file_name + " -query single_seq.fasta -outfmt 7 >> blast.out"	
						# call shell to do the blast 
						subprocess.call(blastn,shell=True)
						# trim all lines with '#' or self hit, then take the top 10 hit and write into a tsv file 
						grep = "grep -v \"#\\|" + str(acc) + '" ' + " blast.out | head -n" + str(hit_number) + ">> blast.tsv"
						subprocess.call(grep,shell=True)
						# replace the query with accession of the query sequence
						change = "sed -i \'s/Query_1/" + str(acc) + "/g\' blast.tsv"
						subprocess.call(change,shell=True)
					# set df1 to none type for later to check for change
					df1 = None
					# read tsv to df1 with no header and tab-delimiter  
					df1 = pd.read_csv('./blast.tsv',sep="\t",header=None)
					# set header
					df1.columns=['query acc.', 'subject acc.', '% identity', 'alignment length', 'mismatches', 'gap opens', 'q. start', 'q. end', 's. start', 's. end', 'e_value', 'bit score']
					# reorder the data frame by e-value and bit score
					df1.sort_values(['e_value','bit score'], ascending=True, inplace=True)
					df1.to_csv(r'./hitresult.csv',sep='\t')
					# check df1 is changed and print success messages
					if df1 is not None:
						print("\n Dataframe containing all hit results is ready in hitresult.csv file. The similarity of sequences is the highest at the top. \n")
						print("\n Showing first 10 lines of the hitresult.csv file: \n ")
						# show first 10 lines of the blast output
						show = " head -n10 hitresult.csv "
						subprocess.call(show,shell=True)
				# user input for top hit number is valid
				else:
					print("\n Processing blastp analysis... Please wait....\n ") 
					# loop over the downloaded file					
					for i in records:
						acc = i.id	# accessions are ids of records
						single_seq = open('single_seq.fasta', 'w')	# open a new file and write; next loop will overwrite the file 
						single_seq.write(str(i.seq))	# write a single sequence into the file
						single_seq.close()	# close the file connection
						# protein sequence against protein database use blastp
						blastn = "blastp -db " + dow_file_name + " -query single_seq.fasta -outfmt 7 >> blast.out"	
						# call shell to do the blast 
						subprocess.call(blastn,shell=True)
						# trim all lines with '#' or self hit, then take the top 10 hit and write into a tsv file
						grep = "grep -v \"#\\|" + str(acc) + '" ' + " blast.out | head -n" + str(int(hit_number)) + ">> blast.tsv"
						subprocess.call(grep,shell=True)
						# replace the query with accession of the query sequence
						change = "sed -i \'s/Query_1/" + str(acc) + "/g\' blast.tsv"
						subprocess.call(change,shell=True)
					# set df1 to none type for later to check for change
					df1 = None
					# read tsv to df1 with no header and tab-delimiter  
					df1 = pd.read_csv('./blast.tsv',sep="\t",header=None)
					# set header
					df1.columns=['query acc.', 'subject acc.', '% identity', 'alignment length', 'mismatches', 'gap opens', 'q. start', 'q. end', 's. start', 's. end', 'e_value', 'bit score']
					# reorder the data frame by e-value and bit score
					df1.sort_values(['e_value','bit score'], ascending=True, inplace=True)
					df1.to_csv(r'./hitresult.csv',sep='\t')
					# check df1 is changed and print success messages
					if df1 is not None:
						print("\n Dataframe containing all hit results is ready in hitresult.csv file. The similarity is the highest for comparison sequences at the top. \n")			
						print("\n Showing first 10 lines of the hitresult.csv file: \n ")
						# show first 10 lines of the blast output
						show = " head -n10 hitresult.csv "
						subprocess.call(show,shell=True)		

		
		# ask user if they wish to use local sequence for blast	
		answer = input("\n------\n Do you wish to do a BLAST analysis with your local sequence file? Please note only FASTA file is accepted.\n Please answer yes or no ")
		# yes_not function; if yes, carry on
		if yes_no(answer):
			blast_ex = input("\n What is you local file name? Please note the file has to be in the same directory where script is running. \n ")
			# check the end of the file name for file format
			if blast_ex.endswith(('fasta','fa')):
				# check the file exists in the same directory as the script
				if os.path.isfile('./' + blast_ex): 
					# if nucleotide database
					if ff == ".nuc.fa":
						print("\n Your local file will be blast against the nucleotide database made in earlier step.\n")
						# read user file 
						file_contents = open(blast_ex).read()
						# count sequence number in the file
						count = file_contents.count('>')
						if int(count) > 0:
							print("\n File loaded successfully, your file contains " + count + " sequences \n " )
							# ask user if the file contains nucleotide sequences or protein sequences
							p_or_n = input("\n Do your file contains protein sequences or nucleotide sequences? \n Please answer p or n ")
							# if contains protein sequences
							if protein_nt(p_or_n):
								print("\n Your protein sequences will be blast against the nucleotide database using tblastn \n ")
								#### similar to above
								#### skipped for easier reading
								
							# contains nucleotides 	
							else:
								print("\n Your nucleotide sequences will be blast against the nucleotide database using blastn \n ")
								#### similar to above
								#### skipped for easier reading
																
						# no sequence found in the file, error message
						else:
							print("\n There is no sequence in your file. Analysis cannot be done, sorry! \n ")
					# if protein database
					else:
						print("\n Your local file will be blast against the protein database made in earlier step.\n")
						# read user file
						file_contents = open(blast_ex).read()
						# count sequence number in the file
						count = file_contents.count('>')
						if int(count) > 0:
							print("\n File loaded successfully, your file contains " + count + " sequences \n ")
							# ask user if the file contains nucleotide sequences or protein sequences
							p_or_n = input("\n Do your file contains protein sequences or nucleotide sequences? \n Please answer p or n ")
							# if contains protein sequences
							if protein_nt(p_or_n):
								print("\n Your protein sequences will be blast against the protein database using blastp \n ")
								#### similar to above
								#### skipped for easier reading
								
							# contains nucleotides 	
							else:
								print("\n Your nucleotide sequences will be blast against the protein database using blastx \n ")
								#### similar to above
								#### skipped for easier reading
								
						# no sequence found in the file, error message
						else:
							print("\n There is no sequence in your file. Analysis cannot be done, sorry! \n ")							
				# file not in the directory, error message
				else:
					print("\n Your file is not in the same directory with the script. Analysis cannot be done, sorry! \n ")  	
			# file has wrong format, error message 
			else:
				print("\n The file is wrong format. Analysis cannot be done, sorry! \n ")
		# user do not want to use this function
		else:
			print("\n Thank you for searching! Bye! \n")
			quit()

# call search function and pass multiple arguments from dictionary 
search(*list(details.values()))
# call function to count the species number if protein sequences are downloaded; pass the first argument of the dictionary to the function
species_number(list(details.values())[0])
# call blast make database function and pass multiple arguments from dictionary 
blastdb(*list(details.values()))
# call blast analysis function and pass multiple arguments from dictionary 
blast(*list(details.values()))

