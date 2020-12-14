#!/usr/bin/python3

import os, sys, subprocess, shutil
import re
import numpy as np
from pandas.core.common import flatten
import string
import glob
from Bio.Seq import Seq

from Bio import Entrez


## user input will be put in a library 
details={}
details["organism"]   = input(" Which organism do you want to make a BLAST database for? ")
details["gene"]   = input("\n Which gene do you want to search for? (leave empty if no specific gene is wanted) ")
details["db"]    = input("\n Do you want to make a BLAST database for protein or nucleotide?\n Please answer p or n ")
details["complete"] = input ("\n Do you want to include partial/incomplete result?\n Please answer yes or no ")

## yes and no function which return True or False to use for conditions 
def yes_no(answer):
	yes = set(['yes','y'])		# both yes or y will return True
	no = set(['no','n']) 		# both no or n will return False
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
	yes = set(['protein','p'])		# both protein or p will return True
	no = set(['nucleotide','n']) 		# both nucleotide or n will return False
	choice = answer.lower()		# set all types of answer to lower case
	# loop forever until something is returned
	while True:		
		if choice in yes:		# both yes or y will return True
			return True
		elif choice in no:		# both no or n will return False
			return False
		else:					# error trap, all other input will cause the function to ask "yes or no" over and over until a desired answer is received
			print ("\n Please respond with 'p' or 'n' \n")
			choice = input(" protein or nucleotide?").lower()
			
# pass the library values variables
organism,gene,db,complete = list(details.values())

# check user inputs are valid and define the search term
if len(organism)==0 or str.isspace(organism):
	print("\n Sorry, no organism was chosen, please try again..\n") 
elif len(gene)==0 or str.isspace(gene):
	print("\n Ok, no specific gene will be searched for. \n")
	if protein_nt(db): 
		if yes_no(complete):
			print("\n Search will be done on following information:\n\tOrganism:",organism,"\n\tDatabase: protein","\n\tInclude partial/incomplete: Yes")
			es = "esearch -db protein -query \" "+ organism +"[organism] \" "
			ef = "|efetch -db protein -format fasta"
			ff = ".prot.fa"
		else:
			print("\n Search will be done on following information:\n\tOrganism:",organism,"\n\tDatabase: protein","\n\tInclude partial/incomplete: No")
			es = "esearch -db protein -query \" "+ organism +"[organism] Not partial \" "
			ef = "|efetch -db protein -format fasta"
			ff = ".prot.fa"
	else:
		if yes_no(complete):
			print("\n Search will be done on following information:\n\tOrganism:",organism,"\n\tDatabase: nucleotide","\n\tInclude partial/incomplete: Yes")
			es = "esearch -db nucleotide -query \" "+ organism +"[organism] \" "
			ef = "|efetch -db nucleotide -format fasta"
			ff = ".nuc.fa"
		else:
			print("\n Search will be done on following information:\n\tOrganism:",organism,"\n\tDatabase: nucleotide","\n\tInclude partial/incomplete: No")
			es = "esearch -db nucleotide -query \" "+ organism +"[organism] Not partial \" "
			ef = "|efetch -db nucleotide -format fasta"
			ff = ".nuc.fa"
else:
	if protein_nt(db): 
		if yes_no(complete):
			print("\n Search will be done on following information:\n\tOrganism:",organism,"\n\tGene:",gene,"\n\tDatabase: protein","\n\tInclude partial/incomplete: Yes")
			es = "esearch -db protein -query \" "+ organism +"[organism] AND " + gene +"[Gene name] \" "
			ef = "|efetch -db protein -format fasta"
			ff = ".prot.fa"
		else:
			print("\n Search will be done on following information:\n\tOrganism:",organism,"\n\tGene:",gene,"\n\tDatabase: protein","\n\tInclude partial/incomplete: No")
			es = "esearch -db protein -query \" "+ organism +"[organism] AND " + gene +"[Gene name] Not partial \" "
			ef = "|efetch -db protein -format fasta"
			ff = ".prot.fa"
	else:
		if yes_no(complete):
			print("\n Search will be done on following information:\n\tOrganism:",organism,"\n\tGene:",gene,"\n\tDatabase: nucleotide","\n\tInclude partial/incomplete: Yes")
			es = "esearch -db nucleotide -query \" "+ organism +"[Organism] AND " + gene +"[Gene name] Not partial \" "
			ef = "|efetch -db nucleotide -format fasta"
			ff = ".nuc.fa"
		else:
			print("\n Search will be done on following information:\n\tOrganism:",organism,"\n\tGene:",gene,"\n\tDatabase: nucleotide","\n\tInclude partial/incomplete: No")
			es = "esearch -db nucleotide -query \" "+ organism +"[Organism] AND " + gene +"[Gene name] Not partial \" "
			ef = "|efetch -db nucleotide -format fasta"
			ff = ".nuc.fa"

def search(organism,gene,db,complete):
	if es:
		es_number = es + "|grep -i \"count\"|awk \'{split($0,a,\"<|>\");print a[3];}\'"
		print("\n This is what is running in shell: \n\n " + es + "\n\n Please wait... \n")
		seq_number = subprocess.check_output(es_number,shell=True)
		if int(seq_number) > 1000:
			print("\n ** Warning: Over 1000 sequences found, continue is not recommended, please narrow down your search,"+ \
				"\n otherwise very slow processing speed and probably taking too much space! Thank you! \n")
			quit()
		if int(seq_number) == 0:
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

def blastdb(organism,gene,db,complete):
	dow_file_name = '_'.join(organism.split())
	if os.path.isfile(dow_file_name + ff):
		blast_q = input(" Do you wish to make a BLAST database based on downloaded sequences? \n\n Please respond yes or no. ")
		if yes_no(blast_q):
			if ff == ".nuc.fa":
				print(" Preparing a nucleotide database... ") 
				database= "makeblastdb -in " + dow_file_name + ff +" -dbtype nucl -out " + dow_file_name
				subprocess.call(database,shell=True)
				if os.path.isfile(dow_file_name + ".nhr"):
					print(" Nucleotide database has been successfully made! \n Nucleotide sequence headers are in " + dow_file_name + ".nhr \n Nucleotide indexes are in " + dow_file_name + ".nin \n Compressed nucleotide sequences are in " + dow_file_name + ".nsq \n")
			else:
				print(" Preparing a protein database... ") 
				database= "makeblastdb -in " + dow_file_name + ff +" -dbtype prot -out " + dow_file_name
				subprocess.call(database,shell=True)
				if os.path.isfile(dow_file_name + ".phr"):
					print(" Protein database has been successfully made! \n Protein sequence headers are in " + dow_file_name + ".phr \n Protein indexes are in " + dow_file_name + ".pin \n Compressed protein sequences are in " + dow_file_name + ".psq \n")


search(*list(details.values()))
blastdb(*list(details.values()))
