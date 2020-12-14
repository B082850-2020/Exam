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
details["gene"]   = input(" Which gene do you want to search for? (leave empty if no specific gene is wanted) ")
details["db"]    = input(" Do you want to make a database for protein or nucleotide? Please answer p or n")
details["complete"] = input (" Do you want to include partial/incomplete result? Please answer yes or no ")

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
	yes = set(['protein','p'])		# both yes or y will return True
	no = set(['nucleotide','n']) 		# both no or n will return False
	choice = answer.lower()		# set all types of answer to lower case
	# loop forever until something is returned
	while True:		
		if choice in yes:		# both yes or y will return True
			return True
		elif choice in no:		# both no or n will return False
			return False
		else:					# error trap, all other input will cause the function to ask "yes or no" over and over until a desired answer is received
			print ("\n Please respond with 'p' or 'n'. \n")
			choice = input(" protein or nucleotide?").lower()
			
## function to check the validility of user input
def check_details(organism,gene,db,complete):
	if len(organism)==0 or str.isspace(organism):
		print("\n Sorry, no organism was chosen, please try again..\n") 
	elif len(gene)==0 or str.isspace(gene):
		print("\n Ok, no specific gene will be searched for. \n")
		if protein_nt(db): 
			if yes_no(complete):
				print("\n Search will be done on:\n\tOrganism:",organism,"\n\tDatabase: protein","\n\tInclude partial/incomplete: Yes")
			else:
				print("\n Search will be done on:\n\tOrganism:",organism,"\n\tDatabase: protein","\n\tInclude partial/incomplete: No")
		else:
			if yes_no(complete):
				print("\n Search will be done on:\n\tOrganism:",organism,"\n\tDatabase: nucleotide","\n\tInclude partial/incomplete: Yes")
			else:
				print("\n Search will be done on:\n\tOrganism:",organism,"\n\tDatabase: protein","\n\tInclude partial/incomplete: No")
	else:
		if protein_nt(db): 
			if yes_no(complete):
				print("\n Search will be done on:\n\tOrganism:",organism,"\n\tGene:",gene,"\n\tDatabase: protein","\n\tInclude partial/incomplete: Yes")
			else:
				print("\n Search will be done on:\n\tOrganism:",organism,"\n\tGene:",gene,"\n\tDatabase: protein","\n\tInclude partial/incomplete: No")
		else:
			if yes_no(complete):
				print("\n Search will be done on:\n\tOrganism:",organism,"\n\tGene:",gene,"\n\tDatabase: nucleotide","\n\tInclude partial/incomplete: Yes")
			else:
				print("\n Search will be done on:\n\tOrganism:",organism,"\n\tGene:",gene,"\n\tDatabase: protein","\n\tInclude partial/incomplete: No")
	
	
check_details(*list(details.values()))
