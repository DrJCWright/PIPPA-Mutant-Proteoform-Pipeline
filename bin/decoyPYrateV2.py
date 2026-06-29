#!/software/bin/python3.2

# Decoy Database Generator - v2.0
# James Wright 2022


#used to get cmd line arguments 
import argparse
#used to shuffle peptides
import random
#used to rename/delete tmp file
import os
#used to extract none specifically cleaved peptides - Removed as faster method implemented
#from itertools import combinations


#Read command line arguments and create help documentation using argparse
parser = argparse.ArgumentParser(
	description='''Create decoy protein sequences. Each protein is reversed and the cleavage sites switched with preceding amino acid. 
		Peptides are checked for existence in target sequences if found the tool will attempt to shuffle them. Version 2.0 allows none specific cleavage, amino acid substitutions and peptide FASTA formatting.

		Recommended settings for tryptic decoy database:
		python decoyPYrateV2.py -c KR -a P -p c -ml 6 -mx 100 -n 100 -u 10 target.fasta

		Recommended settings for none specific cleavage of peptides length 9:
		python decoyPYrateV2.py -c X -a X -p X -ml 9 -mx 9 -n 100 -u 10 -m target.fasta

		James.Wright@icr.ac.uk 2022''')

parser.add_argument('fasta', metavar='*.fasta|*.fa', help='FASTA file of target proteins sequences for which to create decoys')
parser.add_argument('--cleavage_sites', '-c', dest='csites', default='KR', help='A list of amino acids at which to cleave during digestion. Use X for unspecified cleavage. Default = KR')
parser.add_argument('--anti_cleavage_sites', '-a', dest='noc', default='P', help='A list of amino acids at which not to cleave if following cleavage site ie. Proline. Default = P')
parser.add_argument('--cleavage_position', '-p', dest='cpos', default='c', choices=['c', 'n', 'x'], help='Set cleavage to be c or n terminal of specified cleavage sites. Use x for none specific. Default = c')
parser.add_argument('--min_peptide_length', '-ml', dest='minlen', default=5, type=int, help='Set minimum length of peptides to compare between target and decoy. Default = 5')
parser.add_argument('--max_peptide_length', '-mx', dest='maxlen', default=100, type=int, help='Set maximum length of peptides to compare between target and decoy. Default = 100')
parser.add_argument('--do_not_switch', '-s', dest='noswitch', default=False, action='store_true', help='Turn OFF switching of cleavage site with preceding amino acid. Default=false')
parser.add_argument('--do_not_shuffle', '-x', dest='noshuf', default=False, action='store_true', help='Turn OFF shuffling of decoy peptides that are in the target database. Default=false')
parser.add_argument('--max_iterations', '-n', dest='maxit', default=100, type=int, help='Set maximum number of times to shuffle a peptide to make it non-target. Default=100')
parser.add_argument('--do_not_substitute', '-r', dest='nosub', default=False, action='store_true', help='Turn OFF subsituting random amino acids when shuffling fails. Default=false')
parser.add_argument('--max_substitutions', '-u', dest='maxsub', default=100, type=int, help='Set maximum number of times to substitute a peptide amino acid to make it non-target. Default=100')
parser.add_argument('--decoy_prefix', '-d', dest='dprefix', default='DECOY_', help='Set accesion prefix for decoy proteins in output. Default=DECOY_')
parser.add_argument('--output_targets', '-ot', dest='tinc', default=False, action='store_true', help='Include target peptides in output fasta. Default=false')
parser.add_argument('--output_fasta', '-o', dest='dout', default='decoy.fa', help='Set file to write decoy proteins to. Default=decoy.fa')
parser.add_argument('--temp_file', '-t', dest='tout', default='tmp.fa', help='Set temporary file to write decoys prior to shuffling. Default=tmp.fa')
parser.add_argument('--no_isobaric', '-i', dest='iso', default=False, action='store_true', help='Do not make decoy peptides isobaric. Default=false')
parser.add_argument('--peptide_fasta', '-pf', dest='pfasta', default=False, action='store_true', help='Output decoys a individual peptides sequences. Default=false')
parser.add_argument('--memory_save', '-m', dest='mem', default=False, action='store_true', help='Slower but uses less memory (does not store decoy peptide list). Default=false')
args = parser.parse_args()

##################################################################################################

def digest(protein, sites, pos, no, min, max):
	"""Return a list of cleaved peptides with minimum length in protein sequence.
		protein = sequence
		sites = string of amino acid cleavage sites
		pos = n or c for n-terminal or c-terminal cleavage
		no = amino acids following site that would prevent cleavage ie proline
		min = minimum length of peptides returned
		max = maximum length of peptides returned"""

	#for each possible cleavage site insert a comma before or after depending on pos
	for s in sites:
		r = s + ','
		if pos == 'n':
			r = ',' + s
		protein = protein.replace(s, r)
	
	#for each possible cleavage and all none cleavage remove comma
	for s in sites:
		for n in no:
			a = s + ',' + n
			if pos == 'n':
				a = ',' + s + n
			r = s + n
			protein = protein.replace(a, r)
	
	#filter peptides into list by minimum size
	return list(filter(lambda x: len(x)>=min and len(x)<=max, (protein.split(','))))

##################################################################################################

def nsdigest(protein, min, max):
	"""Return a list of none specific cleaved peptides with minimum length and maximum length in protein sequence.
		protein = sequence
		min = minimum length of peptides returned
		max = maximum length of peptides returned"""

	#Use combinations function to efficiently generate all substr peptides of lengths between min and max (Great but very slow for large proteins like TTN)
	#return [protein[x:y] for x, y in combinations(range(len(protein) + 1), r = 2) if (y-x) >= min and  (y-x) <= max ] 

	#Much faster to just use sliding window
	nspeps = list()
	l = len(protein)
	for i in range(min, (max+1)):
		for x in range(l+1):
			y = x + i
			if y <= l:
				nspeps.append(protein[x:y])
	return nspeps



##################################################################################################

def revswitch (protein, sites, noswitch):
	"""Return a reversed protein sequence with cleavage residues switched with preceding residue"""
	#reverse protein sequence with a reverse splice convert to list
	revseq = list(protein[::-1])
	
	if noswitch == False and 'X' not in sites:
		#loop sequence list
		for i, c in enumerate(revseq):
			#if value is in sites switch with previous amino acid
			if c in sites:
				aa = revseq[i-1]
				revseq[i-1] == revseq[i]
				revseq[i] == aa
	#return reversed with/without switched proteins as string
	return ''.join(revseq)	


##################################################################################################
	
def shuffle(peptide, pos):
	"""shuffle peptide without moving c-terminal amino acid cleavage site"""

	n = ''
	c = ''

	#convert peptide to list 
	l = list(peptide)

	if pos == 'n':
		#extract n terminal aa
		n = l.pop(0)

	elif pos == 'c':
		#extract c terminal aa
		c = l.pop()
	
	random.shuffle(l)
	#return new peptide
	return n + ''.join(l) + c 
	

##################################################################################################	

def substitute(peptide, pos):
	"""Substitute a random amino acid in a peptide, fix n or c termini"""

	n = ''
	c = ''

	AminoAcids = ["A","C","D","E","F","G","H","K","L","M","N","P","Q","R","S","T","V","W","Y"]

	#convert peptide to list 
	l = list(peptide)

	if pos == 'n':
		#extract n terminal aa
		n = l.pop(0)

	elif pos == 'c':
		#extract c terminal aa
		c = l.pop()
	
	l[random.choice(range(len(l)))] = random.choice(AminoAcids)

	#return new peptide
	return n + ''.join(l) + c 


###################			MAIN 		#####################
#Create empty sets to add all target and decoy peptides
upeps = set()	
dpeps = set()
	
#Counter for number of decoy sequences
dcount = 0;	

#empty protein sequence
seq = ''	
	
	
#open temporary decoy FASTA file
outfa = open(args.tout, 'w')	
	
#Open FASTA file using first cmd line argument
fasta = open(args.fasta, 'r')
#loop each line in the file
for line in fasta:
	#if this line starts with ">" then process sequence if not empty
	if line[0] == '>':
		if seq != '':
		
			#make sequence isobaric (check args for switch off)
			if args.iso == False:
				seq = seq.replace('I', 'L')
		
			if 'X' in args.csites:
				#nonspecific digest sequence add peptides to set
				upeps.update( nsdigest(seq, args.minlen, args.maxlen) )
			else:
				#digest sequence add peptides to set
				upeps.update( digest(seq, args.csites, args.cpos, args.noc, args.minlen, args.maxlen) )
				
			#reverse and switch protein sequence
			decoyseq = revswitch(seq, args.csites, args.noswitch)
			
			#do not store decoy peptide set in reduced memory mode
			if args.mem == False:

				if 'X' in args.csites:
					#nonspecific digest sequence add peptides to set
					dpeps.update( nsdigest(seq, args.minlen, args.maxlen) )
				else:
					#update decoy peptide set
					dpeps.update( digest(decoyseq, args.csites, args.cpos, args.noc, args.minlen, args.maxlen) )
			
			#write decoy protein accession and sequence to temporary file
			dcount += 1
			outfa.write('>' + args.dprefix + '_' + str(dcount) + '\n')
			outfa.write(decoyseq + '\n')
			
		seq = '';
	
	#if not accession line then append aa sequence (with no newline or white space) to seq string
	else:
		seq+=line.rstrip()
		
#Close files
fasta.close()

#Final Sequence Processing
if seq != '':
		
	#make sequence isobaric (check args for switch off)
	if args.iso == False:
		seq = seq.replace('I', 'L')
		
	if 'X' in args.csites:
		#nonspecific digest sequence add peptides to set
		upeps.update( nsdigest(seq, args.minlen, args.maxlen) )
	else:
		#digest sequence add peptides to set
		upeps.update( digest(seq, args.csites, args.cpos, args.noc, args.minlen, args.maxlen) )
		
	#reverse and switch protein sequence
	decoyseq = revswitch(seq, args.csites, args.noswitch)
	
	#do not store decoy peptide set in reduced memory mode
	if args.mem == False:

		if 'X' in args.csites:
			#nonspecific digest sequence add peptides to set
			dpeps.update( nsdigest(seq, args.minlen, args.maxlen) )
		else:
			#update decoy peptide set
			dpeps.update( digest(decoyseq, args.csites, args.cpos, args.noc, args.minlen, args.maxlen) )
	
	#write decoy protein accession and sequence to temporary file
	dcount += 1
	outfa.write('>' + args.dprefix + '_' + str(dcount) + '\n')
	outfa.write(decoyseq + '\n')

outfa.close()



#Summarise the numbers of target and decoy peptides and their intersection
nonDecoys = set()
print ("proteins:" + str(dcount))
print ("target peptides:" + str(len(upeps)))

#Reloop decoy file in reduced memory mode to store only intersecting decoys 
if args.mem:
	#open temp decoys
	with open(args.tout, "rt") as fin:
		for line in fin:

			#if line is not accession 
			if line[0] != '>':

				if 'X' in args.csites:
					#nonspecific digest protein
					for p in nsdigest(line.rstrip(), args.minlen, args.maxlen):
						#check if in target peptides if true then add to nonDecoys
						if p in upeps:
							nonDecoys.add(p)

				else:
					#digest protein
					for p in digest(line.rstrip(), args.csites, args.cpos, args.noc, args.minlen, args.maxlen):
						#check if in target peptides if true then add to nonDecoys
						if p in upeps:
							nonDecoys.add(p)
	fin.close()
	print ("decoy peptides: !Memory Saving Made!")
else:
	#can only report total number in normal memory mode
	print ("decoy peptides:" + str(len(dpeps)))
	#find intersecting peptides
	nonDecoys = upeps.intersection(dpeps)

#print ("intersection:" + str(nonDecoys))
print ("#intersection:" +  str(len(nonDecoys)))

#if there are decoy peptides that are in the target peptide set
if len(nonDecoys) > 0 or args.noshuf == False or args.pfasta == True:

	#create empty dictionary with bad decoys as keys
	#dAlternative = dict.fromkeys(nonDecoys, '')
	#Group bad decoys by length
	dAlternative = dict()
	for p in nonDecoys:
		plen = len(p)
		if plen not in dAlternative:
			#add peptide length to dictionary
			dAlternative[plen] = dict()
		dAlternative[plen][p] = ''

	noAlternative = set()
	
	#loop bad decoys / dictionary keys
	for plength in sorted(dAlternative.keys(), reverse=True):
		for dPep in dAlternative[plength]:
			subs = 0
			aPep = dPep

			while aPep in upeps and subs < args.maxit:

				shuffles = 0
				# shuffle until aPep is not in target set (maximum of 10 iterations)
				while aPep in upeps and shuffles < args.maxit:
				
					#increment iteration counter
					shuffles += 1
			
					#shuffle peptide
					aPep = shuffle(dPep, args.cpos)
					
					#check if shuffling has an effect if not end iterations
					if (aPep == dPep):
						shuffles=args.maxit

				if aPep in upeps:
					#Substitute random amino acid
					aPep = substitute(dPep, args.cpos)
					subs += 1
				
				#print (str(i) + '\t' + dPep + '\t' + aPep)
			
			#warn if peptide has no suitable alternative, add to removal list
			if aPep in upeps:
				noAlternative.add(dPep)
			else:
				#update dictionary with alternative shuffled peptide
				dAlternative[plength][dPep] = aPep

				#test sub peptides to see if this solution is vaild for them too
				#only if not already minimum length and none specific cleavage
				if plength > args.minlen and 'X' in args.csites:

					#generate all sub peptides of length minlen to plength-1	
					cmax = plength-1
					for i in range(args.minlen, cmax):
						# loop to get each substr of dPep of length i
						for x in range(plength - i + 1):
							subdPep = dPep[x:x + i]
							# if sub peptide is in target set, replace it with the corresponding alternative
							if subdPep in dAlternative[i]:
								subaPep = aPep[x:x + i]
								#check if alternative is in target set
								if subaPep not in upeps:
									del dAlternative[i][subdPep]

		
	
	print ( str(len(noAlternative)) + ' have no alternative peptide')
	#remove peptides with no alternative
	for p in noAlternative:
		#print (p + ' Warning: No Alternative found') 
		plen = len(p)
		del dAlternative[plen][p]
	
	#Free up memory by clearing large sets of peptides
	upeps.clear()
	dpeps.clear()
	#open second decoy file
	with open(args.dout, "wt") as fout:
		#open original decoy file
		with open(args.tout, "rt") as fin:
			#loop each line of original decoy fasta
			for line in fin:
				#if line is not accession replace peptides in dictionary with alternatives
				if line[0] != '>':
				
					opeps = set()

					if 'X' in args.csites:
						#nonspecific digest protein
						for p in nsdigest(line.rstrip(), args.minlen, args.maxlen):
								
							plen = len(p)

							#check if length in dictionary
							if plen not in dAlternative:
								dpeps.add(p)
								continue

							#if decoy peptide is in dictionary append with alternative
							if p in dAlternative[plen]:
								if args.pfasta == False:
									line = dAlternative[plen][p] + line
								dpeps.add(dAlternative[plen][p])
							else:
								dpeps.add(p)

					else:
						#digest protein
						for p in digest(line.rstrip(), args.csites, args.cpos, args.noc, args.minlen, args.maxlen):
								
							#if decoy peptide is in dictionary replace with alternative
							if p in dAlternative:
								plen = len(p)
								if args.pfasta == False:
									line = line.replace(p, dAlternative[plen][p])
								dpeps.add(dAlternative[plen][p])
							else:
								dpeps.add(p)
							
					#write decoy protein to output file
					if not args.pfasta:
						fout.write(line)


				else:
					#write accession line to output file
					if not args.pfasta:
						fout.write(line)		
		fin.close()

		if args.pfasta == True:
			for i, p in enumerate(dpeps):
				fout.write('>' + args.dprefix + '_' + str(i) + '\n')
				fout.write(p + '\n')

		#Append target proteins to decoy file if requested
		if args.tinc:
			#Open FASTA file using first cmd line argument
			fasta = open(args.fasta, 'r')
			#loop each line in the file
			for line in fasta:
				
				if line[0] != '>':
					#make sequence isobaric (check args for switch off)
					if args.iso == False:
						line = line.replace('I', 'L')

				fout.write(line)


	fout.close()
	
	#delete temporary file
	os.remove(args.tout)
else:
	os.rename(args.tout, args.dout)
	
print ("final decoy peptides:" + str(len(dpeps)))




#Possible function to calculate number of permutations
"""
import operator
from collections import Counter
from math import factorial
from functools import reduce
def npermutations(l):
    num = factorial(len(l))
    mults = Counter(l).values()
    den = reduce(operator.mul, (factorial(v) for v in mults), 1)
    return num / den"""	