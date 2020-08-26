#!/usr/bin/env python3
import sys
import argparse
import concurrent.futures
import re

def read_fasta(fasta):
	with open(fasta, 'r') as inf:
		raw_fasta=inf.readlines()
	if not raw_fasta[0].startswith('>'):
		print('The specified alignment file is not a valid Fasta file.')
		sys.exit()
	for i in range(len(raw_fasta)):
		raw_fasta[i]=raw_fasta[i].strip()
	alignment={}
	for i in range(len(raw_fasta)):
		if '>' in raw_fasta[i]:
			h=i
			raw_fasta[h]=raw_fasta[h].replace('>', '')
			raw_fasta[h]=raw_fasta[h].replace(' ', '_')
			alignment[raw_fasta[h]]=[]
			continue
		else:
			alignment[raw_fasta[h]]+=[raw_fasta[i].upper()]
			continue
	for i in alignment.keys():
		alignment[i]=''.join(alignment[i])
	return alignment

def repeat_finder(key, start, length, sequence):
	print('Start '+key)
	repeats={}
	positions=[]
	sequence=sequence[start:-1]
	sequence=alignment[key].replace('-', '')
	osequence=sequence
	l=len(sequence)//2
	n=0
	while l<len(sequence)-length:
		while l>length:
			cand=sequence[0:l]
			if sequence.count(cand)>1:
				if cand not in repeats.keys():
					n+=1
					repeats[cand]=[n,[]]
				sequence=sequence[l:-1]
				l=len(sequence)//2
				continue
			else:
				l-=1
				continue
		sequence=sequence[1:-1]
		l=len(sequence)//2
	for k in repeats.keys():
		repeats[k][1]=[m.start()+start+1 for m in re.finditer(k, osequence)]
		positions+=repeats[k][1]
	positions=list(set(positions))
	positions.sort()
	for i in range(len(positions)):
		for k in reversed(repeats.keys()):
			if positions[i] in repeats[k][1]:
				positions[i]=repeats[k][0]
	print('\n--- '+key+' ---\n')
	for k in repeats.keys():
		print(str(repeats[k][0])+'  '+k)
		print(repeats[k][1])
	for p in positions:
		print(u"\u001b[48;5;"+str(p)+'m '+str(p), end=' ')
	print(u'\u001b[0m\n')
	
parser=argparse.ArgumentParser(description='The script finds repeated motifs of seqences in multifasta file. The output is a sequnce name followed by motifs and their positions in the input alignment. The colored blocks  below show the order of repeats in the sequence (ignoring spaces and overlaps).')
parser.add_argument('-a', '--alignment', type=str, metavar='', required=True, help='Name of (or path to) the input Fasta file.')
parser.add_argument('-l', '--length', type=int, metavar='', default=10, required=False, help='Minimum length of motifs. Default=10.')
parser.add_argument('-s', '--start', type=int, metavar='', default=0, required=False, help='Analysis starting position. Default=0.')
arguments=parser.parse_args()

alignment=read_fasta(arguments.alignment)

with concurrent.futures.ProcessPoolExecutor() as executor:
	for k in alignment.keys():
		executor.submit(repeat_finder, k, arguments.start, arguments.length, alignment[k])
