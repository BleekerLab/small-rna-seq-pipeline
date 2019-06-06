#! /usr/bin/python3

"""
This script is to describe sRNA sequences e.g. mature miRNA based on their 5' and 3' nucleotides.
"""

import sys
import pandas as pd
from Bio.Seq import Seq
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
print('sys.argv is', sys.argv)

def main():
	script = sys.argv[0]
	ph = sys.argv[1]
	filename = sys.argv[2]
	
	ph = ph
  
  # Data is selected, two dataframes are formed one with all the MajorRNAs, and one with MajorRNAs which are miRNA
	data = pd.read_csv(ph, sep= '\t')
	data = data[['Name', 'MajorRNA', 'DicerCall', 'MIRNA']]
	df_miRNA = data[(data.MIRNA == 'Y')] 
	
	#Lists are formed from all the MajorRNAs and from all the MIRNA MajorRNAs
	data_all = list(data.MajorRNA)
	data_miRNA = list(data[(data.MIRNA == 'Y')].MajorRNA)


	# In this step all MajorRNAs are filtered on DicerCall=20-24, then these MajorRNAs are written in a list

	data_DC24 = list(data[(data.DicerCall == '24')].MajorRNA)
	data_DC23 = list(data[(data.DicerCall == '23')].MajorRNA)
	data_DC22 = list(data[(data.DicerCall == '22')].MajorRNA)
	data_DC21 = list(data[(data.DicerCall == '21')].MajorRNA)
	data_DC20 = list(data[(data.DicerCall == '20')].MajorRNA)
	#allDCrest =

	# In this step all MIRNA MajorRNAs are filtered on DicerCall=20-24, then these MajorRNAs are written in a list
	data_MIRDC24 = list(df_miRNA[(df_miRNA.DicerCall == '24')].MajorRNA)
	data_MIRDC23 = list(df_miRNA[(df_miRNA.DicerCall == '23')].MajorRNA)
	data_MIRDC22 = list(df_miRNA[(df_miRNA.DicerCall == '22')].MajorRNA)
	data_MIRDC21 = list(df_miRNA[(df_miRNA.DicerCall == '21')].MajorRNA)
	data_MIRDC20 = list(df_miRNA[(df_miRNA.DicerCall == '20')].MajorRNA)
	# allMIRDCrest = []

	# all.. is for all MajorRNAs selected on length of the RNA, allMIR.. is for all miRNAs selected on length of the RNA
	data_size20 = [] 
	data_size21 = []
	data_size22 = []
	data_size23 = []
	data_size24 = []
	data_sizerest = []
	data_sizeMIR20 = []
	data_sizeMIR21 = []
	data_sizeMIR22 = []
	data_sizeMIR23 = []
	data_sizeMIR24 = []
	data_sizeMIRrest = []


	for s in data_all:
		if len(s) == 20:
			data_size20.append(s)
		elif len(s) == 21:
			data_size21.append(s)
		elif len(s) == 22:
			data_size22.append(s)
		elif len(s) == 23:
			data_size23.append(s)
		elif len(s) == 24:
			data_size24.append(s)
		else:
			data_sizerest.append(s)

	for s in data_miRNA:
		if len(s) == 20:
			data_sizeMIR20.append(s)
		elif len(s) == 21:
			data_sizeMIR21.append(s)
		elif len(s) == 22:
			data_sizeMIR22.append(s)
		elif len(s) == 23:
			data_sizeMIR23.append(s)
		elif len(s) == 24:
			data_sizeMIR24.append(s)
		else:
			data_sizeMIRrest.append(s)

	# The start and end nts are counted and written in a list [[5' A-U-C-G], [3' A-U-C-G]]         
	nt_all = []
	nt_miRNA = []
	nt_DC24 = []
	nt_DC23 = []
	nt_DC22 = []
	nt_DC21 = []
	nt_DC20 = []
	nt_size24 = []
	nt_size23 = []
	nt_size22 = []
	nt_size21 = []
	nt_size20 = []
	nt_MIRDC24 = []
	nt_MIRDC23 = []
	nt_MIRDC22 = []
	nt_MIRDC21 = []
	nt_MIRDC20 = []
	nt_sizeMIR24 = []
	nt_sizeMIR23 = []
	nt_sizeMIR22 = []
	nt_sizeMIR21 = []
	nt_sizeMIR20 = []

	dataset = ['data_all', 'data_miRNA', 
			   'data_DC24', 'data_DC23', 'data_DC22', 'data_DC21', 'data_DC20',
			   'data_MIRDC24', 'data_MIRDC23', 'data_MIRDC22', 'data_MIRDC21', 'data_MIRDC20', 
			   'data_size24', 'data_size23', 'data_size22', 'data_size21', 'data_size20',
			   'data_sizeMIR24', 'data_sizeMIR23', 'data_sizeMIR22', 'data_sizeMIR21', 'data_sizeMIR20']

	for d in dataset:
		count5A = 0
		count5U = 0
		count5C = 0
		count5G = 0
		name = d.replace('data', 'nt')
		for s in vars()[d]:
			if s[0] == "A":
				count5A += 1
			elif s[0] == "U":
				count5U += 1
			elif s[0] == "C":
				count5C += 1
			elif s[0] == "G":
				count5G += 1
			else:
				print('This not a nucleotide')
				print(s[0])
		vars()[name].append([count5A, count5U, count5C, count5G])
	

		count3A = 0
		count3U = 0
		count3C = 0
		count3G = 0
		name = d.replace('data', 'nt')
		for s in vars()[d]:
			if s[-1] == "A":
				count3A += 1
			elif s[-1] == "U":
				count3U += 1
			elif s[-1] == "C":
				count3C += 1
			elif s[-1] == "G":
				count3G += 1
			else:
				print('This not a nucleotide')
		vars()[name].append([count3A, count3U, count3C, count3G])
		print([name])
		print(vars()[name])

	figure = ['nt_all','nt_miRNA','nt_DC24','nt_DC23','nt_DC22','nt_DC21','nt_DC20', 'nt_size24', 
			  'nt_size23','nt_size22','nt_size21','nt_size20','nt_MIRDC24','nt_MIRDC23','nt_MIRDC22',
			  'nt_MIRDC21','nt_MIRDC20','nt_sizeMIR24','nt_sizeMIR23','nt_sizeMIR22','nt_sizeMIR21','nt_sizeMIR20']


	for f in figure:
		a = []
		for d in vars()[f]:
			a.append(d)
		#print(a)    
		#print('morry')
		a1 = a[0]
		#print(a1)
		a2 = a[1]
		#print(a2)
		df = pd.DataFrame(
			{'5': [a1[0],a1[1],a1[2],a1[3]],
			 '3': [a2[0],a2[1],a2[2],a2[3]]},
			index= ['A', 'U', 'C', 'G']
		)
		#print(df)
		pes = df.plot(y= ['3', '5'],kind = 'pie', subplots=True, title = f)
		plt.savefig(f+'.png', bbox_inches='tight')
		#plt.show(all(pes))
	
	"""pdf = matplotlib.backends.backend_pdf.PdfPages("output.pdf")
	for fig in xrange(1, figure().number): ## will open an empty extra figure :(
		pdf.savefig( fig )
	pdf.close()"""
	
	
if __name__ == '__main__':
   main()
