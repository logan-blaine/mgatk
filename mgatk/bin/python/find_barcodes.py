#!/usr/bin/python

import sys
import re
import pysam
import os
from collections import Counter
from contextlib import contextmanager

bamfile = sys.argv[1]
outfolder = sys.argv[2]
barcodeTag = sys.argv[3]
min_barcodes = int(sys.argv[4])
mtchr = sys.argv[5]
quant_file = sys.argv[6]
passing_file = sys.argv[7]
whitelist_file = sys.argv[8]

base=os.path.basename(bamfile)
basename=os.path.splitext(base)[0]

def getBarcode(read):
	'''
	Parse out the barcode per-read
	'''
	try:
		return read.get_tag(barcodeTag)
	except KeyError:
		return("NA")
	# for tg in intags:
	# 	if(barcodeTag == tg[0]):
	# 		return(tg[1])
	# return("NA")

def quantifyBarcodes(mtchr):
	'''
	Make a giant dictionary of observed barcodes at the mitochondrial chr
	'''
	barcodes_all = dict()
	bam = pysam.AlignmentFile(bamfile,'rb')
	Itr = bam.fetch(str(mtchr),multiple_iterators=False)
	
	for read in Itr:
		read_barcode = getBarcode(read)
		barcodes_all[read_barcode] = barcodes_all.get(read_barcode, 0) + 1
	bam.close()
	return {x : n for x, n in barcodes_all.items() if n >= min_barcodes and x != "NA"}

# Quant barcodes and write it out

barcodes = quantifyBarcodes(mtchr)
# bc = list(barcodes.keys())	
with open(quant_file, "w") as quant_file_o:
	for k, v in barcodes.items():
		quant_file_o.write(k +"\t"+ str(v)+"\n")

# intersect discovered barcodes with whitelist, if present
if whitelist_file:
	with open(whitelist_file) as f:
		whitelist={l.strip() for l in f}
	barcodes = {x : n for x, n in barcodes.items() if x in whitelist}

# write final passing barcodes 
with open(passing_file, "w") as passing_file_o:
	for k, v in barcodes.items():
		passing_file_o.write(k +"\n")
