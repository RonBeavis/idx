#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
# Identifies kernels corresponding to spectra
#
# idX version 2019.08.10.01
#

import ujson
#import json
import time
import re
import gzip
import sys
import itertools
import copy
# import the method that deals with spectrum file formats
from spectraX import load_spectra
# import the method for the output of results to a file
from reportX import report_ids

#
# load information from a kernel into dictionaries and sets
#

def index_kernel(_p,_s):
	ft = _p['fragment mass tolerance']
	# permissive parent tolerance in mDa
	pt = 70
	sp_set = set()
	fname = _p['kernel file']
	# create a set of spectra parent masses
	for s in _s:
		sp_set.add(int(0.5 + s['pm']/pt))
	kindex = {}
	mindex = dict()
	# open a handle for the kernel file (may be gzip'd)
	if fname.find('.gz') == len(fname) - 3:
		f = gzip.open(fname,'rt', encoding='utf-8')
	else:
		f = open(fname,'r', encoding='utf-8')
	skipped = 0
	# run through all lines in the kernel
	# each kernel is indexed by its line number, 'i'
	for i,l in enumerate(f):
		# print progress indicator
		if i % 10000 == 0:
			print('.',end='')
			sys.stdout.flush()
		# convert the line into json
		js = ujson.loads(l)
		# check json for a parent mass value
		if 'pm' not in js:
			continue
		pm = js['pm']
		# test to see if parent mass a possible value, based on the spectra
		mv = int(0.5+pm/pt)
		if not (mv in sp_set or mv-1 in sp_set or mv+1 in sp_set):
			skipped += 1
			continue
		# add parent ion mass to dictionary
		mindex[i] = pm
		# add b ions to index
		bs = js['bs']
		if mv not in kindex:
			kindex[mv] = {}
		cindex = kindex[mv]
		for s in bs:
			sv = int(0.5+s/ft)
			if sv not in cindex:
				cindex[sv] = {i}
			else:
				cindex[sv].add(i)
		# add y ions to index
		ys = js['ys']
		for s in ys:
			sv = int(0.5+s/ft)
			if sv not in cindex:
				cindex[sv] = {i}
			else:
				cindex[sv].add(i)
	# finish up and return indexes
#	print('\n\tskipped = %i' % (skipped))
	print('')
	f.close()
	return (kindex,mindex)
#
# generate identifications of spectra, based on the kernel indexes
#

def create_ids(_ki,_mi,_sp,_p):
	# initialize list of identifications
	ids = []
	z = 1
	pt = 70
	ppm = _p['parent mass tolerance']/1.0e06
	dvals = (-1,0,1)
	count = 0
	# iterate through spectra
	for s in _sp:
		ident = []
		ims = []
		rm = float(s['pm'])
		pm = int(0.5+rm/pt)
		# check _ki for appropriate fragment ion matches
		# and record a list of kernel identifiers for each match
		# and the corresponding fragment ion intensities
		for d in dvals:
			if pm+d not in _ki:
				continue
			ci = _ki[pm+d]
			for b,m in enumerate(s['sms']):
				v = m
				if v in ci:
					ident.append(ci[v])
					ims.append(s['ims'][b])
				elif v-1 in ci:
					ident.append(ci[v-1])
					ims.append(s['ims'][b])
				elif v+1 in ci:
					ident.append(ci[v+1])
					ims.append(s['ims'][b])
		total = float(sum(s['ims']))
		ans = {}
		aint = {}
		aok = set()
		# check kernel identifiers for the exact parent ion mass tolerance
		# and populate aok with identifiers that pass the test
		for b,i in enumerate(ident):
			for a in i:
				if a not in aok and abs(rm-_mi[a]) < ppm*rm:
					aok.add(a)
		# perform an ion count and intensity sum
		# for each kernel identifier
		for b,i in enumerate(ident):
			for a in i:
				if a not in aok:
					ans[a] = 0
					aint[a] = 0
				elif a in ans:
					ans[a] += 1
					aint[a] += ims[b]
				else:
					ans[a] = 1
					aint[a] = ims[b]
		mn = 0
		mv = []
		iv = []
		# order the ion and intensity statistics
		for j in ans:
			if ans[j] > mn:
				mn = ans[j]
				mv = [j]
				iv = [aint[j]]
			elif ans[j] == mn:
				mv.append(j)
				iv.append(aint[j])
		# select identifications based on their ion count and summed intensity
		if mn > 7:
			ks = []
			vis = []
			max_i = max(iv)
			for b,m in enumerate(mv):
				if iv[b] < max_i:
					continue
				ks.append(m)
			if max_i/total >= 0.20:
				ids.append({'sn':z,'peaks':mn,'kernels':ks,'ri': max_i/total,'pm': s['pm'],'pz': s['pz'],'sc': s['sc']})
		z += 1
	return ids

#
# Coordinate the identification process, print job stats and progress 
#

def main():
	if len(sys.argv) != 4:
		print('usage:\n\t>python3 idX.py SPECTRA_FILE KERNEL_FILE OUTPUT_FILE')
		exit()
	start = time.time()
	# record relavent parameters
	param = {}
	#fragment tolerance in millidaltons
	param['fragment mass tolerance'] = float(20)
	# parent tolerance in ppm
	param['parent mass tolerance'] = float(20)
	spectra = []
	# report files named on command line
	print(' spectrum file: %s' % (sys.argv[1]))
	print('   kernel file: %s' % (sys.argv[2]))
	print('   output file: %s' % (sys.argv[3]))

	param['spectrum file'] = sys.argv[1]
	print('load & index spectra')
	print('\t',end='')
	# read the spectrum file and perform all necessary spectrum conditioning
	spectra = load_spectra(param['spectrum file'],param)
	delta = time.time()-start
	start = time.time()
	print('\n\tspectra = %i' % (len(spectra)))
	print('\tspectra loading = %.1f s' % (delta))
	param['kernel file'] = sys.argv[2]
	print('load & index kernel')
	print('\t',end='')
	# read the kernel file and create an index of peptide fragmentation patterns
	(ki,mi) = index_kernel(param,spectra)
	delta = time.time()-start
	start = time.time()
	print('\tkernels = %i' % (len(ki)))
	print('\tkernel indexing = %.1f s' % (delta))
	print('perform ids')
	# generate a list of identifications for the spectra using the kernel index
	ids = create_ids(ki,mi,spectra,param)
	# free memory associated with indexes and spectra
	delta = time.time()-start
	start = time.time()
	print('\tid process = %.3f s' % (delta))
	if len(spectra) > 0:
		print('\ttime per spectrum = %.1f μs' % (1.0e06*delta/len(spectra)))
	else:
		pass
	# simple reporting of the kernels assigned to spectra
	print('release memory')
	ki = None
	spectra = None
	param['output file'] = sys.argv[3]
	print('create report')
	report_ids(ids,param)
	print('done')

if __name__== "__main__":
	main()


