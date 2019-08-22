#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
# Identifies kernels corresponding to spectra
#
# idX version 2019.08.10.02
#

import ujson
import time
import gzip
import sys
import datetime

# import the method that deals with spectrum file formats
from spectraX import load_spectra
# import the method for the output of results to a file
from reportX import report_ids

version = '2019.08.10.2'

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
	hmatched = 0
	# run through all lines in the kernel
	# each kernel is indexed by its line number, 'i'
	for i,l in enumerate(f):
		# print progress indicator
		if i != 0 and i % 2500 == 0:
			print('.',end='')
			sys.stdout.flush()
		if i != 0 and i % 50000 == 0:
			print(' %i' % (i))
			sys.stdout.flush()
		# convert the line into json
		js = ujson.loads(l)
		# check json for a parent mass value
		if 'pm' not in js:
			continue
		if 'u' not in js:
			print('\nERROR: kernel file does not have unique identifiers')
			print('\n\tline %i: check for "u" value in "%s"' % (i+1,fname))
			f.close()
			exit() 
		u = js['u']
		if js['u'] != js['h']:
			hmatched += 1
			continue
		pm = js['pm']
#		if not (js['seq'] == 'THEAQIQEMR'):
#			continue
		# test to see if parent mass a possible value, based on the spectra
		mv = int(0.5+pm/pt)
		if not (mv in sp_set or mv-1 in sp_set or mv+1 in sp_set):
			skipped += 1
			continue
		# add parent ion mass to dictionary
		mindex[u] = pm
		# add b ions to index
		bs = js['bs']
		if mv not in kindex:
			kindex[mv] = {}
		cindex = kindex[mv]
		for s in bs:
			sv = int(0.5+s/ft)
			if sv not in cindex:
				cindex[sv] = [u]
			else:
				cindex[sv].append(u)
		# add y ions to index
		ys = js['ys']
		for s in ys:
			sv = int(0.5+s/ft)
			if sv not in cindex:
				cindex[sv] = [u]
			else:
				cindex[sv].append(u)
	# finish up and return indexes
#	print('\n\t      skipped = %i' % (skipped))
#	print('\t     identical = %i' % (hmatched))
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
	for c,s in enumerate(_sp):
		if c != 0 and c % 50 == 0:
			print('.',end='')
			sys.stdout.flush()
		if c != 0 and c % 1000 == 0:
			print(' %i' % (c))
			sys.stdout.flush()
		ident = []
		ims = []
		rm = float(s['pm'])
		pm = int(0.5+rm/pt)
		# check _ki for appropriate fragment ion matches
		# and record a list of kernel identifiers for each match
		# and the corresponding fragment ion intensities
		idx = None
		idi = None
		for d in dvals:
			if pm+d not in _ki:
				continue
			ci = _ki[pm+d]
			for b,m in enumerate(s['sms']):
				if m in ci:
					idx = ci[m]
					idi = s['ims'][b]
				else:
					continue
				ident.append(idx)
				ims.append(idi)
		total = float(sum(s['ims']))/3
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
			for a in set(i):
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
		if mn > 4:
			ks = []
			vis = []
			max_i = max(iv)
			for b,m in enumerate(mv):
				if iv[b] < max_i:
					continue
				ks.append(m)
			id_values = {'sn':z,'peaks':mn,'kernels':ks,'ri': max_i/total,'pm': s['pm'],
					'pz': s['pz'],'sc': s['sc'],'ions' : len(s['sms'])/3}
			ids.append(id_values)
		z += 1
	print('')
	return ids

#
# Coordinate the identification process, print job stats and progress 
#

def main():
	if len(sys.argv) < 4:
		print('usage:\n\t>python3 idX.py SPECTRA_FILE KERNEL_FILE OUTPUT_FILE (high|medium|low*)')
		exit()
	start = time.time()
	# record relavent parameters
	param = {}
	#fragment tolerance in millidaltons
	param['fragment mass tolerance'] = float(400)
	try:
		if sys.argv[4] == 'high':
			param['fragment mass tolerance'] = float(20)
		elif sys.argv[4] == 'low':
			param['fragment mass tolerance'] = float(400)
		elif sys.argv[4] == 'medium':
			param['fragment mass tolerance'] = float(100)
		else:
			print('ERROR: argument 4 must be high or low, not "%s"'% (sys.argv[4]))
			exit()
	except:
		pass
	param['maximum spectra'] = -1
	try:
		param['maximum spectra'] = int(sys.argv[5])
	except:
		pass
	# parent tolerance in ppm
	param['parent mass tolerance'] = float(20)
	spectra = []
	# report files named on command line
	print('\nstart ...\nidX parameters')
	print('\t   max spectra: %i mDa' % (param['maximum spectra']))
	print('\t  fragment tol: %i mDa' % (param['fragment mass tolerance']))
	print('\t    parent tol: %i ppm' % (param['parent mass tolerance']))
	print('\t spectrum file: %s' % (sys.argv[1]))
	print('\t   kernel file: %s' % (sys.argv[2]))
	print('\t   output file: %s' % (sys.argv[3]))
	print('\t       version: %s' % (version))
	print('\t      run time: %s' % (str(datetime.datetime.now())))
	param['spectrum file'] = sys.argv[1]
	print('load & index spectra')
	# read the spectrum file and perform all necessary spectrum conditioning
	spectra = load_spectra(param['spectrum file'],param)
	if param['maximum spectra'] != -1:
		spectra = spectra[0:param['maximum spectra']]
	if len(spectra) == 0:
		print('exiting: 0 spectra found')
		print('done')
		exit()
	delta = time.time()-start
	start = time.time()
	print('\n\t   spectra = %i' % (len(spectra)))
	print('\tspectra &Delta;T = %.1f s' % (delta))
	param['kernel file'] = sys.argv[2]
	print('load & index kernel')
	# read the kernel file and create an index of peptide fragmentation patterns
	(ki,mi) = index_kernel(param,spectra)
	delta = time.time()-start
	start = time.time()
	print('\n\t   kernels = %i' % (len(ki)))
	print('\t &Delta;T = %.1f s' % (delta))
	print('perform ids')
	# generate a list of identifications for the spectra using the kernel index
	ids = create_ids(ki,mi,spectra,param)
	# free memory associated with indexes and spectra
	delta = time.time()-start
	start = time.time()
	print('\tid &Delta;T = %.3f s' % (delta))
	if len(spectra) > 0:
		print('\t   &delta;T = %.0f microseconds' % (1.0e06*delta/len(spectra)))
	else:
		pass
	# simple reporting of the kernels assigned to spectra
	print('release memory')
	ki = None
	spectra = None
	print('\tdone')
	param['output file'] = sys.argv[3]
	print('create report')
	report_ids(ids,param)
	print('... done')

if __name__== "__main__":
	main()


