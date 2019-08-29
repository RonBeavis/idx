#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
# Identifies kernels corresponding to spectra
#
# idX version 2019.08.10.02
#
from __future__ import print_function
from libcpp cimport bool as bool_t

import ujson
import time
import gzip
import sys
import datetime
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
	c13 = 1003
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
#		if not (js['seq'] == 'LPNLTHLNLSGNK'):
#			continue
		# test to see if parent mass a possible value, based on the spectra
		mv = int(0.5+pm/pt)
		cv = int(0.5+(pm+c13)/pt)
		if mv in sp_set or mv-1 in sp_set or mv+1 in sp_set:
			pass
		elif pm > 1500000 and (cv in sp_set or cv-1 in sp_set or cv+1 in sp_set):
			pass
		else:
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

