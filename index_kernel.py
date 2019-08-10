#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#

#
# Identifies kernels corresponding to spectra
#
# Reads kernels in either plain text or gzip'd text
#

#import ujson
import json
import time
import re
import gzip
import sys
import itertools
import copy
from load_spectra import load_spectra

def index_kernel(_f,_p,_s):
	ft = float(_p['fragment mass tolerance'])
	pt = float(_p['parent mass tolerance'])
	sp_set = set()
	for s in _s:
		sp_set.add(int(0.5 + s['pm']/pt))
	kindex = {}
	if _f.find('.gz') == len(_f) - 3:
		f = gzip.open(_f,'rt', encoding='utf-8')
	else:
		f = open(_f,'r', encoding='utf-8')
	skipped = 0
	for i,l in enumerate(f):
		if i % 10000 == 0:
			print('.',end='')
			sys.stdout.flush()
		js_master = json.loads(l)
		if 'pm' not in js_master:
			continue
		pm = js_master['pm']
		mv = int(0.5+pm/pt)
		if not (mv in sp_set or mv-1 in sp_set or mv+1 in sp_set):
			skipped += 1
			continue
		bs = js_master['bs']
		if mv not in kindex:
			kindex[mv] = {}
		cindex = kindex[mv]
		for s in bs:
			sv = int(0.5+s/ft)
			if sv not in cindex:
				cindex[sv] = {i}
			else:
				cindex[sv].add(i)
		ys = js_master['ys']
		for s in ys:
			sv = int(0.5+s/ft)
			if sv not in cindex:
				cindex[sv] = {i}
			else:
				cindex[sv].add(i)
	print('\n\tskipped = %i' % (skipped))
	f.close()
	return kindex

def create_ids(_ki,_sp,_p):
	ids = []
	z = 1
	pt = float(_p['parent mass tolerance'])
	dvals = (-1,0,1)
	count = 0
	for s in _sp:
		ident = []
		pm = int(0.5+s['pm']/pt)
		for d in dvals:
			if pm+d not in _ki:
				continue
			ci = _ki[pm+d]
			for m in s['sms']:
				v = m
				if v in ci:
					ident.append(ci[v])
				elif v-1 in ci:
					ident.append(ci[v-1])
				elif v+1 in ci:
					ident.append(ci[v+1])
		ans = {}
		for i in ident:
			for a in i:
				if a in ans:
					ans[a] += 1
				else:
					ans[a] = 1
		mn = 0
		mv = []
		for j in ans:
			if ans[j] > mn:
				mn = ans[j]
				mv = [j]
			elif ans[j] == mn:
				mv.append(j)
		if mn >= 5:
			ids.append({'sn':z,'peaks':mn,'kernels':mv})
		z += 1
	return ids

def report_ids(_ids):
	for j in _ids:
		print(j)
	print('\tids = %i' % (len(_ids)))

# Coordinate the identification process
def main():
	start = time.time()
	# record relavent parameters
	param = {}
	param['fragment mass tolerance'] = 20
	param['parent mass tolerance'] = 20
	spectra = []
	_in = '/mnt/ssd1/jsms/ai/spectra/GPM06610040149.cmn.mgf'
	print('load spectra')
	# read the spectrum file and perform all necessary spectrum conditioning
	spectra = load_spectra(_in,param)
	delta = time.time()-start
	start = time.time()
	print('\n\tspectra = %i' % (len(spectra)))
	print('\tspectra loading = %.0f s' % (delta))
	fname = '/mnt/ssd1/jsms/ai/kernels/UP000005640_9606.KR.kernel'
	print('index kernel')
	# read the kernel file and create an index of peptide fragmentation patterns
	ki = index_kernel(fname,param,spectra)
	delta = time.time()-start
	start = time.time()
	print('\tkernels = %i' % (len(ki)))
	print('\tkernel indexing = %.0f s' % (delta))
	print('perform ids')
	# generate a list of identifications for the spectra using the kernel index
	ids = create_ids(ki,spectra,param)
	delta = time.time()-start
	start = time.time()
	# simple reporting of the kernels assigned to spectra
	report_ids(ids)
	print('\tid process = %.0f s' % (delta))

if __name__== "__main__":
	main()


