#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#

#
# Identifies kernels corresponding to spectra
#
# Reads kernels in either plain text or gzip'd text
# idX version 2019.08.10.01
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

def index_kernel(_p,_s):
	ft = float(_p['fragment mass tolerance'])
	pt = float(_p['parent mass tolerance'])
	sp_set = set()
	fname = _p['kernel file']
	for s in _s:
		sp_set.add(int(0.5 + s['pm']/pt))
	kindex = {}
	if fname.find('.gz') == len(fname) - 3:
		f = gzip.open(fname,'rt', encoding='utf-8')
	else:
		f = open(fname,'r', encoding='utf-8')
	skipped = 0
	for i,l in enumerate(f):
		if i % 10000 == 0:
			print('.',end='')
			sys.stdout.flush()
		js = json.loads(l)
		if 'pm' not in js:
			continue
		pm = js['pm']
		mv = int(0.5+pm/pt)
		if not (mv in sp_set or mv-1 in sp_set or mv+1 in sp_set):
			skipped += 1
			continue
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
		ys = js['ys']
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
		ims = []
		pm = int(0.5+s['pm']/pt)
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
		for b,i in enumerate(ident):
			for a in i:
				if a in ans:
					ans[a] += 1
					aint[a] += ims[b]
				else:
					ans[a] = 1
					aint[a] = ims[b]
		mn = 0
		mv = []
		iv = []
		for j in ans:
			if ans[j] > mn:
				mn = ans[j]
				mv = [j]
				iv = [aint[j]]
			elif ans[j] == mn:
				mv.append(j)
				iv.append(aint[j])
		if mn > 7:
			ks = []
			vis = []
			max_i = max(iv)
			for b,m in enumerate(mv):
				if iv[b] < max_i:
					continue
				ks.append(m)
			if max_i/total >= 0.20:
				ids.append({'sn':z,'peaks':mn,'kernels':ks,'ri': max_i/total,'pm': s['pm'],'pz': s['pz']})
		z += 1
	return ids

def report_ids(_ids,_p):
	sdict = {}
	scores = {}
	ivalues = {}
	zvalues = {}
	pmvalues = {}
	for j in _ids:
		scores[j['sn']] = j['peaks']
		ivalues[j['sn']] = j['ri']
		zvalues[j['sn']] = j['pz']
		pmvalues[j['sn']] = j['pm']
		for k in j['kernels']:
			if k in sdict:
				sdict[k].add(j['sn'])
			else:
				sdict[k] = {j['sn']}
	fname = _p['kernel file']
	odict = {}
	if fname.find('.gz') == len(fname) - 3:
		f = gzip.open(fname,'rt', encoding='utf-8')
	else:
		f = open(fname,'r', encoding='utf-8')
	for i,l in enumerate(f):
		if i not in sdict:
			continue
		js = json.loads(l)
		if 'pm' not in js:
			continue
		spectra = sdict[i]
		for s in spectra:
			oline = '%i\t%.3f\t%.3f\t%i\t%s\t%i\t%i\t%s\t%i\t%.2f\t%i' % (s,js['pm']/1000.0,(pmvalues[s]-js['pm'])/1000.0,zvalues[s],js['lb'],js['beg'],js['end'],js['seq'],scores[s],ivalues[s],sum(js['ns']))
			if s in odict:
				odict[s].append(oline)
			else:
				odict[s] = [oline]
			
	f.close()
	o = open(_p['output file'],'w')
	oline = 'Spectrum\tPeptide mass\tDelta\tz\tProtein acc\tStart\tEnd\tSequence\tScore\tRI\tFreq'
	o.write(oline + '\n')
	for a in sorted(odict):
		for t in odict[a]:
			o.write(t + '\n')
	o.close()
	print('\tids = %i' % (len(_ids)))

# Coordinate the identification process
def main():
	start = time.time()
	# record relavent parameters
	param = {}
	param['fragment mass tolerance'] = 20
	param['parent mass tolerance'] = 20
	spectra = []
	param['spectrum file'] = sys.argv[1]
	print('load spectra')
	# read the spectrum file and perform all necessary spectrum conditioning
	spectra = load_spectra(param['spectrum file'],param)
	delta = time.time()-start
	start = time.time()
	print('\n\tspectra = %i' % (len(spectra)))
	print('\tspectra loading = %.1f s' % (delta))
	param['kernel file'] = sys.argv[2]
	print('index kernel')
	# read the kernel file and create an index of peptide fragmentation patterns
	ki = index_kernel(param,spectra)
	delta = time.time()-start
	start = time.time()
	print('\tkernels = %i' % (len(ki)))
	print('\tkernel indexing = %.1f s' % (delta))
	print('perform ids')
	# generate a list of identifications for the spectra using the kernel index
	ids = create_ids(ki,spectra,param)
	# free memory associated with indexes and spectra
	delta = time.time()-start
	start = time.time()
	print('\tid process = %.3f s' % (delta))
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


