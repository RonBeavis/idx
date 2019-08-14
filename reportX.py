#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
# Creates an output file from an idX job
#
# idX version 2019.08.10.02
#

import ujson
import gzip
import re
import sys
import math
from scipy.stats import hypergeom,tmean,tstd
from operator import itemgetter

def apply_model(_res,_seq,_pm,_ions,_lspectrum):
	sfactor = 20
	sadjust = 1
	if _res > 100:
		sfactor = 40
	lseq = list(_seq)
	pmass = int(_pm/1000)
	cells = int(pmass-200)
	if cells > 1500:
		cells = 1500
	total_ions = 2*(len(lseq) - 1)
	if total_ions > sfactor:
		total_ions = sfactor
	if total_ions < _ions:
		total_ions = _ions + 1
	sc = _lspectrum
	if _ions >= sc:
		sc = _ions + 2
	rv = hypergeom(cells,total_ions,sc)
	p = rv.pmf(_ions)
	pscore = -100.0*math.log10(p)*sadjust
	return pscore

#
# formats the output of an idX job into a TSV file
#

def report_ids(_ids,_p):
	sdict = {}
	sv = {}
	# create an index of the results in _ids
	for j in _ids:
		sv[j['sn']] = j
		for k in j['kernels']:
			if k in sdict:
				sdict[k].add(j['sn'])
			else:
				sdict[k] = {j['sn']}
	# load the kernel file corresponding to the results
	fname = _p['kernel file']
	res = _p['fragment mass tolerance']
	odict = {}
	if fname.find('.gz') == len(fname) - 3:
		f = gzip.open(fname,'rt', encoding='utf-8')
	else:
		f = open(fname,'r', encoding='utf-8')
	sub = dict()
	inferred = 0
	score_min = 200.0
	#read through the kernel file and extract required information
	for c,l in enumerate(f):
		# check to see if line has information in _ids
		if c % 20000 == 0:
			print('.',end='')
			sys.stdout.flush()
		js = ujson.loads(l)
		if 'pm' not in js:
			continue
		h = js['h']
		if h not in sdict:
			continue
		if h != js['u']:
			inferred += 1
		spectra = sdict[h]
		last_i = 0
		# create lines of output information
		# and store in odict
		for s in spectra:
			delta = (sv[s]['pm']-js['pm'])/1000.0
			ppm = 1.0e6*(sv[s]['pm']-js['pm'])/js['pm']
			# update the "sub" identifier for an identification
			# used to distiguish ids when there is more than
			# oene kernel identified for a particular spectrum
			if s+1 in sub:
				sub[s+1] += 1
			else:
				sub[s+1] = 1
			score = apply_model(res,js['seq'],js['pm'],sv[s]['peaks'],sv[s]['ions'])
			if score < score_min or sv[s]['ri'] < 0.20 or sv[s]['peaks'] < 8:
				continue
			oline = '%i.%i\t' % (s,sub[s+1])
			oline += '%i\t%.3f\t' % (sv[s]['sc'],js['pm']/1000.0)
			oline += '%.3f\t%.1f\t' % (delta,ppm)
			oline += '%i\t%s\t' % (sv[s]['pz'],js['lb'])
			oline += '%i\t%i\t' % (js['beg'],js['end'])
			oline += '%s\t%s\t' % (js['pre'],js['seq'])
			oline += '%s\t%i\t' % (js['post'],sv[s]['peaks'])
			oline += '%.2f\t%i\t%.0f\t' % (sv[s]['ri'],sum(js['ns']),score)
			if 'mods' in js:
				mod_string = ''
				for m in sorted(js['mods'],key=itemgetter(1)):
					mod_string += '%s%i#%.3f;' % (m[0],m[1],m[2]/1000.0)
				oline += re.sub('\;$','',mod_string)
			last_i = s+1
			if s in odict:
				odict[s].append(oline)
			else:
				odict[s] = [oline]
			
	f.close()
	# open the output file specified in _p
	o = open(_p['output file'],'w')
	# create the header line
	oline = 'Id\tScan\tPeptide mass\tDelta\tppm\tz\tProtein acc\tStart\tEnd\tPre\tSequence\tPost\tIC\tRI\tFreq\tScore\tModifications'
	o.write(oline + '\n')
	# output the lines in odict, sorted by id number
	total = 0
	for a in sorted(odict):
		for t in odict[a]:
			o.write(t + '\n')
			total += 1
	o.close()
	print('\n\tlines = %i' % (total))
	print('\tidenticals = %i' % (inferred))

