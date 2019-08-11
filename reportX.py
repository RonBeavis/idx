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

import ujson
#import json
import time
import re
import gzip
import sys
import itertools
import copy
def report_ids(_ids,_p):
	sdict = {}
	sv = {}
	for j in _ids:
		sv[j['sn']] = j
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
	sub = dict()
	for i,l in enumerate(f):
		if i not in sdict:
			continue
		js = ujson.loads(l)
		if 'pm' not in js:
			continue
		spectra = sdict[i]
		last_i = 0
		for s in spectra:
			delta = (sv[s]['pm']-js['pm'])/1000.0
			ppm = 1.0e6*(sv[s]['pm']-js['pm'])/js['pm']
			if s+1 in sub:
				sub[s+1] += 1
			else:
				sub[s+1] = 1
			oline = '%i.%i\t%i\t%.3f\t%.3f\t%.1f\t%i\t%s\t%i\t%i\t%s\t%i\t%.2f\t%i' % (
				s+1,sub[s+1],sv[s]['sc'],js['pm']/1000.0,delta,ppm,sv[s]['pz'],
				js['lb'],js['beg'],js['end'],js['seq'],sv[s]['peaks'],sv[s]['ri'],sum(js['ns']))
			last_i = s+1
			if s in odict:
				odict[s].append(oline)
			else:
				odict[s] = [oline]
			
	f.close()
	o = open(_p['output file'],'w')
	oline = 'Id\tScan\tPeptide mass\tDelta\tppm\tz\tProtein acc\tStart\tEnd\tSequence\tScore\tRI\tFreq'
	o.write(oline + '\n')
	for a in sorted(odict):
		for t in odict[a]:
			o.write(t + '\n')
	o.close()
	print('\tids = %i' % (len(_ids)))

