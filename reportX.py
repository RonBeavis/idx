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

#
# finds upper and lower limits for identification window
#

def load_mods():
	ls = [x.strip() for x in open('f:/idx/report_mods.txt','r')]
	mods = {}
	for l in ls:
		vs = l.split('\t')
		if len(vs) == 2:
			mods[int(vs[0])] = vs[1]
	return mods
	
def find_window(_ppms):
	vs = {}
	for i in range(-20,21):
		vs[i] = 0
	for s in _ppms:
		i = round(_ppms[s])
		if i >= -20 and i <= 20:
			vs[i] += 1
	center = max(vs, key=vs.get)
	ic = float(vs[center])
	if ic < 100:
		return (-20,20)
	low = center
	for i in range(-20,center):
		if vs[i]/ic >= 0.01:
			low = i
			break
	high = center
	for i in range(20,center,-1):
		if vs[i]/ic >= 0.01:
			high = i
			break
	return (low,high)
#
# applys the current statistical model for evaluating an
# identification
#

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
	mods_translation = load_mods()
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
	inferred = 0
	score_min = 200.0
	#read through the kernel file and extract required information
	ppms = {}
	for c,l in enumerate(f):
		# check to see if line has information in _ids
		if c != 0 and c % 2500 == 0:
			print('.',end='')
			sys.stdout.flush()
		if c != 0 and c % 50000 == 0:
			print(' %i' % (c))
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
			score = apply_model(res,js['seq'],js['pm'],sv[s]['peaks'],sv[s]['ions'])
			if score < score_min or sv[s]['ri'] < 0.20 or sv[s]['peaks'] < 8:
				continue
			oline = '%i\t%.3f\t' % (sv[s]['sc'],js['pm']/1000.0)
			oline += '%.3f\t%.1f\t' % (delta,ppm)
			oline += '%i\t%s\t' % (sv[s]['pz'],js['lb'])
			oline += '%i\t%i\t' % (js['beg'],js['end'])
			oline += '%s\t%s\t' % (js['pre'],js['seq'])
			oline += '%s\t%i\t' % (js['post'],sv[s]['peaks'])
			oline += '%.2f\t%i\t%.0f\t' % (sv[s]['ri'],sum(js['ns']),score)
			if 'mods' in js:
				vmods = []
				for m in sorted(js['mods'],key=itemgetter(1,0)):
					if m[2] in mods_translation:
						vmods.append('%s%i~%s;' % (m[0],m[1],mods_translation[m[2]]))
					else:
						vmods.append('%s%i#%.3f;' % (m[0],m[1],float(m[2])/1000))
				mod_string = ''
				for v in vmods:
					if v.find('[') != -1:
						mod_string += v
				for v in vmods:
					if v.find('[') == -1 and v.find(']') == -1:
						mod_string += v
				for v in vmods:
					if v.find(']') != -1:
						mod_string += v
				oline += re.sub('\;$','',mod_string)
			last_i = s+1
			if s in odict:
				odict[s].append(oline)
			else:
				odict[s] = [oline]
			ppms[s] = ppm
			
	f.close()
	# open the output file specified in _p
	o = open(_p['output file'],'w')
	# create the header line
	oline = 'Id\tScan\tPeptide mass\tDelta\tppm\tz\tProtein acc\t'
	oline += 'Start\tEnd\tPre\tSequence\tPost\tIC\tRI\tFreq\tScore\tModifications\n'
	o.write(oline)
	# output the lines in odict, sorted by id number
	tot = 0
	(low,high) = find_window(ppms)
	err = 0
	for a in sorted(odict):
		sub = 1
		for t in odict[a]:
			ps = t.split('\t')
			if high >= round(float(ps[3])) >= low:
				o.write('%i:%i\t%s\n' % (a,sub,t))
				sub += 1
				tot += 1
			else:
				err += 1
	o.close()
	print('\n\t     lines = %i' % (tot))
	if low != -20 and high != 20:
		ble = 100.0*((high-low)/41.0)*float(err)/float(tot)
		print('\t       ble = %.1f%%' % (ble))
	print('\tppm window = %i,%i' % (low,high))

