#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
# Creates an output file from an idX job
#
# idX version 2019.08.10.02
#
#from __future__ import print_function
#from libcpp cimport bool as bool_t

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

distribution = {800:[1.3,1.0,0.6],
		900:[3.2,2.0,0.8],
		1000:[4.9,2.7,1.0],
		1100:[5.6,3.1,1.0],
		1200:[6.2,3.3,1.1],
		1300:[6.6,3.5,1.1],
		1400:[7.1,3.8,1.2],
		1500:[7.5,4.0,1.2],
		1600:[7.8,4.1,1.3],
		1700:[8.1,4.3,1.3],
		1800:[8.2,4.4,1.3],
		1900:[8.4,4.5,1.3],
		2000:[8.5,4.6,1.3],
		2100:[8.4,4.7,1.4],
		2200:[8.6,4.7,1.4],
		2300:[8.3,4.7,1.4],
		2400:[8.2,4.6,1.4],
		2500:[8.2,4.7,1.4],
		2600:[8.1,4.8,1.4],
		2700:[7.8,4.7,1.4],
		2800:[7.5,4.7,1.4],
		2900:[7.3,4.6,1.4],
		3000:[7.3,4.7,1.5],
		3100:[6.9,4.5,1.4],
		3200:[6.5,4.4,1.4],
		3300:[6.4,4.4,1.4],
		3400:[5.8,4.1,1.4],
		3500:[5.0,3.8,1.4],
		3600:[5.5,4.1,1.4],
		3700:[4.7,3.6,1.4],
		3800:[4.2,3.3,1.4],
		3900:[3.8,3.1,1.3],
		4000:[3.7,3.0,1.3]	}

def load_mods():
	ls = [x.strip() for x in open('report_mods.txt','r')]
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

def get_cells(_pm,_res):
	global distribution
	pm = 100*int(_pm/100000)
	r = 2
	if _res == 50:
		r = 1
	elif _res == 20:
		r = 0
	if pm < 800:
		return int(pm*distribution[800][r])
	if pm > 4000:
		return int(pm*distribution[4000][r])
	return int(pm*distribution[pm][r])

#
# applys the current statistical model for evaluating an
# identification
#

def apply_model(_res,_seq,_pm,_ions,_lspectrum):
	sfactor = 20
	sadjust = 3
	if _res > 100:
		sfactor = 40
	lseq = list(_seq)

	cells = get_cells(_pm,_res)
	total_ions = 2*(len(lseq) - 1)
	if total_ions > sfactor:
		total_ions = sfactor
	if total_ions < _ions:
		total_ions = _ions + 1
	sc = _lspectrum * sadjust
	if _ions >= sc:
		sc = _ions + 2
	hyper = hypergeom(cells,total_ions,sc)
	p = hyper.pmf(_ions)
	pscore = -100.0*math.log(p)/2.3
#	print('%s: p=%.1e, cells=%i, total=%i, sc=%i, ions=%i, score=%.0f' % (_seq,p,cells,total_ions,sc,_ions,pscore))
	return (pscore,p)

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
	specs = _p['spectra']
	score_min = 200.0
	if specs > 0:
		score_min += 100*math.log(specs)/2.3
	#read through the kernel file and extract required information
	total_prob = 0.0
	ppms = {}
	min_c = 8
	c13 = 1003
	if res == 50:
		min_c = 7
	elif res == 20:
		min_c = 6
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
		max_prob = 0
		prob = 0
		for s in spectra:
			delta = (sv[s]['pm']-js['pm'])/1000.0
			ppm = 1.0e6*(sv[s]['pm']-js['pm'])/js['pm']
			if delta > 0.5:
				ppm = 1.0e6*(sv[s]['pm']-c13-js['pm'])/js['pm']
			# update the "sub" identifier for an identification
			# used to distiguish ids when there is more than
			# oene kernel identified for a particular spectrum
			(score,prob) = apply_model(res,js['seq'],js['pm'],sv[s]['peaks'],sv[s]['ions'])
			if score < score_min or sv[s]['ri'] < 0.20 or sv[s]['peaks'] < min_c:
				continue
			if prob > max_prob:
				max_prob = prob
			oline = '%i\t%.3f\t' % (sv[s]['sc'],js['pm']/1000.0)
			oline += '%.3f\t%.1f\t' % (delta,ppm)
			oline += '%i\t%s\t' % (sv[s]['pz'],js['lb'])
			oline += '%i\t%i\t' % (js['beg'],js['end'])
			oline += '%s\t%s\t' % (js['pre'],js['seq'])
			oline += '%s\t%i\t' % (js['post'],sv[s]['peaks'])
			if sum(js['ns']) > 0:
				oline += '%.2f\t%.1f\t%.1f\t' % (sv[s]['ri'],math.log(sum(js['ns']))/2.3,-0.01*score)
			else:
				oline += '%.2f\t-\t%.1f\t' % (sv[s]['ri'],-0.01*score)
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
			oline = re.sub(r'\tinf\t','\t5000\t',oline)
			if s in odict:
				odict[s].append(oline)
			else:
				odict[s] = [oline]
			ppms[s] = ppm
		total_prob += max_prob
			
	f.close()
	# open the output file specified in _p
	o = open(_p['output file'],'w')
	# create the header line
	oline = 'Id\tScan\tPeptide mass\tDelta\tppm\tz\tProtein acc\t'
	oline += 'Start\tEnd\tPre\tSequence\tPost\tIC\tRI\tlog(f)\tlog(p)\tModifications\n'
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
	print('\n     lines = %i' % (tot))
	if total_prob > 0:
		print('     fpr = %.1e' % (total_prob*len(spectra)))
	if low != -20 and high != 20:
		ble = 100.0*((high-low)/41.0)*float(err)/float(tot)
		print('     baseline error = %.1f%%' % (ble))
	else:
		print('     baseline error = n/a')
	print('     parent ion tolerance = %i,%i' % (low,high))

