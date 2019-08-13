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
	odict = {}
	if fname.find('.gz') == len(fname) - 3:
		f = gzip.open(fname,'rt', encoding='utf-8')
	else:
		f = open(fname,'r', encoding='utf-8')
	sub = dict()
	inferred = 0
	#read through the kernel file and extract required information
	for c,l in enumerate(f):
		# check to see if line has information in _ids
		if c % 10000 == 0:
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
			# one kernel identified for a particular spectrum
			if s+1 in sub:
				sub[s+1] += 1
			else:
				sub[s+1] = 1
			oline = '%i.%i\t%i\t%.3f\t%.3f\t%.1f\t%i\t%s\t%i\t%i\t%s\t%i\t%.2f\t%i\t' % (
				s+1,sub[s+1],sv[s]['sc'],js['pm']/1000.0,delta,ppm,sv[s]['pz'],
				js['lb'],js['beg'],js['end'],js['seq'],sv[s]['peaks'],sv[s]['ri'],sum(js['ns']))
			if 'mods' in js:
				mod_string = ''
				for m in js['mods']:
					mod_string += '%s%i#%.3f,' % (m[0],m[1],m[2]/1000.0)
				oline += re.sub('\,$','',mod_string)
			last_i = s+1
			if s in odict:
				odict[s].append(oline)
			else:
				odict[s] = [oline]
			
	f.close()
	# open the output file specified in _p
	o = open(_p['output file'],'w')
	# create the header line
	oline = 'Id\tScan\tPeptide mass\tDelta\tppm\tz\tProtein acc\tStart\tEnd\tSequence\tScore\tRI\tFreq\tModifications'
	o.write(oline + '\n')
	# output the lines in odict, sorted by id number
	for a in sorted(odict):
		for t in odict[a]:
			o.write(t + '\n')
	o.close()
	print('\n\tids = %i' % (len(_ids)))
	print('\tidenticals = %i' % (inferred))

