#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
# Updates kernels with residue modification information
#
# idX version 2019.08.10.01
#

import ujson
import time
import gzip
import sys
import copy
import itertools
import hashlib
import re

def sites_process(_os,_of,_mhash,_sdict,_type,_limit):
	mods = {'pho': 79966,'ubi': 114043, 'ace': 42011}
	not_last= {'ace','ubi'}
	mod = mods[_type]
	seq = list(_os['seq'])
	is_ace = False
	if _type == 'ace':
		is_ace = True
	is_ubi = False
	if _type == 'ubi':
		is_ubi = True
	depth = 1
	csites = {}
	if 'mods' in _os:
		csites = {m[1] for m in _os['mods']}
	spos = []
	if _type not in _sdict:
		return 0
	jsites = _sdict[_type]
	beg = _os['beg']
	last = _os['end']+1
	if _type in not_last:
		last -= 1
	for i in range(beg,last):
		if str(i) in jsites and str(i) not in csites:
			if jsites[str(i)] > _limit:
				spos.append(i-beg)
	if len(spos) == 0:
		return 0
	site_list = []
	dlimit = 1
	for d in range(1,dlimit + 1):
		site_list += list(itertools.combinations(spos,d))
	count = 0
	for site in site_list:
		js = copy.deepcopy(_os)
		seq = list(_os['seq'])
		lres = seq.pop()
		bions = js['bs']
		dm = 0
		mods = []
		if 'mods' in js:
			mods = js['mods']
		for j,r in enumerate(seq):
			if j in site:
				if is_ace and j+js['beg'] < 4 and r != 'K':
					mods.append(['[',j+js['beg'],mod])
					dm += mod
				elif is_ace and r != 'K':
					pass
				elif is_ubi and r != 'K':
					pass
				else:
					mods.append([r,j+js['beg'],mod])
					dm += mod
			bions[j] += dm
		j = len(js['seq']) - 1
		if j in site:
			dm += mod
			mods.append([lres,j+js['beg'],mod])
		js['pm'] += dm
		js['bs'] = bions
		yions = js['ys']
		dm = 0
		len_seq = len(seq)
		ysite = [len_seq-s for s in site]
		for j,r in enumerate(seq):
			if j in ysite:
				dm += mod
			yions[j] += dm
		js['ys'] = yions
		js['mods'] = mods
		_mhash.update(ujson.dumps(js).strip().encode())
		_of.write(ujson.dumps(js) + '\n')
		count += 1
	return count


def load_sites(_info):
	f = open(_info['site file'])
	sites = {}
	for i,l in enumerate(f):
		js = ujson.loads(l)
		sites[js['label']] = js['ptms']
	f.close()
	return sites

def update_kernels(_info):
	depth = 0
	mhash = hashlib.sha256()
	in_file = _info['input file']
	out_file = _info['output file']
	if in_file.find('.gz') == len(in_file) - 3:
		f = gzip.open(in_file,'rt', encoding='utf-8')
	else:
		f = open(in_file,'r', encoding='utf-8')
	of = open(out_file,'w',encoding='utf-8')
	site_dict = load_sites(info)
	scount = {'pho': 0, 'ubi': 0 ,'ace': 0}
	types = ['pho','ace','ubi']
	for i,l in enumerate(f):
		# print progress indicator
		if i % 20000 == 0:
			print('.',end='')
			sys.stdout.flush()
		# convert the line into json
		os = ujson.loads(l)
		if 'validation' in os:
			continue
		if 'pm' not in os:
			of.write(l)
			mhash.update(l.strip().encode())
			continue
		of.write(l)
		mhash.update(l.strip().encode())
		lb = re.sub('\w\w\|(\w+)\|',r"\1",os['lb'])
		if lb not in site_dict or site_dict[lb] == None:
			continue
		for t in types:
			scount[t] += sites_process(os,of,mhash,site_dict[lb],t,4)

	js = {"validation" : "sha256", "value" : mhash.hexdigest()}
	string = ujson.dumps(js) + '\n'
	of.write(string)
	f.close()
	of.close()
	print('')
	print(scount)

if len(sys.argv) != 4:
	print('usage:\n\t>python3 special_modify.py INPUT OUTPUT SITES')
	exit()

print('Additional rare modifications')
info = dict()
info['input file'] = sys.argv[1]
info['output file'] = sys.argv[2]
info['site file'] = sys.argv[3]
update_kernels(info)

