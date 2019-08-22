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

def collagen_process(os,of,mhash):
	seq = list(os['seq'])
	seq_plus = os['pre']+os['seq']+os['post']
	depth = len(re.findall('(?=(G.PG))',seq_plus))
	if depth == 0:
		return 0
	depth = len(re.findall('(?=(G.[PK]G))',seq_plus))
	if depth > 10:
		depth = 10
	iter = re.finditer(r"G.PG",seq_plus)
	spos = [m.start(0)+1 for m in iter]
	iter = re.finditer(r"G.KG",seq_plus)
	spos.extend([m.start(0)+1 for m in iter])
	site_list = []
	dlimit = depth
	if len(spos) < depth:
		dlimit = len(spos)
	for d in range(1,dlimit + 1):
		site_list += list(itertools.combinations(spos,d))
	mod = 15995
	count = 0
	for site in site_list:
		js = copy.deepcopy(os)
		seq = list(os['seq'])
		lres = seq.pop()
		bions = js['bs']
		dm = 0
		mods = []
		if 'mods' in js:
			mods = js['mods']
		for j,r in enumerate(seq):
			if j in site:
				dm += mod
				mods.append([r,j+js['beg'],mod])
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
		mhash.update(ujson.dumps(js).strip().encode())
		of.write(ujson.dumps(js) + '\n')
		count += 1
	return count

def ubi_process(os,of,mhash):
	ubi_proteins = {'sp|P0CG47|','sp|P0CG49|','YLL039C'}
	if os['lb'] not in ubi_proteins:
		return 0
	seq = list(os['seq'])
	seq_plus = os['seq']
	depth = len(re.findall('(?=(K.))',seq_plus))
	if depth == 0:
		return 0
	iter = re.finditer(r"K.",seq_plus)
	spos = [m.start(0) for m in iter]
	site_list = []
	dlimit = 1
	for d in range(1,dlimit + 1):
		site_list += list(itertools.combinations(spos,d))
	mod = 114043
	count = 0
	for site in site_list:
		js = copy.deepcopy(os)
		seq = list(os['seq'])
		lres = seq.pop()
		bions = js['bs']
		dm = 0
		mods = []
		if 'mods' in js:
			mods = js['mods']
		for j,r in enumerate(seq):
			if j in site:
				dm += mod
				mods.append([r,j+js['beg'],mod])
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
		mhash.update(ujson.dumps(js).strip().encode())
		of.write(ujson.dumps(js) + '\n')
		count += 1
	return count

def ng_process(os,of,mhash):
	seq = list(os['seq'])
	seq_plus = os['seq']+os['post']
	depth = len(re.findall('(?=(NG))',seq_plus))
	if depth == 0:
		return 0
	depth = 1
	iter = re.finditer(r"NG",seq_plus)
	spos = [m.start(0) for m in iter]
	site_list = []
	dlimit = 1
	for d in range(1,dlimit + 1):
		site_list += list(itertools.combinations(spos,d))
	mod = 984
	count = 0
	for site in site_list:
		js = copy.deepcopy(os)
		seq = list(os['seq'])
		lres = seq.pop()
		bions = js['bs']
		dm = 0
		mods = []
		if 'mods' in js:
			mods = js['mods']
		for j,r in enumerate(seq):
			if j in site:
				dm += mod
				mods.append([r,j+js['beg'],mod])
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
		mhash.update(ujson.dumps(js).strip().encode())
		of.write(ujson.dumps(js) + '\n')
		count += 1
	return count

def rg_process(os,of,mhash):
	seq = list(os['seq'])
	seq_plus = os['seq']
	depth = len(re.findall('(?=(RG))',seq_plus))
	if depth == 0:
		return 0
	depth = 1
	iter = re.finditer(r"RG",seq_plus)
	spos = [m.start(0) for m in iter]
	site_list = []
	dlimit = 1
	for d in range(1,dlimit + 1):
		site_list += list(itertools.combinations(spos,d))
	mod = 28031
	count = 0
	for site in site_list:
		js = copy.deepcopy(os)
		seq = list(os['seq'])
		lres = seq.pop()
		bions = js['bs']
		dm = 0
		mods = []
		if 'mods' in js:
			mods = js['mods']
		for j,r in enumerate(seq):
			if j in site:
				dm += mod
				mods.append([r,j+js['beg'],mod])
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
		mhash.update(ujson.dumps(js).strip().encode())
		of.write(ujson.dumps(js) + '\n')
		count += 1
	return count

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
	hpro = 0
	ng = 0
	ug = 0
	rg = 0
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
		hpro += collagen_process(os,of,mhash)
		ng += ng_process(os,of,mhash)
#		ug += ubi_process(os,of,mhash)
		rg += rg_process(os,of,mhash)

	js = {"validation" : "sha256", "value" : mhash.hexdigest()}
	string = ujson.dumps(js) + '\n'
	of.write(string)
	f.close()
	of.close()
	print('')
	print('hpro = %i, ng = %i, ug = %i, rg = %i' % (hpro,ng,ug,rg))

if len(sys.argv) != 3:
	print('usage:\n\t>python3 special_modify.py INPUT OUTPUT')
	exit()

print('Additional rare modifications')
info = dict()
info['input file'] = sys.argv[1]
info['output file'] = sys.argv[2]

update_kernels(info)

