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

def update_index(_info,_skip):
	mhash = hashlib.sha256()
	in_file = _info['input file']
	out_file = _info['output file']
	if in_file.find('.gz') == len(in_file) - 3:
		f = gzip.open(in_file,'rt', encoding='utf-8')
	else:
		f = open(in_file,'r', encoding='utf-8')
	of = open(out_file,'w',encoding='utf-8')
	res = _info['residue']
	mod_list = _info['modification mass']
	acetyl_skip = _skip
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
		if os['seq'].find(res) == -1:
			of.write(l)
			mhash.update(l.strip().encode())
			continue
		for mod in mod_list:
			js = copy.deepcopy(os)
			seq = list(js['seq'])
			ll = seq.pop()
			bions = js['bs']
			dm = 0
			mods = []
			if 'mods' in js:
				mods = js['mods']
			skip = {m[1]-js['beg'] for m in mods if m[0] == res and m[2] != 114043}
			for j,r in enumerate(seq):
				if r == res:
					if acetyl_skip:
						if j not in skip:
							dm += mod
							mods.append([res,j+js['beg'],mod])

					else:
						dm += mod
						mods.append([res,j+js['beg'],mod])
				bions[j] += dm
			if ll == res:
				if acetyl_skip:
					if len(seq) not in skip:
						dm += mod	
						mods.append([res,js['end'],mod])
				else:
					dm += mod	
					mods.append([res,js['end'],mod])
			js['pm'] += dm
			js['bs'] = bions
			yions = js['ys']
			dm = 0
			seq = list(js['seq'])
			seq.reverse()
			seq.pop()
			lv = len(seq)
			for j,r in enumerate(seq):
				if r == res:
					if acetyl_skip:
						if lv-j not in skip:
							dm += mod
					else:
						dm += mod
				yions[j] += dm
			js['ys'] = yions
			js['mods'] = mods
			string = ujson.dumps(js) + '\n'
			mhash.update(string.strip().encode())
			of.write(string)
	js = {"validation" : "sha256", "value" : mhash.hexdigest()}
	string = ujson.dumps(js) + '\n'
	of.write(string)
	f.close()
	of.close()
	print('')

def update_bnt_index(_info,_v):
	mhash = hashlib.sha256()
	in_file = _info['input file']
	out_file = _info['output file']
	if in_file.find('.gz') == len(in_file) - 3:
		f = gzip.open(in_file,'rt', encoding='utf-8')
	else:
		f = open(in_file,'r', encoding='utf-8')
	of = open(out_file,'w',encoding='utf-8')
	res = _info['residue']
	res = re.sub('\[','',res)
	mod = _info['modification mass'][0]
	for i,l in enumerate(f):
		# print progress indicator
		if i % 20000 == 0:
			print('.',end='')
			sys.stdout.flush()
		# convert the line into json
		js = ujson.loads(l)
		if 'validation' in js:
			continue
		if 'pm' not in js:
			of.write(l)
			mhash.update(l.strip().encode())
			continue
		if _v == 1:
			of.write(l)	
			mhash.update(l.strip().encode())
		seq = list(js['seq'])
		if seq[0] != res:
			continue
		seq.pop()
		bions = js['bs']
		dm = mod
		mods = []
		if 'mods' in js:
			mods = js['mods']
		mods.append(['[',js['beg'],mod])
		for j,r in enumerate(seq):
			bions[j] += dm
		js['bs'] = bions
		js['pm'] += dm
		js['mods'] = mods
		mhash.update(ujson.dumps(js).strip().encode())
		of.write(ujson.dumps(js) + '\n')
	js = {"validation" : "sha256", "value" : mhash.hexdigest()}
	string = ujson.dumps(js) + '\n'
	of.write(string)
	f.close()
	of.close()
	print('')

def update_pnt_index(_info,_v):
	mhash = hashlib.sha256()
	in_file = _info['input file']
	out_file = _info['output file']
	if in_file.find('.gz') == len(in_file) - 3:
		f = gzip.open(in_file,'rt', encoding='utf-8')
	else:
		f = open(in_file,'r', encoding='utf-8')
	of = open(out_file,'w',encoding='utf-8')
	res = _info['residue']
	res = re.sub('\[','',res)
	mod = _info['modification mass'][0]
	for i,l in enumerate(f):
		# print progress indicator
		if i % 20000 == 0:
			print('.',end='')
			sys.stdout.flush()
		# convert the line into json
		js = ujson.loads(l)
		if 'validation' in js:
			continue
		if 'pm' not in js:
			of.write(l)
			mhash.update(l.strip().encode())
			continue
		if _v == 1:
			of.write(l)	
			mhash.update(l.strip().encode())
		seq = list(js['seq'])
		if js['beg'] > 3 or 'FWILHRKEQ'.find(seq[0]) != -1:
			continue
		seq.pop()
		bions = js['bs']
		dm = mod
		mods = []
		if 'mods' in js:
			mods = js['mods']
		mods.append(['[',js['beg'],mod])
		for j,r in enumerate(seq):
			bions[j] += dm
		js['bs'] = bions
		js['pm'] += dm
		js['mods'] = mods
		mhash.update(ujson.dumps(js).strip().encode())
		of.write(ujson.dumps(js) + '\n')
	js = {"validation" : "sha256", "value" : mhash.hexdigest()}
	string = ujson.dumps(js) + '\n'
	of.write(string)
	f.close()
	of.close()
	print('')

def update_nt_index(_info,_v):
	mhash = hashlib.sha256()
	in_file = _info['input file']
	out_file = _info['output file']
	if in_file.find('.gz') == len(in_file) - 3:
		f = gzip.open(in_file,'rt', encoding='utf-8')
	else:
		f = open(in_file,'r', encoding='utf-8')
	of = open(out_file,'w',encoding='utf-8')
	res = _info['residue']
	mod = _info['modification mass'][0]
	for i,l in enumerate(f):
		# print progress indicator
		if i % 20000 == 0:
			print('.',end='')
			sys.stdout.flush()
		# convert the line into json
		js = ujson.loads(l)
		if 'validation' in js:
			continue
		if 'pm' not in js:
			of.write(l)
			mhash.update(l.strip().encode())
			continue
		if _v == 1:
			of.write(l)	
			mhash.update(l.strip().encode())
		seq = list(js['seq'])
		seq.pop()
		bions = js['bs']
		dm = mod
		mods = []
		if 'mods' in js:
			mods = js['mods']
		nt_exists = False
		for m in mods:
			if m[0] == '[':
				nt_exists = True
		if nt_exists:
			if _v != 1:
				of.write(l)
				mhash.update(l.strip().encode())
			continue
		mods.append([res,js['beg'],mod])
		for j,r in enumerate(seq):
			bions[j] += dm
		js['bs'] = bions
		pm = js['pm']
		js['pm'] += dm
		js['mods'] = mods
		mhash.update(ujson.dumps(js).strip().encode())
		of.write(ujson.dumps(js) + '\n')
	js = {"validation" : "sha256", "value" : mhash.hexdigest()}
	string = ujson.dumps(js) + '\n'
	of.write(string)
	f.close()
	of.close()
	print('')

def update_ct_index(_info,_v):
	mhash = hashlib.sha256()
	in_file = _info['input file']
	out_file = _info['output file']
	if in_file.find('.gz') == len(in_file) - 3:
		f = gzip.open(in_file,'rt', encoding='utf-8')
	else:
		f = open(in_file,'r', encoding='utf-8')
	of = open(out_file,'w',encoding='utf-8')
	res = _info['residue']
	mod = _info['modification mass'][0]
	for i,l in enumerate(f):
		# print progress indicator
		if i % 20000 == 0:
			print('.',end='')
			sys.stdout.flush()
		# convert the line into json
		js = ujson.loads(l)
		if 'validation' in js:
			continue
		if 'pm' not in js:
			of.write(l)
			mhash.update(l.strip().encode())
			continue
		if _v == 1:
			of.write(l)	
			mhash.update(l.strip().encode())
		seq = list(js['seq'])
		seq.reverse()
		seq.pop()
		yions = js['ys']
		dm = mod
		mods = []
		ct_exists = False
		if 'mods' in js:
			mods = js['mods']
		for m in mods:
			if m[0] == ']':
				ct_exists = True
		if ct_exists:
			if _v != 1:
				of.write(l)
				mhash.update(l.strip().encode())
			continue
		mods.append([res,js['end'],mod])
		for j,r in enumerate(seq):
			yions[j] += dm
		js['ys'] = yions
		js['pm'] += dm
		js['mods'] = mods
		mhash.update(ujson.dumps(js).strip().encode())
		of.write(ujson.dumps(js) + '\n')
	js = {"validation" : "sha256", "value" : mhash.hexdigest()}
	string = ujson.dumps(js) + '\n'
	of.write(string)
	f.close()
	of.close()
	print('')

def update_vindex(_info,_d):
	depth = _d
	print('depth = %i' % (depth))
	mhash = hashlib.sha256()
	in_file = _info['input file']
	out_file = _info['output file']
	if in_file.find('.gz') == len(in_file) - 3:
		f = gzip.open(in_file,'rt', encoding='utf-8')
	else:
		f = open(in_file,'r', encoding='utf-8')
	of = open(out_file,'w',encoding='utf-8')
	res = _info['residue']
	mod = _info['modification mass'][0]
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
		if os['seq'].find(res) == -1:
			of.write(l)
			mhash.update(l.strip().encode())
			continue
		of.write(l)
		mhash.update(l.strip().encode())
		seq = list(os['seq'])
		spos = []
		for j,r in enumerate(seq):
			if r == res:
				spos.append(j)
		site_list = []
		dlimit = depth
		if len(spos) < depth:
			dlimit = len(spos)
		for d in range(1,dlimit + 1):
			site_list += list(itertools.combinations(spos,d))
		for site in site_list:
			js = copy.deepcopy(os)
			seq = list(os['seq'])
			seq.pop()
			bions = js['bs']
			dm = 0
			mods = []
			if 'mods' in js:
				mods = js['mods']
			for j,r in enumerate(seq):
				if j in site:
					dm += mod
					mods.append([res,j+js['beg'],mod])
				bions[j] += dm
			j = len(js['seq']) - 1
			if j in site:
				dm += mod
				mods.append([res,j+js['beg'],mod])
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
	js = {"validation" : "sha256", "value" : mhash.hexdigest()}
	string = ujson.dumps(js) + '\n'
	of.write(string)
	f.close()
	of.close()
	print('')

if len(sys.argv) != 5:
	print('''usage:\n\t>python3 kernel_modifier.py INPUT OUTPUT mass@X -(vfs)
	-v variable mod
	-f fixed mod
	-s fixed mod, skip if already modified, unless it is K-ubiquitin''')
	exit()

info = dict()
info['input file'] = sys.argv[1]
info['output file'] = sys.argv[2]
v = sys.argv[3].split('@')

ms = v[0].split(',')
info['modification mass'] = []
for m in ms:
	info['modification mass'].append(int(0.5 + float(m)*1000.0))
info['residue'] = v[1]
for i in info:
	print('%s:\t%s' % (i,str(info[i])))
if sys.argv[4] == '-f':
	if info['residue'] == '[':
		update_nt_index(info,0)
	elif len(info['residue']) == 2 and info['residue'].find('[') == 0:
		update_bnt_index(info,0)
	elif info['residue'] == '^':
		update_pnt_index(info,0)
	elif info['residue'] == ']':
		update_ct_index(info,0)
	else:
		update_index(info,False)
if sys.argv[4] == '-s':
	if info['residue'] == '[':
		update_nt_index(info,0)
	elif len(info['residue']) == 2 and info['residue'].find('[') == 0:
		update_bnt_index(info,0)
	elif info['residue'] == '^':
		update_pnt_index(info,0)
	elif info['residue'] == ']':
		update_ct_index(info,0)
	else:
		update_index(info,True)
elif sys.argv[4] == '-v':
	if info['residue'] == '[':
		update_nt_index(info,1)
	elif len(info['residue']) == 2 and info['residue'].find('[') == 0:
		update_bnt_index(info,1)
	elif info['residue'] == '^':
		update_pnt_index(info,1)
	elif info['residue'] == ']':
		update_ct_index(info,1)
	else:
		vs = list(info['residue'])
		if len(vs) == 1:
			update_vindex(info,1)
		else:
			info['residue'] = vs[1]
			update_vindex(info,int(vs[0]))


