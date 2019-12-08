#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
# Updates kernels with residue modification information
#
# idX version 2019.08.10.01
#

import ujson
import hashlib
import re
import sys
import datetime

uids = {}

def get_hvalue(_js):
	global uids
	u = _js['u']
	ms = []
	if 'mods' in _js:
		ms = _js['mods']
	mstring = _js['seq'] + ' '
	beg = _js['beg']
	for m in ms:
		mstring += '%s %i %i;' % (m[0],m[1]-beg,m[2])
	if mstring in uids:
		return uids[mstring]
	uids[mstring] = u
	return u

if len(sys.argv) < 3:
	print('usage:\n\t>python3 OUTPUT IN1 IN2 IN3 ... (-f)')
	exit()

ofile = sys.argv[1]
i = 2
files = []
finalize = False
while i < len(sys.argv):
	if sys.argv[i] == '-f':
		finalize = True
		print('finalize kernel set')
		i += 1
		continue
	files.append(sys.argv[i])
	i += 1

mhash = hashlib.sha256()
of = open(ofile,'w')
js = {'format': 'jsms 1.0', 'source': ','.join(files), 'created': str(datetime.datetime.now())}
l = ujson.dumps(js)
of.write(l + '\n')
mhash.update(l.strip().encode())

u = 1
hcount = 0
for f in files:
	print('Adding "%s" to "%s"' % (f,ofile))
	jf = open(f,'r')
	for i,l in enumerate(jf):
		if i % 20000 == 0:
			print('.',end='')
			sys.stdout.flush()
		js = ujson.loads(l)
		if 'pm' not in js:
			continue
		cnts = {}
		for s in js['seq']:
			if s in cnts:
				cnts[s] += 1
			else:
				cnts[s] = 1
		ls = len(js['seq'])
		bail = False
		for s in cnts:
			if cnts[s] == ls:
				bail = True
				break
		if bail:
			print(js['seq'])
			continue
		js['u'] = u
		h = get_hvalue(js)
		js['h'] = h
		if u != h:
			if finalize:
				del js['bs']
				del js['ys']
			hcount += 1
		li = ujson.dumps(js) + '\n'
		of.write(li)
		mhash.update(li.strip().encode())
		u += 1
	print('')
	jf.close()

js = {"validation" : "sha256", "value" : mhash.hexdigest()}
l = ujson.dumps(js)
of.write(l + '\n')
of.close()
print('total = %i, duplicates = %i' % (u-1,hcount))
