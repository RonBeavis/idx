#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
# Identifies kernels corresponding to spectra
#
# idX version 2019.08.10.02
#
from __future__ import print_function
from libcpp cimport bool as bool_t

import ujson
import time
import gzip
import sys
import datetime
import struct
#
# generate identifications of spectra, based on the kernel indexes
#

def create_ids(_ki,_mi,_sp,_p):
	# initialize list of identifications
	ids = []
	z = 1
	pt = 70
	ppm = _p['parent mass tolerance']/1.0e06
	dvals = (-1,0,1)
	count = 0
	c13 = 1003
	# iterate through spectra
	for c,s in enumerate(_sp):
		if c != 0 and c % 50 == 0:
			print('.',end='')
			sys.stdout.flush()
		if c != 0 and c % 1000 == 0:
			print(' %i' % (c))
			sys.stdout.flush()
		ident = []
		ims = []
		rm = float(s['pm'])
		pm = int(0.5+rm/pt)
		pz = int(s['pz'])
		# check _ki for appropriate fragment ion matches
		# and record a list of kernel identifiers for each match
		# and the corresponding fragment ion intensities
		idx = None
		idi = None
		use2 = (pm > 1500000 and pz == 2)
		use3 = (pz > 2)
		ms = []
		for d in dvals:
			if pm+d in _ki:
				ci = _ki[pm+d]
				pass
			else:
				continue
			for b,m in enumerate(s['sms']):
				m2 = m * 2
				if m in ci:
#					idx = list(struct.unpack('!%dI' % (len(ci[m])/4),ci[m]))
					idx = ci[m].split(' ')
					idi = s['ims'][b]
				elif use2 and m2 > 1000000 and m2 in ci:
#					idx = list(struct.unpack('!%dI' % (len(ci[m2])/4),ci[m2]))
					idx = ci[m2].split(' ')
					idi = s['ims'][b]
				elif use3 and m2 in ci:
#					idx = list(struct.unpack('!%dI' % (len(ci[m2])/4),ci[m2]))
					idx = ci[m2].split(' ')
					idi = s['ims'][b]
				else:
					continue
				idx.pop()
				ident.append([int(q) for q in idx])
				ims.append(idi)
				ms.append(rm)
		if rm > 1500000:
			rm -= c13
			cpm = int(0.5+rm/pt)
			for d in dvals:
				if cpm+d in _ki:
					ci = _ki[cpm+d]
					pass
				else:
					continue
				for b,m in enumerate(s['sms']):
					m2 = m * 2
					if m in ci:
#						idx = list(struct.unpack('!%dI' % (len(ci[m])/4),ci[m]))
						idx = ci[m].split(' ')
						idi = s['ims'][b]
					elif use2 and m2 > 1000000 and m2 in ci:
#						idx = list(struct.unpack('!%dI' % (len(ci[m2])/4),ci[m2]))
						idx = ci[m2].split(' ')
						idi = s['ims'][b]
					elif use3 and m2 in ci:
#						idx = list(struct.unpack('!%dI' % (len(ci[m2])/4),ci[m2]))
						idx = ci[m2].split(' ')
						idi = s['ims'][b]
					else:
						continue
					idx.pop()
					ident.append([int(q) for q in idx])
					ims.append(idi)
					ms.append(rm)
		if len(ident) == 0:
			continue
		total = float(sum(s['ims']))/3
		ans = {}
		aint = {}
		aok = set()
		# check kernel identifiers for the exact parent ion mass tolerance
		# and populate aok with identifiers that pass the test
		for b,i in enumerate(ident):
			for a in i:
				if a not in aok and abs(ms[b]-_mi[a]) < ppm*ms[b]:
					aok.add(a)
		# perform an ion count and intensity sum
		# for each kernel identifier
		for b,i in enumerate(ident):
			for a in set(i):
				if a not in aok:
					ans[a] = 0
					aint[a] = 0
				elif a in ans:
					ans[a] += 1
					aint[a] += ims[b]
				else:
					ans[a] = 1
					aint[a] = ims[b]
		mn = 0
		mv = []
		iv = []
		# order the ion and intensity statistics
		for j in ans:
			if ans[j] > mn:
				mn = ans[j]
				mv = [j]
				iv = [aint[j]]
			elif ans[j] == mn:
				mv.append(j)
				iv.append(aint[j])
		# select identifications based on their ion count and summed intensity
		if mn > 4:
			ks = []
			vis = []
			max_i = max(iv)
			for b,m in enumerate(mv):
				if iv[b] < max_i:
					continue
				ks.append(m)
			id_values = {'sn':z,'peaks':mn,'kernels':ks,'ri': max_i/total,'pm': s['pm'],
					'pz': s['pz'],'sc': s['sc'],'ions' : len(s['sms'])/3}
			ids.append(id_values)
		z += 1
	print('')
	return ids

