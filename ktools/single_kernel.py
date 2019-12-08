import json
import re
import datetime
import hashlib
import mysql.connector
import sys
import random

mhash = hashlib.sha256()
dhash = hashlib.sha256()

isotopes = {'p' : 1.007276,
'H' : 1.007825,
'C' : 12.0,
'N' : 14.003074,
'O': 15.994915 }

a_to_m = {"A":71.037114,
"R":156.101111,
"B":114.042927,
"N":114.042927,
"D":115.026943,
"C":103.009185,
"E":129.042593,
"Q":128.058578,
"Z":128.058578,
"G":57.021464,
"O":237.147727,
"H":137.058912,
"I":113.084064,
"J":113.084064,
"L":113.084064,
"K":128.094963,
"M":131.040485,
"F":147.068414,
"P":97.052764,
"S":87.032028,
"T":101.047679,
"U":150.95363,
"W":186.079313,
"Y":163.06332,
"V":99.068414 }

js= {}
js['lv'] = 0
if len(sys.argv) == 1 or sys.argv[1] == "?" or sys.argv[1] == "-h":
	print('usage: >single_kernel.py "lb=sp|ALBU_BOVIN|" seq=RDTHKSEIAHR pre=R post=F beg=24 end=34 mod=28,42.011 mod=24,6.000')
	exit()
mods = []
for l in sys.argv[1:]:
	eq = l.find('=')
	if eq == -1:
		print('item %s improperly formatted' % (l))
		exit()
	eq += 1
	if l.find('lb=') == 0:
		js['lb'] = l[eq:]
	elif l.find('seq=') == 0:
		js['seq'] = l[eq:]
	elif l.find('pre=') == 0:
		js["pre"] = l[eq:]
	elif l.find('post=') == 0:
		js["post"] = l[eq:]
	elif l.find('beg=') == 0:
		js["beg"] = int(l[eq:])
	elif l.find('end=') == 0:
		js["end"] = int(l[eq:])
	elif l.find('mod=') == 0:
		ms = l[eq:].split(',')
		mods.append([int(ms[0]),float(ms[1])])
js['mods'] = []
msum = {}
for m in mods:
	t = []
	t.append(js['seq'][m[0]-js['beg']])
	t.append(m[0])
	t.append(int(round(1000.0*m[1],0)))
	js['mods'].append(t)
	if m[0] in msum:
		msum[m[0]] += m[1]
	else:
		msum[m[0]] = m[1]

dB = 0.0
dY = 2*isotopes['H'] + isotopes['O']
water = 2*isotopes['H'] + isotopes['O']
proton = isotopes['p']

frags = []
ion = dB
bions = []
parent = water
mr = js['beg']
pep = js['seq']
for p in pep:
	ptm = 0
	if mr in msum:
		ptm = msum[mr]
	ion += a_to_m[p] + ptm
	parent += a_to_m[p] + ptm
	bions.append(int(round(1000.0*ion,0)))
	mr += 1
bions = bions[:-1]
yions = []
ion = dY
mr = js['end']
for p in pep[::-1]:
	ptm = 0
	if mr in msum:
		ptm = msum[mr]
	ion += a_to_m[p] + ptm
	yions.append(int(round(1000.0*ion,0)))
	mr -= 1
yions = yions[:-1]
ns = [1,1,1,1]
js['ns'] = ns
js['ys'] = yions
js['bs'] = bions
js['pm'] = int(round(1000.0*parent,0))
print(json.dumps(js))

