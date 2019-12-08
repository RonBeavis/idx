import os
import sys

output = 'k/hs.KR.kernel'
base = 'hs.KR'
use_special = True
cys_mod = [57.021464]
modify_u = True

ok = False
while not ok:
	ipath = input('\ninput kernel: ')
	output = input('output file (e.g. k/hs.mhc1.kernel): ')
	base = input('base name (e.g. hs.mhc1): ')
	value = input('use special mods (yes/no): ')
	if value == 'no':
		use_special = False
	value = input('cysteine modification (e.g. 57.021464,119.004099): ')
	cys = value.split(',')
	cys_mod = []
	for c in cys:
		cys_mod.append(c)
	value = input('modify U same as C (yes/no): ')
	if value == 'no':
		modify_u = False

	print('\ninput path: %s' % (str(ipath)))
	print('output file: %s' % (str(output)))
	print('base name: %s' % (str(base)))
	print('use special mods: %s' % (str(use_special)))
	print('cysteine modification: %s' % (str(cys_mod)))
	print('modify U same as C: %s' % (str(modify_u)))
	value = input('\n\tok (yes/no/exit): ')
	if value == 'exit':
		exit()
	elif value == 'no':
		ok = False
	else:
		ok = True

alt = '.iac'
out_file = 'k/%s%s.kernel' % (base,alt)
tfiles = set({out_file})
os.system('python3 ktools/modify.py %s %s %s@C -f' % (ipath,out_file,','.join(cys_mod)))

if modify_u:
	base += alt
	alt = '.iau'
	out_file = 'k/%s%s.kernel' % (base,alt)
	tfiles.add(out_file)
	if '119.004099' in cys_mod:
		cys_mod.remove('119.004099')
	if len(cys_mod) > 0:
		os.system('python3 ktools/modify.py k/%s.kernel %s %s@U -f' % (base,out_file,','.join(cys_mod)))

base += alt
alt = '.pQ'
out_file = 'k/%s%s.kernel' % (base,alt)
tfiles.add(out_file)
os.system('python3 ktools/modify.py k/%s.kernel %s -17.026549@[Q -v' % (base,out_file))

if 57.021464 in cys_mod:
	base += alt
	alt = '.pC'
	out_file = 'k/%s%s.kernel' % (base,alt)
	tfiles.add(out_file)
	os.system('python3 ktools/modify.py k/%s.kernel %s -17.026549@[C -v' % (base,out_file))


base += alt
alt = '.MO'
out_file = 'k/%s%s.kernel' % (base,alt)
tfiles.add(out_file)
os.system('python3 ktools/modify.py k/%s.kernel %s 15.994915@2M -v' % (base,out_file))

base += alt
alt = '.WO'
out_file = 'k/%s%s.kernel' % (base,alt)
tfiles.add(out_file)
os.system('python3 ktools/modify.py k/%s.kernel %s 31.98983@1W -v' % (base,out_file))

base += alt
alt = '.ac'
out_file = 'k/%s%s.kernel' % (base,alt)
os.system('python3 ktools/modify.py k/%s.kernel %s 42.010565@^ -v' % (base,out_file))

if use_special:
	tfiles.add(out_file)
	base += alt
	alt = '.sp'
	out_file = 'k/%s%s.kernel' % (base,alt)
	os.system('python3 ktools/special_modify.py k/%s.kernel %s' % (base,out_file))

os.rename(out_file,output)

for t in tfiles:
	os.remove(t)

