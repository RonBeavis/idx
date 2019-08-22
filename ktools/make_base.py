import os

output = 'k/hs_sp.KR.kernel'
base = 'hs.1'
use_special = True

alt = '.iac'
out_file = 'k/%s%s.kernel' % (base,alt)
tfiles = set({out_file})
os.system('python3 ktools/modify.py %s.kernel %s 57.021464@C -f' % (base,out_file))

base += alt
alt = '.iau'
out_file = 'k/%s%s.kernel' % (base,alt)
tfiles.add(out_file)
os.system('python3 ktools/modify.py k/%s.kernel %s 57.021464@U -f' % (base,out_file))

base += alt
alt = '.pQ'
out_file = 'k/%s%s.kernel' % (base,alt)
tfiles.add(out_file)
os.system('python3 ktools/modify.py k/%s.kernel %s -17.026549@[Q -v' % (base,out_file))

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

