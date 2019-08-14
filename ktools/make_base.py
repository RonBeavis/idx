import os

output = 'k/mm.KR.kernel'
base = 'mm.KR'

alt = '.ia'
out_file = 'k/%s%s.kernel' % (base,alt)
tfiles = set({out_file})
os.system('python3 ktools/modify.py ../kernels/%s.kernel %s 57.021464@C -f' % (base,out_file))

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
os.system('python3 ktools/modify.py k/%s.kernel %s 15.994915@M -v' % (base,out_file))

base += alt
alt = '.ac'
out_file = 'k/%s%s.kernel' % (base,alt)
os.system('python3 ktools/modify.py k/%s.kernel %s 42.010565@^ -v' % (base,out_file))

os.rename(out_file,output)

for t in tfiles:
	os.remove(t)

