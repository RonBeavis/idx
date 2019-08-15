import json
import sys

f = open(sys.argv[1],'r')
for j in f:
	js = json.loads(j)
	if 'seq' in js and js['seq'] == 'GPSGPQGPGGPPGPK':
		print(js['seq'])
		if 'mods' in js:
			print(js['mods'])
f.close()

	
