#!c:/python36/python.exe

#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
# Part of the idx web wrapper
#

import cgi,cgitb
import os
import time

form = cgi.FieldStorage()
fn = form['fn'].value
ifile = 'f:\idx\o\%s' % (fn)
if os.path.isfile(ifile):
	f = open(ifile,'r')
else:
	print('Content-type: text/plain; charset=utf8;\n')
	print('File "%s" no longer on server.\nOutput files are deleted on retrieval.' % (fn))
	exit()

print('Content-type: application/vnd.ms-excel; charset=utf8;\nContent-Disposition: attachment; filename=%s;\n' % (fn))
for l in f:
	print(l.rstrip())
f.close()
try:
	os.remove(ifile)
	fs = os.listdir('f:\idx\o')
	t = time.time()
	for f in fs:
		cfile = 'f:\idx\o\%s' % (f)
		if os.path.isfile(cfile) and cfile.find('.tsv') != -1:
			if t - os.path.getctime(cfile) > 3600*6:
				os.remove(cfile)
	fs = os.listdir('f:\idx\s')
	for f in fs:
		cfile = 'f:\idx\s\%s' % (f)
		if os.path.isfile(cfile) and cfile.find('.mgf') != -1:
			if t - os.path.getctime(cfile) > 3600*6:
				os.remove(cfile)

except:
	pass
