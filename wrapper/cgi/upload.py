#!c:/python36/python.exe

#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
# Part of the idx web wrapper
#

import cgi,cgitb
from shutil import copyfile

cgitb.enable()
form = cgi.FieldStorage()
print('Content-type: text/plain\n')
try:
	fileitem = form['filename']
	nn = form['newname'].value
	f = open('f:/idx/s/%s' % (nn),'wb')
	f.write(fileitem.file.read())
	f.close()
	print('ok')
except:
	print('error')

	
