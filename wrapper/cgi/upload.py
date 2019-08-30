#!c:/python36/python.exe

import cgi,cgitb
import os
from shutil import copyfile
import time

cgitb.enable()
form = cgi.FieldStorage()
ipaddress = cgi.escape(os.environ["REMOTE_ADDR"])
print('Content-type: text/plain\n')
log = open('f:/cidx/logs/upload.log','a')

fileitem = form['filename']
nn = form['newname'].value
fn = 'f:/cidx/s/%s' % (nn)

try:
	log.write('***\nstart:\tfile=%s\tip=%s\ttime=%f\n' % (fn,ipaddress,time.time()))
	f = open(fn,'wb')
	f.write(fileitem.file.read())
	f.close()
	log.write('ok:\tfile=%s\tip=%s\tsize=%i\ttime=%f\n' % (fn,ipaddress,os.path.getsize(fn),time.time()))
except:
	pass
	log.write('error:\tfile=%s\tip=%s\ttime=%f\n' % (fn,ipaddress,time.time()))
log.close()
	
