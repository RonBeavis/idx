#!c:/python36/python.exe

import cgi,cgitb
import os
import sys
import subprocess
import time
import sqlite3
import re

cgitb.enable()

def start_page():
	print('''<!DOCTYPE html>
<html lang="en" class="no-js">
<head>
<meta http-equiv="X-UA-Compatible" content="IE=edge">
<meta charset="utf-8">
<title>idX MS/MS to PSM processor</title>
<meta name="viewport" content="width=device-width,initial-scale=1" />
<style>
	body {
		font-family: sans-serif;
	}
	a {
		color: #369;
	}
	div.display	{
		background-color: #ddeedd;
		padding: 10px;
		position: absolute;
		bottom: 10px;
		top: 10px;
		right: 10px;
		left: 200px;
	}
	p {
		display: block;
		margin-top: 0.5em;
		margin-bottom: 0.5em;
		margin-left: 5 em;
		margin-right: auto;
		padding-left: -200px;
	}
	.d1 {
		display: block;
		margin-top: 0.5em;
		margin-bottom: 0.5em;
		margin-left: 5 em;
		margin-right: auto;
		padding-left: -200px;
	}
	.d2 {
		display: block;
		margin-top: 0.5em;
		margin-bottom: 0.5em;
		margin-left: 5 em;
		margin-right: auto;
		padding-left: -180px;
		font-style: italic;
	}
	hr.top {
		width: 50%;
		display: block;
		margin-top: 0.5em;
		margin-bottom: 0.5em;
		margin-left: 0;
		margin-right: auto;
		border-style: inset;
		border-width: 1px;
	}
	hr.nottop {
		width: 50%;
		display: block;
		margin-top: 0.5em;
		margin-bottom: 0.5em;
		margin-left: 200px;
		margin-right: auto;
		border-style: inset;
		border-width: 1px;
	}
	img.floatLeft {
		float: left; 
		margin-top: 4px;
		margin-bottom: 4px;
		margin-left: 4px;
		margin-right: 10px;
		width : 200 px;
	}
	div.floatLeft {
		float: left; 
		margin-top: 10px;
		margin-bottom: 2px;
		margin-left: 2px;
		margin-right: 0px;
		width : 200 px;
	}
}
</style>
</head>
<body>
		''')

def get_mod_text(_t):
	if _t == 'hs-plain':
		return '''Checks for:<ol><li>IAA (C)</li><li>acetyl (N-T)</li><li>acetyl (K)</li><li>phosphoryl (S)</li><li>phosphoryl (T)</li><li>
		phosphoryl (Y)</li><li>ubiquitinyl (K)</li><li>deamidated (N)</li><li>pyro (n-t Q)</li><li>pyro (n-t C)</li><li>oxidized (M)</li><li>hydroxyproline</li>
		<li>hydroxylysine</li><li>dimethyl (R)</li><li>protein n-t</li><li>er signal</li><li>mito signal</li></ol>'''
	if _t == 'hs-tmt':
		return '''Checks for:<ol><li>TMT6+ (n-t)</li><li>TMT6+ (K)</li><li>IAA (C)</li><li>acetyl (N-T)</li><li>acetyl (K)</li><li>phosphoryl (S)</li>
		<li>phosphoryl (T)</li><li>phosphoryl (Y)</li><li>ubiquitinyl (K)</li><li>deamidated (N)</li><li>pyro (n-t Q)</li><li>pyro (n-t C)</li><li>oxidized (M)</li>
		<li>hydroxyproline</li><li>hydroxylysine</li><li>dimethyl (R)</li><li>protein n-t</li><li>er signal</li><li>mito signal</li></ol>'''
	if _t == 'hs-hla1':
		return '''Checks for:<ol><li>acetyl (N-T)</li><li>cystine (C)</li><li>deamidated (N)</li><li>pyro (n-t Q)</li><li>oxidized (M)</li>
		<li>hydroxyproline</li><li>hydroxylysine</li><li>dimethyl (R)</li></ol>'''
	if _t == 'hs-hla2':
		return '''Checks for:<ol><li>acetyl (N-T)</li><li>cystine (C)</li><li>deamidated (N)</li><li>pyro (n-t Q)</li><li>oxidized (M)</li>
		<li>hydroxyproline</li><li>hydroxylysine</li><li>dimethyl (R)</li></ol>'''
	if _t == 'hs-hla1-2':
		return '''Checks for:<ol><li>acetyl (N-T)</li><li>cystine (C)</li><li>deamidated (N)</li><li>pyro (n-t Q)</li><li>oxidized (M)</li>
		<li>hydroxyproline</li><li>hydroxylysine</li><li>dimethyl (R)</li></ol>'''
	return '''<br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br />'''
	
def end_page():
	print('''</body>\n</html>\n''')
	
kernels = {'hs-plain': 'f:\idx\k\hs_plain.KR.kernel', 'hs-tmt': 'f:\idx\k\hs_tmt.KR.kernel',
			'hs-hla1': 'f:\idx\k\hla_class1.kernel', 'hs-hla2': 'f:\idx\k\hla_class2.kernel',
			'hs-hla1-2': 'f:\idx\k\hla_class1+2.kernel'}
	
form = cgi.FieldStorage()
print('Content-type: text/html\n\n')
start_page()
print('<div class="floatLeft"><a href="/idx"><img src="/idx/idx.png" border="0" /></a>')
mod_text = get_mod_text(form['etype'].value)
print('<br />%s</div>' % (mod_text))
print('<div class="display">')
print('<p>Starting idX session ...</p>')
sys.stdout.flush()
mgf = 'f:\idx\s\%s' % (form['nn'].value)
display_mgf = form['fn'].value
print('<p>Upload file: %s</p>' % (display_mgf))

ofn = form['nn'].value + '.tsv'
ofpath = 'f:\idx\o\%s' % (ofn)
idxpath = 'f:\idx\idX.py'
pythonpath = 'c:\python36\python.exe'
if mgf.find('.raw') != -1:
	param = [pythonpath, 'f:\idx\convert.py',form['nn'].value]
	param = [pythonpath, 'f:\idx\convert.py','AD1_01.raw']
	param = ['C:/Program Files/ProteoWizard/ProteoWizard 3.0.18172.8eb65f19f/msconvert.exe',mgf,'--outdir','f:\idx\s','--mgf','--filter','peakPicking true 2','--filter','msLevel 2',]
	x = subprocess.Popen(param, stdout=subprocess.PIPE)
	print('<hr class="top"/>')
	print('<div class="d1" id="2">processing .raw file (sec) = 0</div>')
	sys.stdout.flush()
	raw_ok = False
	secs = 0
	while x.poll() == None:
		print('<script>document.getElementById("2").innerHTML="processing .raw file (sec) = %i" </script>' % (secs))
		sys.stdout.flush()
		secs += 5
		time.sleep(5)
	print('<hr class="top"/>')
	if secs == 0:
		print('Could not process raw file')
		print('<hr class="top" /><br /><p><a href="javascript: window.history.back();">Run another search</a></p>')
		exit()
	try:
		os.remove(mgf)
	except:
		print('<p>Error deleting "%s"</p>' % (mgf))
	mgf = re.sub('\.raw$','.mgf',mgf)
	ofn = re.sub('\.raw','.mgf',ofn)
	ofpath = re.sub('\.raw','.mgf',ofpath)

etype = form['etype'].value
resolution = form['res'].value
ress = {'high': 20, 'medium': 50, 'low': 400}
ipaddress = cgi.escape(os.environ["REMOTE_ADDR"])
slimit = 15000
try:
	conn = sqlite3.connect('f:\idx\db\idx.db')
	curses = conn.cursor()
except:
	print('SQLITE3 connection failure')
	exit()
curses.execute("select count(*) from session where active=1 and stime > (strftime('%s', 'now')-1800)*1000")
rows = curses.fetchall()
max_threads = 2
outfile = True
if rows[0][0]+1 > max_threads:
	curses.close()
	conn.close()
	print('Too many idX searches current active: try again later')
	exit()
current_sid = None
if os.path.isfile(mgf):
	print('<p>Processed data : %s B</p>' % (format(os.path.getsize(mgf), ',d')))
	print('<!--%s-->' % (mgf))
	print('<p>Total spectra limit: %i</p>' % (slimit))
	dkernel = re.sub(r'.+\\','',kernels[etype])
	dkernel = re.sub(r'\.kernel','',dkernel)
	print('<p>Kernel: %s</p>' % (dkernel))
	if resolution in ress:
		print('<p>Fragment tolerance: %i mDa</p>' % (ress[resolution]))
	param = [pythonpath, '-u', idxpath, mgf, kernels[etype], ofpath ,resolution,'%i' % (slimit)]
	print('<hr class="top"/>')
	x = subprocess.Popen(param, stdout=subprocess.PIPE)
	sql = 'INSERT INTO session(stime,etime,pid,ofile,ifile,etype,userid,active,ip) values(?,?,?,?,?,?,?,?,?)'
	vals = (int(time.time()*1000),None,x.pid,ofn,display_mgf,etype,None,1,ipaddress)
	curses.execute(sql,vals)
	conn.commit()
	curses.execute('select last_insert_rowid()')
	rows = curses.fetchall()
	current_sid = rows[0][0]
	m = x.stdout.readline().decode('utf8')
	t0 = time.time()
	print('<div class="d1" id="100">loading spectra (sec) = 0</div>')
	ok = True
	while ok and len(m) > 0:
		now = time.time()
		print('<script>document.getElementById("100").innerHTML="loading spectra (sec) = %.0f" </script>' % (now-t0))
		sys.stdout.flush()
		m = x.stdout.readline().decode('utf8')
		if m.find('spectra = ') != -1:
			print('<div class="d2">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;%s</div>' % (m[m.find('spectra = '):]))
		if m.find('exiting:') != -1:
			print('<div class="d2">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;%s</div>' % (m[m.find('exiting: '):]))
			ok = False
		if m.find('load & index kernel') != -1:
			break;
	t0 = time.time()
	if ok:
		print('<div class="d1" id="200">loading kernels (sec) = 0</div>')
	while ok and len(m) > 0:
		now = time.time()
		print('<script>document.getElementById("200").innerHTML="loading kernels (sec) = %.0f" </script>' % (now-t0))
		sys.stdout.flush()
		m = x.stdout.readline().decode('utf8')
		if m.find('kernels = ') != -1:
			print('<div class="d2">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;%s</div>' % (m[m.find('kernels = '):]))
		if m.find('perform ids') != -1:
			break;
	t0 = time.time()
	if ok:
		print('<div class="d1" id="300">perform ids (sec) = 0</div>')
	while ok and len(m) > 0:
		now = time.time()
		print('<script>document.getElementById("300").innerHTML="perform ids (sec) = %.0f" </script>' % (now-t0))
		sys.stdout.flush()
		m = x.stdout.readline().decode('utf8')
		if m.find('create report') != -1:
			break;
	t0 = time.time()
	if ok:
		print('<div class="d1" id="400">create report (sec) = 0</div>')
	while ok and len(m) > 0:
		now = time.time()
		print('<script>document.getElementById("400").innerHTML="create report (sec) = %.0f" </script>' % (now-t0))
		sys.stdout.flush()
		if m.find('lines = 0') != -1:
			print('<div class="d2">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;%s</div>' % (m[m.find('lines = '):]))
			ok = False
			outfile = False
			break
		elif m.find('lines = ') != -1:
			print('<div class="d2">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;%s</div>' % (m[m.find('lines = '):]))
		if m.find('ble = ') != -1:
			print('<div class="d2">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;%s</div>' % (m[m.find('ble = '):]))
		if m.find('ppm window = ') != -1:
			print('<div class="d2">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;%s</div>' % (m[m.find('ppm window = '):]))
		m = x.stdout.readline().decode('utf8')
	
	try:
		os.remove(mgf)
		print('<div class="d2" id="500"><i>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Data file "%s" cleared</i></div>' % (display_mgf))
	except:
		print('<div class="d2" id="500"><i>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Data file "%s" could not be cleared</i></div>' % (display_mgf))
else:
	ok = False
	print('<p>Data file "%s" not on server: process could not start</p>' % (mgf))
print('<hr class="top" />')
print('<p>idX session complete: </p>')
url = '/thegpm-cgi/results.py?fn=%s' % (ofn)
if ok:
	print('<p>&nbsp;&nbsp;&nbsp;your results: <a href="%s" target="_blank">%s</a></p>' % (url,ofn))
else:
	if outfile:
		print('<p>&nbsp;&nbsp;&nbsp;no result file generated because of an error.</p>')
	else:
		print('<p>&nbsp;&nbsp;&nbsp;no PSMs assigned.</p>')
print('<hr class="top" /><br /><p><a href="javascript: window.history.back();">Run another search</a></p>')
print('</div>')

if current_sid:
	sql = 'UPDATE session set etime=?,active=? where sid=?'
	vals = (int(time.time()*1000),0,current_sid)
	curses.execute(sql,vals)
	conn.commit()

curses.close()
conn.close()
end_page()

