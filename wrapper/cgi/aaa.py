#!c:/python36/python.exe

import cgi,cgitb
import os
os.environ[ 'HOME' ] = 'c:/temp'
import sys
import re
import time
import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib as mpl

cgitb.enable()

def clean_temp(_r):
	fs = [os.path.join(_r,f) for f in os.listdir(_r) if os.path.isfile(os.path.join(_r,f))]
	now = time.time()
	for f in fs:
		if os.stat(f).st_mtime < now - 3600:
			try:
				os.remove(f)
			except:
				print(f)

def start_page():
	return '''<!DOCTYPE html>
<html lang="en" class="no-js">
<head>
<meta http-equiv="X-UA-Compatible" content="IE=edge">
<meta charset="utf-8">
<title>idX MS/MS to PSM processor</title>
<meta name="viewport" content="width=device-width,initial-scale=1" />
<link rel="stylesheet" href="/css/scb.css">
<style>
	body {
		font-family: Roboto;
		font-size: 12pt;
	}
	a {
		color: #0000FF;
		text-decoration: none;
	}
	td	{
		white-space: nowrap;
		overflow: hidden;
		text-align: center;
	}
	td.left	{
		white-space: nowrap;
		overflow: hidden;
		color: white;
		background-color: black;
		text-align: right;
		width: 120px;
	}
	td.header	{
		white-space: nowrap;
		overflow: hidden;
		text-align: center;
		color: white;
		background-color: black;
	}
	.linkBar {
		font-size: 11pt;
		color: rgb(102,102,102);
		position: fixed;
		z-index: 100;
		text-align: center;
		top: 0px;
		left: 0px;
		margin: 0px;
		width: 100%;
		height: 24px;
		line-height: 24px;
		padding-top: 0px;
		background: none repeat scroll 0% 0% rgb(0, 0,150);
		border-top: 1px solid rgb(50,0,0);
		box-shadow: 0px 0px 10px rgb(0, 0, 0, 0.5);
	}
	.linkBar a {color: white;font-family: Roboto,Arial, Helvetica, sans-serif;text-decoration: none;font-size: 11pt;}
	.linkBar a:hover {color: white;font-family: Roboto,Arial, Helvetica, sans-serif;text-decoration: underline;font-size: 11pt;}
	.linkBar a:visited {color: white;font-family: Roboto,Arial, Helvetica, sans-serif;text-decoration: none;font-size: 11pt;}
	.linkBar a:active {color: white;font-family: Roboto,Arial, Helvetica, sans-serif;text-decoration: none;font-size: 11pt;}
</style>
</head>
<body>
		'''

	
def end_page():
	return '''</body>\n</html>\n'''

def draw_graphs(_h,_ls,_r):
	htm = '<br /> <hr style="width: 1000px; margin-left: 5px;" />'
	aas = ['A','C','D','E','F','G','H','I','K','L','M','N','O','P','Q','R','S','T','U','V','W','Y','[',']']
	ipre = _h.index('Pre')
	ipost = _h.index('Post')
	iseq = _h.index('Sequence')
	pre = {}
	post = {}
	aa1 = {}
	aaN = {}
	pep = {}
	for a in aas:
		pre[a] = 0
		post[a] = 0
		aa1[a] = 0
		aaN[a] = 0
		pep[a] = 0
	for l in _ls:
		pre[l[ipre]] += 1
		post[l[ipost]] += 1
		seq = l[iseq]
		aa1[seq[0]] += 1
		aaN[seq[-1]] += 1
		for a in range(0,len(seq)):
			pep[seq[a]] += 1
	tot = sum(pre.values())
	ys = []
	xs = []
	for i,a in enumerate(aas):
		xs.append(i)
		ys.append(100*pre[a]/tot)
	mpl.style.use('seaborn-notebook')
	ms = 6
	plt.bar(xs,ys,color=(0.8,0.25,0.25,.8),width=0.5)
	plt.ylabel('AAA (percent)')
	plt.xlabel('residue')
	plt.title('pre residue AAA')
	ax = plt.gca()
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
	ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
	ax.set_xticks(xs)
	ax.set_xticklabels(aas)
	plt.xlim(left=-1,right=len(aas))
#	plt.ylim(bottom=0,top=1.05)
	fig = plt.gcf()
	_dims = (10,2.5)
	fig.set_size_inches(_dims[0], _dims[1])
	fig.savefig('..\\temp\\%s_pre.png' % (_r), dpi=400, bbox_inches='tight')
	plt.close()
	htm += '<img src="/temp/%s_pre.png" width="1000"/> <br/>\n' % (_r)
	
	tot = sum(aa1.values())
	ys = []
	xs = []
	for i,a in enumerate(aas):
		xs.append(i)
		ys.append(100*aa1[a]/tot)
	mpl.style.use('seaborn-notebook')
	ms = 6
	plt.bar(xs,ys,color=(0.25,0.8,0.25,.8),width=0.5)
	plt.ylabel('AAA (percent)')
	plt.xlabel('residue')
	plt.title('first residue AAA')
	ax = plt.gca()
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
	ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
	ax.set_xticks(xs)
	ax.set_xticklabels(aas)
	plt.xlim(left=-1,right=len(aas))
#	plt.ylim(bottom=0,top=1.05)
	fig = plt.gcf()
	_dims = (10,2.5)
	fig.set_size_inches(_dims[0], _dims[1])
	fig.savefig('..\\temp\\%s_aa1.png' % (_r), dpi=400, bbox_inches='tight')
	plt.close()
	htm += '<img src="/temp/%s_aa1.png" width="1000"/> </br>\n' % (_r)

	tot = sum(aaN.values())
	ys = []
	xs = []
	for i,a in enumerate(aas):
		xs.append(i)
		ys.append(100*aaN[a]/tot)
	mpl.style.use('seaborn-notebook')
	ms = 6
	plt.bar(xs,ys,color=(0.25,0.25,0.8,.8),width=0.5)
	plt.ylabel('AAA (percent)')
	plt.xlabel('residue')
	plt.title('last residue AAA')
	ax = plt.gca()
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
	ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
	ax.set_xticks(xs)
	ax.set_xticklabels(aas)
	plt.xlim(left=-1,right=len(aas))
#	plt.ylim(bottom=0,top=1.05)
	fig = plt.gcf()
	_dims = (10,2.5)
	fig.set_size_inches(_dims[0], _dims[1])
	fig.savefig('..\\temp\\%s_aaN.png' % (_r), dpi=400, bbox_inches='tight')
	plt.close()
	htm += '<img src="/temp/%s_aaN.png" width="1000"/> </br>\n' % (_r)

	tot = sum(post.values())
	ys = []
	xs = []
	for i,a in enumerate(aas):
		xs.append(i)
		ys.append(100*post[a]/tot)
	mpl.style.use('seaborn-notebook')
	ms = 6
	plt.bar(xs,ys,color=(0.7,0.7,0.25,.8),width=0.5)
	plt.ylabel('AAA (percent)')
	plt.xlabel('residue')
	plt.title('post residue AAA')
	ax = plt.gca()
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
	ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
	ax.set_xticks(xs)
	ax.set_xticklabels(aas)
	plt.xlim(left=-1,right=len(aas))
#	plt.ylim(bottom=0,top=1.05)
	fig = plt.gcf()
	_dims = (10,2.5)
	fig.set_size_inches(_dims[0], _dims[1])
	fig.savefig('..\\temp\\%s_post.png' % (_r), dpi=400, bbox_inches='tight')
	plt.close()
	htm += '<img src="/temp/%s_post.png" width="1000"/> </br>\n' % (_r)

	tot = sum(pep.values())
	ys = []
	xs = []
	for i,a in enumerate(aas):
		xs.append(i)
		ys.append(100*pep[a]/tot)
	mpl.style.use('seaborn-notebook')
	ms = 6
	plt.bar(xs,ys,color=(0.7,0.25,0.7,.8),width=0.5)
	plt.ylabel('AAA (percent)')
	plt.xlabel('residue')
	plt.title('PSM residues AAA')
	ax = plt.gca()
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
	ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
	ax.set_xticks(xs)
	ax.set_xticklabels(aas)
	plt.xlim(left=-1,right=len(aas))
#	plt.ylim(bottom=0,top=1.05)
	fig = plt.gcf()
	_dims = (10,2.5)
	fig.set_size_inches(_dims[0], _dims[1])
	fig.savefig('..\\temp\\%s_psm.png' % (_r), dpi=400, bbox_inches='tight')
	plt.close()
	htm += '<img src="/temp/%s_psm.png" width="1000"/> </br>\n' % (_r)

	return htm

def draw_AAA(_h,_ls):
	htm = ''
	htm += '<table cellspacing="0" cellpadding="6" width="1000">'
	aas = ['A','C','D','E','F','G','H','I','K','L','M','N','O','P','Q','R','S','T','U','V','W','Y','[',']']
	ipre = _h.index('Pre')
	ipost = _h.index('Post')
	iseq = _h.index('Sequence')
	pre = {}
	post = {}
	aa1 = {}
	aaN = {}
	pep = {}
	for a in aas:
		pre[a] = 0
		post[a] = 0
		aa1[a] = 0
		aaN[a] = 0
		pep[a] = 0
	for l in _ls:
		pre[l[ipre]] += 1
		post[l[ipost]] += 1
		seq = l[iseq]
		aa1[seq[0]] += 1
		aaN[seq[-1]] += 1
		for a in range(0,len(seq)):
			pep[seq[a]] += 1
	htm += '<tr><td colspan="%i" align="center">AAA analysis of PSM peptide sequences (percent)</td>' % (len(aas)+1)
	htm += '<tr><td class="left">Residue</td>'
	for a in aas:
		htm += '<td class="header">%s</td>' % (a)
	htm += '</tr>'
	htm += '<tr><td class="left">Pre</td>'
	tot = sum(pre.values())
	for a in aas:
		htm += '<td>%.1f</td>' % (100*pre[a]/tot)
	htm += '</tr>'
	htm += '<tr><td class="left">Post</td>'
	tot = sum(post.values())
	for a in aas:
		htm += '<td>%.1f</td>' % (100*post[a]/tot)
	htm += '</tr>'
	htm += '<tr><td class="left">N-terminal</td>'
	tot = sum(aa1.values())
	for a in aas:
		htm += '<td>%.1f</td>' % (100*aa1[a]/tot)
	htm += '</tr>'
	htm += '<tr><td class="left">C-terminal</td>'
	tot = sum(aaN.values())
	for a in aas:
		htm += '<td>%.1f</td>' % (100*aaN[a]/tot)
	htm += '</tr>'
	htm += '<tr><td class="left">Peptide</td>'
	tot = sum(pep.values())
	for a in aas:
		htm += '<td>%.1f</td>' % (100*pep[a]/tot)
	htm += '</tr>'
	htm += '</table>'
	return htm

def draw_mods(_h,_ls):
	htm = ''
	htm += '<table cellspacing="0" cellpadding="6" width="1000">'
	aas = ['A','C','D','E','F','G','H','I','K','L','M','N','O','P','Q','R','S','T','U','V','W','Y','[',']']
	imods = _h.index('Modifications')
	mods = {}
	for l in _ls:
		vs = l[imods].split(';')
		for v in vs:
			if len(v) == 0:
				continue
			ms = re.split('[~#]',v)
			if len(ms) == 0:
				continue
			c = ms[0][0]
			if ms[1] in mods:
				if c in mods[ms[1]]:
					mods[ms[1]][c] += 1
				else:
					mods[ms[1]][c] = 1
			else:
				mods[ms[1]] = {}
				mods[ms[1]][c] = 1
				
	htm += '<tr><td colspan="%i" align="center"><br /><br />Modification analysis of PSM peptide sequences (residues)</td>' % (len(aas)+1)
	htm += '<tr><td class="left">Residue</td>'
	for a in aas:
		htm += '<td class="header">%s</td>' % (a)
	htm += '</tr>'
	for mtype in sorted(mods):
		htm += '<tr><td class="left">%s</td>' % (mtype)
		for a in aas:
			if a in mods[mtype]:
				htm += '<td>%i</td>' % (mods[mtype][a])
			else:
				htm += '<td>-</td>'
			
		htm += '</tr>'
	htm += '</table>'
	return htm

def link_bar(_fn):
	htm = '<div class="linkBar" style="color: white">'
	htm +='&nbsp;|&nbsp;<a href="/i/index.html" target="_blank">search</a>'
	htm += '&nbsp;|&nbsp;<a href="/a/look.py?fn=%s">table</a>' % (_fn)
	htm += '&nbsp;|&nbsp;<a href="/a/aaa.py?fn=%s">aaa</a>' % (_fn)
	htm += '&nbsp;|&nbsp;<a href="/a/images.py?fn=%s">graphs</a>' % (_fn)
	htm += '&nbsp;|&nbsp;'
	htm += '</div><br /><br />'
	return htm

form = cgi.FieldStorage()
print('Content-type: text/html\n\n')
clean_temp('D:/somecrazyblogger-org/temp')
html = start_page()
try:
	rname = form['fn'].value
except:
	print(html + 'Error: requires a value for "fn"')
	print(end_page())
	exit()
html += link_bar(rname)
fname = 'f:\\cidx\\o\\' + rname
try:
	lines = [(l.rstrip('\n')).split('\t') for l in open(fname,'r')]
except:
	print(html + 'Error: "%s" is not available' % (rname))
	print(end_page())
	exit()
if len(lines) < 2:
	print(html + 'Error: "%s" has no data available' % (rname))
	print(end_page())
	exit()
headers = lines[:1][0]
lines = lines[1:]
html += draw_AAA(headers,lines)
html += draw_mods(headers,lines)
html += draw_graphs(headers,lines,rname)
html += end_page()
print(html)