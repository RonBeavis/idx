#!c:/python36/python.exe

#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#

import cgi,cgitb
import os
os.environ[ 'HOME' ] = 'c:/temp'
import time
import sys
import re
import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator

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
		font-size: 14pt;
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
<script type="text/javascript">
function next()	{
	var d = 1;
	var dobj = 1;
	var dlimit = 1;
	for(d = 1; d < 20; d++)	{
		obj = document.getElementById(d);
		if(obj)	{
			if(obj.style.display != "none")	{
				dobj = d;
			}
		}
		else	{
			dlimit = d;
			break;
		}
	}
	if(dobj + 1 < dlimit)	{
		dobj += 1;
	}
	else	{
		dobj = 1;
	}
	for(d = 1; d < dlimit; d++)	{
		obj = document.getElementById(d);
		if(obj)	{
			if(d == dobj)	{
				obj.style.display = "block";
				obj.style.visibility = "visible";
				document.getElementById("pos").innerHTML = d.toString() + "/" + (dlimit-1).toString();
			}
			else	{
				obj.style.display = "none";
				obj.style.visibility = "hidden";
			}
		}
	}
}
function back()	{
	var d = 1;
	var dobj = 1;
	var dlimit = 1;
	for(d = 1; d < 20; d++)	{
		obj = document.getElementById(d);
		if(obj)	{
			if(obj.style.display != "none")	{
				dobj = d;
			}
		}
		else	{
			dlimit = d;
			break;
		}
	}
	if(dobj - 1 < 1)	{
		dobj = dlimit - 1;
	}
	else	{
		dobj -= 1;
	}
	for(d = 1; d < dlimit; d++)	{
		obj = document.getElementById(d);
		if(obj)	{
			if(d == dobj)	{
				obj.style.display = "block";
				obj.style.visibility = "visible";
				document.getElementById("pos").innerHTML = d.toString() + "/" + (dlimit-1).toString();
			}
			else	{
				obj.style.display = "none";
				obj.style.visibility = "hidden";
			}
		}
	}
}
</script>
</head>
<body>
		'''

def link_bar(_fn):
	htm = '<div class="linkBar" style="color: white">'
	htm +='&nbsp;|&nbsp;<a href="/i/index.html" target="_blank">search</a>'
	htm += '&nbsp;|&nbsp;<a href="/a/look.py?fn=%s">table</a>' % (_fn)
	htm += '&nbsp;|&nbsp;<a href="/a/aaa.py?fn=%s">aaa</a>' % (_fn)
	htm += '&nbsp;|&nbsp;<a href="/a/images.py?fn=%s">graphs</a>' % (_fn)
	htm += '&nbsp;|&nbsp;'
	htm += '</div><br /><br />'
	return htm


def end_page():
	return '''</body>\n</html>\n'''

def draw_ppm(_h,_ls,_r,_id,_dims):
	ppm = {}
	ippm = _h.index('ppm')
	for i in range(-20,21):
		ppm[i] = 0
	for l in _ls:
		ppm[round(float(l[ippm])+0.5)] += 1

	xs = []
	ys = []
	for i in range(-20,21):
		xs.append(i)
		ys.append(ppm[i])
	mpl.style.use('seaborn-notebook')
	ms = 6
	plt.bar(xs,ys,color=(0.25,0,1,.8))
	plt.yscale('linear')
	plt.ylabel('PSMs')
	plt.xlabel('ppm')
	plt.legend(loc='best')
	plt.title('Parent ion mass accuracy')
	ax = plt.gca()
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
	ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
	ax.xaxis.set_major_locator(MaxNLocator(integer=True))
	fig = plt.gcf()
	fig.set_size_inches(_dims[0], _dims[1])
	fig.savefig('..\\temp\\%s_ppms.png' % (_r), dpi=400, bbox_inches='tight')
	plt.close()
	return '<div id="%i" style="display: inline;"><img src="/temp/%s_ppms.png" width="800"/></div>\n' % (_id,_r)

def draw_length(_h,_ls,_r,_id,_dims):
	ls = {}
	sl = 0;
	seq = set()
	iseq = _h.index('Sequence')
	for l in _ls:
		s = l[iseq]
		if s in seq:
			pass
		else:
			seq.add(s)
			sl = len(s)
			if sl in ls:
				ls[sl] += 1
			else:
				ls[sl] = 1
	xs = []
	ys = []
	for i in sorted(ls):
		xs.append(i)
		ys.append(ls[i])
	mpl.style.use('seaborn-notebook')
	ms = 6
	plt.bar(xs,ys,color=(1,0,0.25,.8))
	plt.yscale('linear')
	plt.ylabel('peptides')
	plt.xlabel('length')
	plt.legend(loc='best')
	plt.title('Peptide length distribution')
	plt.xlim(left=0)
	ax = plt.gca()
	box = ax.get_position()
	ax.xaxis.set_major_locator(MaxNLocator(integer=True))
	ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
	ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
	fig = plt.gcf()
	fig.set_size_inches(_dims[0], _dims[1])
	fig.savefig('..\\temp\\%s_len.png' % (_r), dpi=400, bbox_inches='tight')
	plt.close()
	return '<div id="%i" style="display: none;"><img src="/temp/%s_len.png" width="800"/></div>\n' % (_id,_r)

def draw_ions(_h,_ls,_r,_id,_dims):
	ics = {}
	sl = 0;
	iIC = _h.index('IC')
	for l in _ls:
		s = int(l[iIC])
		if s in ics:
			ics[s] += 1
		else:
			ics[s] = 1
		
	xs = []
	ys = []
	for i in sorted(ics):
		xs.append(i)
		ys.append(ics[i])
	mpl.style.use('seaborn-notebook')
	ms = 6
	plt.bar(xs,ys,color=(0.5,0.5,0.25,.8))
	plt.yscale('linear')
	plt.ylabel('PSMs')
	plt.xlabel('ion count')
	plt.legend(loc='best')
	plt.title('Ion count distribution')
	plt.xlim(left=0)
	ax = plt.gca()
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
	ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
	ax.xaxis.set_major_locator(MaxNLocator(integer=True))
	fig = plt.gcf()
	fig.set_size_inches(_dims[0], _dims[1])
	fig.savefig('..\\temp\\%s_ics.png' % (_r), dpi=400, bbox_inches='tight')
	plt.close()
	return '<div id="%i" style="display: none;"><img src="/temp/%s_ics.png" width="800"/></div>\n' % (_id,_r)

def draw_ns(_h,_ls,_r,_id,_dims):
	ics = {}
	sl = 0;
	ilf = _h.index('log(f)')
	for l in _ls:
		s = round(float(l[ilf]))
		if s in ics:
			ics[s] += 1
		else:
			ics[s] = 1
		
	xs = []
	ys = []
	for i in sorted(ics):
		xs.append(i)
		ys.append(ics[i])
	mpl.style.use('seaborn-notebook')
	ms = 6
	plt.bar(xs,ys,color=(0.25,0.5,1,.8))
	plt.ylabel('PSMs')
	plt.xlabel('log(f)')
	plt.title('Historical frequency distribution')
	ax = plt.gca()
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
	ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
	ax.xaxis.set_major_locator(MaxNLocator(integer=True))
	fig = plt.gcf()
	fig.set_size_inches(_dims[0], _dims[1])
	fig.savefig('..\\temp\\%s_ns.png' % (_r), dpi=400, bbox_inches='tight')
	plt.close()
	return '<div id="%i" style="display: none;"><img src="/temp/%s_ns.png" width="800"/></div>\n' % (_id,_r)

def draw_ionf(_h,_ls,_r,_id,_dims):
	ics = {}
	sl = 0;
	ilf = _h.index('RI')
	m = 10.0
	for l in _ls:
		s = round(float(l[ilf])*m)
		if s in ics:
			ics[s] += 1
		else:
			ics[s] = 1
		
	xs = []
	ys = []
	for i in sorted(ics):
		xs.append(float(i)/m)
		ys.append(ics[i])
	mpl.style.use('seaborn-notebook')
	ms = 6
	plt.bar(xs,ys,color=(0.25,0.75,0.5,.8),width=0.08)
	plt.ylabel('PSMs')
	plt.xlabel('RI')
	plt.title('Relative intensity')
	ax = plt.gca()
	ax.set_xlim(0,1)
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
	ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
	fig = plt.gcf()
	fig.set_size_inches(_dims[0], _dims[1])
	fig.savefig('..\\temp\\%s_ri.png' % (_r), dpi=400, bbox_inches='tight')
	plt.close()
	return '<div id="%i" style="display: none;"><img src="/temp/%s_ri.png" width="800"/></div>\n' % (_id,_r)

def draw_a0(_h,_ls,_r,_id,_dims):
	htm = ''
	ics = {'A0':0,'A1':0,'A2':0}
	sl = 0;
	ilf = _h.index('Delta')
	for l in _ls:
		s = float(l[ilf])
		if s < 0.5 :
			ics['A0'] += 1
		elif s > 0.5 and s < 1.5:
			ics['A1'] += 1
		elif s > 1.5:
			ics['A2'] += 1
		
	xs = [0,1,2]
	tot = ics['A0']+ics['A1']+ics['A2']
	ys = [ics['A0']/tot,ics['A1']/tot,ics['A2']/tot]
	mpl.style.use('seaborn-notebook')
	ms = 6
	plt.bar(xs,ys,color=(1.0,0.75,0.25,.8),width=0.08)
	plt.ylabel('f(PSMs)')
	plt.xlabel('peak')
	plt.title('A0, A1, A2 assignments')
	ax = plt.gca()
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
	ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
	ax.set_xticks([0,1,2])
	ax.set_xticklabels(['A0','A1','A2'])
	plt.xlim(left=-1,right=3)
	plt.ylim(bottom=0,top=1.05)
	fig = plt.gcf()
	fig.set_size_inches(_dims[0], _dims[1])
	fig.savefig('..\\temp\\%s_a0.png' % (_r), dpi=400, bbox_inches='tight')
	plt.close()
	return '<div id="%i" style="display: none;"><img src="/temp/%s_a0.png" width="800"/></div>\n' % (_id,_r)

def draw_zs(_h,_ls,_r,_id,_dims):
	ics = {1:0,2:0,3:0,4:0,5:0}
	sl = 0;
	ilf = _h.index('z')
	for l in _ls:
		z = int(l[ilf])
		if z in ics:
			ics[z] += 1
		else:
			ics[z] += 1
		
	xs = []
	ys = []
	for i in sorted(ics):
		xs.append(i)
		ys.append(ics[i])
	mpl.style.use('seaborn-notebook')
	ms = 6
	plt.bar(xs,ys,color=(1.0,0.25,0.75,.8),width=0.08)
	plt.ylabel('PSMs')
	plt.xlabel('z')
	plt.title('Parent ion charge distribution')
	ax = plt.gca()
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
	ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
	ax.xaxis.set_major_locator(MaxNLocator(integer=True))
	fig = plt.gcf()
	fig.set_size_inches(_dims[0], _dims[1])
	fig.savefig('..\\temp\\%s_zs.png' % (_r), dpi=400, bbox_inches='tight')
	plt.close()
	return '<div id="%i" style="display: none;"><img src="/temp/%s_zs.png" width="800"/></div>\n' % (_id,_r)

def draw_scans(_h,_ls,_r,_id,_dims):
	ss = {}
	scan = 0
	sub = 0;
	iscan = _h.index('Scan')
	isub = _h.index('Sub')
	for l in _ls:
		if l[isub] != '1':
			continue
		scan = round(float(l[iscan])/100)
		if scan in ss:
			ss[scan] += 1
		else:
			ss[scan] = 1
	xs = []
	ys = []
	psms = 0;
	bins = 1
	for i in sorted(ss):
		xs.append(i*100 + 50)
		bins = 1
		psms = 0
		for r in range(i-1,i+2):
			if r in ss:
				psms += ss[r]
				bins += 1
		ys.append(psms/(bins*100.0))
	mpl.style.use('seaborn-notebook')
	ms = 0
	plt.plot(xs,ys,color=(0,0.5,0.25,.8),markersize=ms,marker='o',linestyle='solid',label='')
	plt.yscale('linear')
	plt.ylabel('frequency')
	plt.xlabel('scan')
	plt.legend(loc='best')
	plt.title('Scan id frequency')
	plt.xlim(left=0)
	ax = plt.gca()
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
	ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
	fig = plt.gcf()
	fig.set_size_inches(_dims[0], _dims[1])
	fig.savefig('..\\temp\\%s_freq.png' % (_r), dpi=400, bbox_inches='tight')
	plt.close()
	return '<div id="%i" style="display: none;"><img src="/temp/%s_freq.png" width="800"/></div>\n' % (_id,_r)
	
def draw_hscan(_h,_ls,_r,_id,_dims):
	htm = ''
	ss = {}
	scan = 0
	sub = 0;
	iscan = _h.index('Scan')
	isub = _h.index('Sub')
	iseq = _h.index('Sequence')
	hydrophilic = 'STEDQNKRPYH'
	hydrophobic = 'LIVMFWVA'
	plus = 'KRH'
	minus = 'ED'
	philic_ratio = {}
	phobic_ratio = {}
	pcharge = {}
	mcharge = {}
	lratio = {}
	for l in _ls:
		if l[isub] != '1':
			continue
		scan = int(l[iscan])
		seq = l[iseq]
		philic = 0
		phobic = 0
		p = 0
		m = 0
		for i in range(0,len(seq)):
			s = seq[i]
			if hydrophilic.find(s) != -1:
				philic += 1
			if hydrophobic.find(s) != -1:
				phobic += 1
			if plus.find(s) != -1:
				p += 1
			if minus.find(s) != -1:
				m += 1
		philic_ratio[scan] = philic/len(seq)
		phobic_ratio[scan] = phobic/len(seq)
		pcharge[scan] = p/len(seq)
		mcharge[scan] = m/len(seq)
		lratio[scan] = len(seq)
	xs1 = []
	ys1 = []
	save = 0
	rave = 0
	bin = 0
	bin_size = 100
	
	for i in sorted(phobic_ratio):
		save += i
		rave += phobic_ratio[i]
		bin += 1
		if bin == bin_size:
			xs1.append(save/bin_size)
			ys1.append(rave/bin_size)
			bin = 0
			save = 0
			rave = 0
	xs2 = []
	ys2 = []
	save = 0
	rave = 0
	bin = 0
	for i in sorted(philic_ratio):
		save += i
		rave += philic_ratio[i]
		bin += 1
		if bin == bin_size:
			xs2.append(save/bin_size)
			ys2.append(rave/bin_size)
			bin = 0
			save = 0
			rave = 0
	xs3 = []
	ys3 = []
	save = 0
	rave = 0
	bin = 0
	
	for i in sorted(pcharge):
		save += i
		rave += pcharge[i]
		bin += 1
		if bin == bin_size:
			xs3.append(save/bin_size)
			ys3.append(rave/bin_size)
			bin = 0
			save = 0
			rave = 0
	xs4 = []
	ys4 = []
	save = 0
	rave = 0
	bin = 0
	
	for i in sorted(mcharge):
		save += i
		rave += mcharge[i]
		bin += 1
		if bin == bin_size:
			xs4.append(save/bin_size)
			ys4.append(rave/bin_size)
			bin = 0
			save = 0
			rave = 0
	xs5 = []
	ys5 = []
	save = 0
	rave = 0
	bin = 0
	ml = 0
	for i in sorted(lratio):
		save += i
		rave += lratio[i]
		bin += 1
		if bin == bin_size:
			if rave/bin_size > ml:
				ml = rave/bin_size
			bin = 0
			save = 0
			rave = 0
	save = 0;
	rave = 0;
	bin = 0;
	for i in sorted(lratio):
		save += i
		rave += lratio[i]
		bin += 1
		if bin == bin_size:
			xs5.append(save/(bin_size))
			ys5.append(rave/(ml*bin_size))
			bin = 0
			save = 0
			rave = 0
	mpl.style.use('seaborn-notebook')
	ms = 2
	plt.scatter(xs5,ys5,color=(1.0,0.7,0.3,.8),label="length/maximum")
	plt.scatter(xs1,ys1,color=(0,0.7,0.25,.8),label="hydrophobic")
	plt.scatter(xs2,ys2,color=(0.1,0.0,0.9,.8),label="hydrophilic")
	plt.scatter(xs3,ys3,color=(0.9,0.0,0.1,.8),label="positive")
	plt.scatter(xs4,ys4,color=(0.0,0.0,0.0,.8),label="negative")
	plt.yscale('linear')
	plt.ylabel('fraction of residues')
	plt.xlabel('scan')
	plt.title('Scan vs sequence correlations')
	plt.xlim(left=0)
	plt.ylim(bottom=0,top=1.05)
	ax = plt.gca()
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
	ax.legend(loc='upper left')
	fig = plt.gcf()
	fig.set_size_inches(_dims[0], _dims[1])
	fig.savefig('..\\temp\\%s_corr.png' % (_r), dpi=400, bbox_inches='tight')
	plt.close()
	return '<div id="%i" style="display: none;"><img src="/temp/%s_corr.png" width="800"/></div>\n' % (_id,_r)

form = cgi.FieldStorage()
print('Content-type: text/html\n\n')
clean_temp('D:/somecrazyblogger-org/temp')
html = start_page()
display_ratio = 1/2.0
try:
	display_ratio = 1.0/float(form['r'].value)
except:
	pass
try:
	rname = form['fn'].value
except:
	print(html + 'Error: requires a value for "fn"')
	print(end_page())
	exit()
fname = 'f:\\cidx\\o\\' + rname
html += link_bar(rname)
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
id = 1
html += '<table border="0">\n<tr><td>\n'

html += draw_ppm(headers,lines,rname,id,(10,10*display_ratio))
id += 1

html += draw_length(headers,lines,rname,id,(10,10*display_ratio))
id += 1

html += draw_scans(headers,lines,rname,id,(10,10*display_ratio))
id += 1

html += draw_ions(headers,lines,rname,id,(10,10*display_ratio))
id += 1

html += draw_ns(headers,lines,rname,id,(10,10*display_ratio))
id += 1

html += draw_ionf(headers,lines,rname,id,(10,10*display_ratio))
id += 1

html += draw_hscan(headers,lines,rname,id,(10,10*display_ratio))
id += 1

html += draw_a0(headers,lines,rname,id,(10,10*display_ratio))
id += 1

html += draw_zs(headers,lines,rname,id,(10,10*display_ratio))
id += 1

html += '\n</td>\n</tr>\n<tr><td align="center" valign="top">'
html += '''<a href="javascript: back();">
		<img src="/i/left.png" height="14" /></a>&nbsp;&nbsp;
		<span id="pos">1/%i</span>&nbsp;&nbsp;
		<a href="javascript: next();">
		<img src="/i/right.png" height="14" /></a></td>''' % (id-1)
html += '</tr>\n</table>\n'

html += end_page()
print(html)
