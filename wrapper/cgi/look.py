#!c:/python36/python.exe

import cgi,cgitb
import os
import sys
import re

cgitb.enable()

			
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
		font-size: 10pt;
	}
	a {
		color: #FFDDFF;
	}
	input[type=submit] {
		padding:5px 15px; 
		background:#000;
		color: #FFF;
		border:0 none;
		cursor:pointer;
		-webkit-border-radius: 5px;
		border-radius: 5px; 
	}	
	.filter_form	{
		margin: 10px 10px 10px 10px;
		text-align:left;
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

def add_form(_fn,_fil):
	ret = '<div class="filter_form"><form style="display: inline" action="/a/look.py" method="post">\n'
	ret += '<input type="hidden" name="fn" value="%s" />\n' % (_fn)
	ret += '<input type="text" name="fil" value="%s" size="40" placeholder="show only lines containing this text"/>\n' % (_fil)
	ret += '&nbsp;<input type="submit" name="submit" value="filter" />'
	ret += '</form></div>\n'
	return ret
	
def end_page():
	return '''</body>\n</html>\n'''
	
form = cgi.FieldStorage()
print('Content-type: text/html\n\n')
html = start_page()
try:
	rname = form['fn'].value
	fname = 'f:\\cidx\\o\\' + rname
except:
	print(html + 'Error: requires a value for "fn"')
	print(end_page())
	exit()
html += link_bar(rname)
filter = ''
try:
	filter=form['fil'].value
except:
	filter = ''
try:
	if len(filter) == 0:
		lines = [l.rstrip('\n') for l in open(fname,'r')]
	else:
		tfilter = filter
		if filter.find('\\') == -1 and filter.find('^') == -1 and filter.find('$') == -1  and filter.find('*') == -1 and filter.find('+') == -1 and len(re.findall(r'\[.+\]',filter)) == 0:
			tfilter = re.escape(filter)
		lines = [l.rstrip('\n') for l in open(fname,'r') if len(re.findall(r'%s' % tfilter,l,re.IGNORECASE)) != 0 or l.find('Id') == 0]
except:
	print(html + 'Error: "%s" is not available' % (fname))
	print(end_page())
	exit()
if len(lines) < 1:
	print(html + 'Error: "%s" has no data available' % (fname))
	print(end_page())
	exit()

l1 = lines[:1]
lines = lines[1:]
vs = l1[0].split('\t')
headers = []
html += add_form(form['fn'].value,filter)
html += '<table cellpadding="4" cellspacing="0">'
html += '<tr style="background-color: black;color: white;">'
for v in vs:
	if v == 'Modifications' or v == 'Sequence':
		html += '<td style="text-align: left;"><b>%s</b></td>' % (v)
	else:
		html += '<td><b>%s</b></td>' % (v)
	headers.append(v)
html += '</tr>'
tfilter = filter
if filter.find('\\') == -1 and filter.find('^') == -1 and filter.find('$') == -1 and len(re.findall(r'\[.+\]',filter)) == 0 and filter.find('*') == -1 and filter.find('+') == -1:
	tfilter = re.escape(filter)
for i,l in enumerate(lines):
	if len(tfilter) == 0:
		pass
	else:
		l = re.sub(r'(%s)' % tfilter,r'<u>\1</u>',l,flags=re.IGNORECASE)

	vs = l.split('\t')
	if i % 2 == 1:
		html += '<tr style="background-color: #cccccc;">'
	else:
		html += '<tr>'
	for j,v in enumerate(vs):
		if headers[j] == 'Sequence' or headers[j] == 'Modifications':
			html += '<td style="text-align: left;">%s</td>' % (v)
		else:
			html += '<td>%s</td>' % (v)
	html += '</tr>'
html += '</table>'
print(html)
print(end_page())

