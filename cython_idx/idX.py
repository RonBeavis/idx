#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
# Identifies kernels corresponding to spectra
#
# idX version 2019.08.10.02
#

import ujson
import time
import gzip
import sys
import datetime

# import the method that deals with spectrum file formats
from spectraX import load_spectra
# import the method for the output of results to a file
from reportX import report_ids
from kernelX import index_kernel
from createX import create_ids

version = '2019.09.01'
#
# Coordinate the identification process, print job stats and progress 
#

def main():
	if len(sys.argv) < 4:
		print('usage:\n\t>python3 idX.py SPECTRA_FILE KERNEL_FILE OUTPUT_FILE (high|medium|low*)')
		exit()
	start = time.time()
	# record relavent parameters
	param = {}
	#fragment tolerance in millidaltons
	param['fragment mass tolerance'] = float(400)
	try:
		if sys.argv[4] == 'high':
			param['fragment mass tolerance'] = float(20)
		elif sys.argv[4] == 'low':
			param['fragment mass tolerance'] = float(400)
		elif sys.argv[4] == 'medium':
			param['fragment mass tolerance'] = float(100)
		else:
			print('ERROR: argument 4 must be high or low, not "%s"'% (sys.argv[4]))
			exit()
	except:
		pass
	param['maximum spectra'] = -1
	try:
		param['maximum spectra'] = int(sys.argv[5])
	except:
		pass
	# parent tolerance in ppm
	param['parent mass tolerance'] = float(20)
	spectra = []
	# report files named on command line
	print('\nstart ...\nidX parameters')
	if param['maximum spectra'] != -1:
		print('\t   max spectra: %i mDa' % (param['maximum spectra']))
	else:
		print('\t   max spectra: unlimited')
	
	print('\t  fragment tol: %i mDa' % (param['fragment mass tolerance']))
	print('\t    parent tol: %i ppm' % (param['parent mass tolerance']))
	print('\t spectrum file: %s' % (sys.argv[1]))
	print('\t   kernel file: %s' % (sys.argv[2]))
	print('\t   output file: %s' % (sys.argv[3]))
	print('\t       version: %s' % (version))
	print('\t      run time: %s' % (str(datetime.datetime.now())))
	param['spectrum file'] = sys.argv[1]
	print('load & index spectra')
	# read the spectrum file and perform all necessary spectrum conditioning
	spectra = load_spectra(param['spectrum file'],param)
	if param['maximum spectra'] != -1:
		spectra = spectra[0:param['maximum spectra']]
	if len(spectra) == 0:
		print('exiting: 0 spectra found')
		print('done')
		exit()
	param['spectra'] = len(spectra)
	delta = time.time()-start
	start = time.time()
	print('\n\t   spectra = %i' % (len(spectra)))
	print('\tspectra &Delta;T = %.1f s' % (delta))
	param['kernel file'] = sys.argv[2]
	print('load & index kernel')
	# read the kernel file and create an index of peptide fragmentation patterns
	(ki,mi) = index_kernel(param,spectra)
	delta = time.time()-start
	start = time.time()
	print('\n\t   kernels = %i' % (len(ki)))
	print('\t &Delta;T = %.1f s' % (delta))
	print('perform ids')
	# generate a list of identifications for the spectra using the kernel index
	ids = create_ids(ki,mi,spectra,param)
	# free memory associated with indexes and spectra
	delta = time.time()-start
	start = time.time()
	print('\tid &Delta;T = %.3f s' % (delta))
	if len(spectra) > 0:
		print('\t   &delta;T = %.0f microseconds' % (1.0e06*delta/len(spectra)))
	else:
		pass
	# simple reporting of the kernels assigned to spectra
	print('release memory')
	ki = None
	spectra = None
	print('\tdone')
	param['output file'] = sys.argv[3]
	print('create report')
	report_ids(ids,param)
	print('... done')

if __name__== "__main__":
	main()


