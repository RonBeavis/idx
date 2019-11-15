/*
#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
# Identifies kernels corresponding to spectra
#
*/
// #include "pch.h"
#include <iostream>
#include <cstdio>
#include <sys/stat.h>
#include <string>
#include <cstring>
#include <map>
#include <unordered_map>
#include <set>
#include <vector>
#include <ctime>
#include <chrono>
using namespace std;
using namespace std::chrono;
#include "load_spectra.hpp"
#include "load_kernel.hpp"
#include "create_results.hpp"
#include "create_output.hpp"

// for convenience
int main(int argc, char* argv[])	{
	if(argc < 4)	{
		cout << "usage:\t>idX SPECTRA_FILE KERNEL_FILE OUTPUT_FILE (high|medium|low*)\n";
		return 0;
	}
	map<string,string> params;
	string version = "idX, 2020.1";
	params["version"] = version;

	long fragment_tolerance = 400;
	if(argc > 4 and strcmp(argv[4],"high") == 0)	{
		fragment_tolerance = 20;
	}
	else if (argc > 4 and strcmp(argv[4],"medium") == 0)	{
		fragment_tolerance = 50;
	}
	char ptemp[64];
	sprintf(ptemp,"%li",fragment_tolerance);
	params["fragment tolerance"] = ptemp;

	string spectrum_file = argv[1];
	struct stat fbuffer;
    	if(stat(spectrum_file.c_str(), &fbuffer) != 0)	{
		cout << "Error (idx:0001): spectrum file \"" << spectrum_file << "\" does not exist\n";
		return 0;
	}
	params["spectrum file"] = spectrum_file;

	string kernel_file = argv[2];
    	if(stat(kernel_file.c_str(), &fbuffer) != 0)	{
		cout << "Error (idx:0002): kernel file \"" << kernel_file << "\" does not exist\n";
		return 0;
	}
	params["kernel file"] = kernel_file;

	string output_file = argv[3];
	params["output file"] = output_file;

	long maximum_spectra = -1;
	if(argc > 5)	{
		maximum_spectra = atoi(argv[5]);
	}
	sprintf(ptemp,"%li",maximum_spectra);
	params["maximum spectra"] = ptemp;

	long parent_tolerance = 20;
	sprintf(ptemp,"%li",parent_tolerance);
	params["parent tolerance"] = ptemp;

	cout << "\nstart ...\nidX parameters" << "\n";
	if(maximum_spectra != -1)	{
		cout << "\t   max spectra: " << maximum_spectra << "\n";
	}
	else	{
		cout << "\t   max spectra: unlimited" << "\n";
	}
	cout << "\t  fragment tol: " << fragment_tolerance << " mDa\n";
	cout << "\t    parent tol: " << params["parent_tolerance"] << " ppm\n";
	cout << "\t spectrum file: " << spectrum_file << "\n";
	cout << "\t   kernel file: " << kernel_file << "\n";
	cout << "\t   output file: " << output_file << "\n";
	cout << "\t       version: " << version << "\n";
	cout << "load & index spectra"  << "\n";
	cout.flush();
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	vector<spectrum> spectra;
	load_spectra ls;
	if(!ls.load(params,spectra))	{
		cout << "Error (idx:0003): failed to load spectrum file \"" << spectrum_file << "\"\n";
		return 0;
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	cout << "	   spectra = " << spectra.size() << "\n";
	cout << "	spectra &Delta;T = " << duration_cast<milliseconds>(t2 - t1).count()/1000.0 << " s\n";
	sprintf(ptemp,"%li",spectra.size());
	params["spectra"] = ptemp;
	t1 = high_resolution_clock::now();
	cout << "load & index kernel"  << "\n";
	cout.flush();
	kernels kindex;
	map<long,long> mindex;
	load_kernel lk;
	if(!lk.load(params,spectra,kindex,mindex))	{
		cout << "Error (idx:0005): failed to load kernel file \"" << spectrum_file << "\"\n";
		return 0;
	}
	t2 = high_resolution_clock::now();
	cout << "	   kernels = " << kindex.size() << "\n";
	cout << "	kernels &Delta;T = " << duration_cast<milliseconds>(t2 - t1).count()/1000.0 << " s\n";
	t1 = high_resolution_clock::now();
	cout << "perform ids"  << "\n";
	cout.flush();
	create_results cr;
	if(!cr.create(params,spectra,kindex,mindex))	{
		cout << "Error (idx:0006): failed to create results " << "\n";
		return 0;
	}
	t2 = high_resolution_clock::now();
	cout << "	   results = " << cr.size() << "\n";
	cout << "	results &Delta;T = " << duration_cast<milliseconds>(t2 - t1).count()/1000.0 << " s\n";
	cout << "create report"  << "\n";
	cout.flush();
	create_output co;
	if(!co.create(params,cr))	{
		cout << "Error (idx:0007): failed to create output " << "\n";
		return 0;
	}
	t2 = high_resolution_clock::now();
	cout << "... done\n";
	return 0;
}

