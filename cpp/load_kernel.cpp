/*
#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
# Loads a spectrum file into a vector of spectrum objects
#
*/

#include "nlohmann/json.hpp"
#include <fstream>
#include <cstdio>
#include <iostream>
#include <sys/stat.h>
#include <map>
#include <string>
#include <set>
#include <vector>
using namespace std;
#include "load_spectra.hpp"
#include "load_kernel.hpp"

load_kernel::load_kernel(void)	{
}

load_kernel::~load_kernel(void)	{
}

bool load_kernel::load(map<string,string>& _params,vector<spectrum>& _spectra,vector<kernel>& _kernels,map<long,long>& _mindex)	{
	ifstream istr;
	istr.open(_params["kernel file"]);
	if(istr.fail())	{
		return false;
	}
	const long ft = (long)atoi(_params["fragment tolerance"].c_str());
	const long pt = 70;
	set<long> sp_set;
	for(size_t a = 0; a < _spectra.size();a++)	{
		sp_set.insert((long)(0.5 + _spectra[a].pm/pt));
	}
	auto itsp = sp_set.end();
	string line;
	using json = nlohmann::json;
	long s = 0;
	kernel k;
	long skipped = 0;
	long hmatched = 0;
	long pm = 0;
	long mv = 0;
	long u = 0;
	while(getline(istr,line))	{
		auto js = json::parse(line);
		auto ij = js.find("pm");
		if(ij == js.end())	{
			continue;
		}
		s++;
		if(s % 10000 == 0)	{
			cout << ".";
			cout.flush();
		}
		if(js["u"] != js["h"])	{
			hmatched++;
			continue;
		}
		pm = (long)js["pm"];
		mv = (long)(0.5+pm/pt);
		
		if(sp_set.find(mv) == itsp && sp_set.find(mv-1) == itsp && sp_set.find(mv+1) == itsp)	{
			skipped += 1;
			continue;
		}
		u = (long)js["u"];
		_mindex[u] = pm;
	}
	cout << "\n";
	cout << "skipped: " << skipped << ", hmatches: " << hmatched << endl;
	cout.flush();
	istr.close();
	return true;
}

