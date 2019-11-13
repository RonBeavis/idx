/*
#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
# Loads a spectrum file into a vector of spectrum objects
#
*/

#include "rapidjson/document.h"
#include <fstream>
#include <cstdio>
#include <iostream>
#include <sys/stat.h>
#include <unordered_map>
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

bool load_kernel::load(map<string,string>& _params,vector<spectrum>& _spectra,kernels& _kernels,map<long,long>& _mindex)	{
	ifstream istr;
	istr.open(_params["kernel file"]);
	if(istr.fail())	{
		return false;
	}
	double ft = 1.0/atof(_params["fragment tolerance"].c_str());
	double pt = 1.0/70.0;
	set<long> sp_set;
	for(size_t a = 0; a < _spectra.size();a++)	{
		sp_set.insert((long)(0.5 + _spectra[a].pm*pt));
	}
	auto itsp = sp_set.end();
	string line;
	using namespace rapidjson;
	long s = 0;
	long skipped = 0;
	long hmatched = 0;
	long pm = 0;
	long mv = 0;
	long cv = 0;
	long u = 0;
	const long c13 = 1003;
	unsigned int val = 0;
	while(getline(istr,line))	{
//		auto js = json::parse(line);
//		auto ij = js.find("pm");
		Document js;
    		js.Parse(line.c_str(),line.length());
		if(!js.HasMember("pm"))	{
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
		pm = (long)js["pm"].GetInt();
		mv = (long)(0.5+pm*pt);
		cv = (long)(0.5+(pm+c13)*pt);
		if(sp_set.find(mv) == itsp or sp_set.find(mv-1) == itsp or sp_set.find(mv+1) == itsp)	{
		}
		else if(pm > 1500000 and (sp_set.find(cv) == itsp or sp_set.find(cv-1) == itsp or sp_set.find(cv+1) == itsp))	{
		}
		else	{
			skipped += 1;
			continue;
		}
		u = (unsigned int)js["u"].GetInt();
		_mindex[u] = pm;
		if(_kernels.kindex.find(mv) == _kernels.kindex.end())	{
			_kernels.kindex.insert(pair<unsigned int,unordered_map<unsigned int,vector<unsigned int> > >(mv,unordered_map<unsigned int,vector<unsigned int> >()));
		}
		const Value& jbs = js["bs"];
		for(SizeType a = 0; a < jbs.Size();a++)	{
			val = (unsigned int)(0.5+jbs[a].GetDouble()*ft);
			if(_kernels.kindex[mv].find(val) == _kernels.kindex[mv].end())	{
				_kernels.kindex[mv].insert(pair<unsigned int,vector<unsigned int> >(val,vector<unsigned int>()));
			}
			_kernels.kindex[mv][val].push_back(u);
		}
		const Value& jys = js["ys"];
		for(size_t a = 0; a < jys.Size();a++)	{
			val = (unsigned int)(0.5+jys[a].GetDouble()*ft);
			if(_kernels.kindex[mv].find(val) == _kernels.kindex[mv].end())	{
				_kernels.kindex[mv].insert(pair<unsigned int,vector<unsigned int> >(val,vector<unsigned int>()));
			}
			_kernels.kindex[mv][val].push_back(u);
		}
	}
	cout << "\n";
	cout << "skipped: " << skipped << ", hmatches: " << hmatched << endl;
	cout.flush();
	istr.close();
	return true;
}

