/*
#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
# Loads a spectrum file into a vector of spectrum objects
#
*/
#include "pch.h"

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
#include "parallel_hashmap/phmap.h"
using namespace std;
#include "load_spectra.hpp"
#include "load_kernel.hpp"
#include "json_loader.hpp"

load_kernel::load_kernel(void)	{
}

load_kernel::~load_kernel(void)	{
}

bool load_kernel::load(map<string,string>& _params,load_spectra& _l,kernels& _kernels,map<long,long>& _mindex)	{
	ifstream istr;
	istr.open(_params["kernel file"]);
	if(istr.fail())	{
		return false;
	}
	string line;
	using namespace rapidjson;
	const double ft = 1.0/atof(_params["fragment tolerance"].c_str());
	const double pt = 1.0/70.0;
	const double ppm = 2.0E-5;
	set<long> sp_set;
	phmap::parallel_flat_hash_set<sPair> spairs;
	for(size_t a = 0; a < _l.spectra.size();a++)	{
		sp_set.insert(_l.spectra[a].pm);
		spairs.insert(_l.spectra[a].spairs.begin(),_l.spectra[a].spairs.end());
	}
	auto itsp = sp_set.end();
	auto itppm = sp_set.end();
	long skipped = 0;
	long hmatched = 0;
	long pm = 0;
	long mv = 0;
	long cv = 0;
	long u = 0;
	const long c13 = 1003;
	long val = 0;
	long lines = 0;
	bool skip = true;
	long lower = 0;
	long delta = 0;
	kPair pr;
	auto itsend = spairs.end();

/*	static const size_t val_buf_sz{ 32768 };
	char val_buf[val_buf_sz];
	rapidjson::MemoryPoolAllocator<> val_alloc{ &val_buf[0], sizeof(val_buf) };

	static const size_t parse_buf_sz{ 32768 };
	char p_buf[parse_buf_sz];
	rapidjson::MemoryPoolAllocator<> parse_alloc{ &p_buf[0], sizeof(p_buf) };

	rapidjson::GenericDocument<
		rapidjson::UTF8<>,
		rapidjson::MemoryPoolAllocator<>,
		rapidjson::MemoryPoolAllocator<> > js{ &val_alloc, sizeof(p_buf), &parse_alloc };
	json_loader jload;	*/
	while(getline(istr,line))	{
		if(lines != 0 and lines % 10000 == 0)	{
			cout << '.';
			cout.flush();
		}
		if(lines != 0 and lines % 200000 == 0)	{
			cout << " " << lines << endl;
			cout.flush();
		}
		lines++;
		Document js;
   		js.Parse(line.c_str(),line.length());
//		js.SetNull();
//		val_alloc.Clear();
//		parse_alloc.Clear();
//		jload.load(line);
//		continue;
		if(!js.HasMember("pm"))	{
			continue;
		}
		if(js["u"].GetInt() != js["h"].GetInt())	{
			hmatched++;
			continue;
		}
		pm = (long)js["pm"].GetInt();
		mv = (long)(0.5+(double)pm*pt);
		cv = (long)(0.5+(double)(pm+c13)*pt);
		delta = (long)(0.5+(double)pm*ppm);
		lower = pm-delta;
		itppm = sp_set.lower_bound(lower);
		skip = true;
		if(itppm != itsp and (*itppm-lower) <= 2*delta)	{
			skip = false;
		}
		lower = pm+c13-delta;
		itppm = sp_set.lower_bound(lower);
		if(itppm != itsp and (*itppm-lower) <= 2*delta)	{
			skip = false;
		}
		if(skip)	{
			skipped++;
			continue;
		}
		u = (long)js["u"].GetInt();
		_mindex[u] = pm;
		const Value& jbs = js["bs"];
		size_t vpos = 0;
		pr.first = (long)mv;
		for(SizeType a = 0; a < jbs.Size();a++)	{
			val = (long)(0.5+jbs[a].GetDouble()*ft);
			pr.second = val;
			if(spairs.find(pr) == spairs.end())	{
				continue;
			}
			if(_kernels.kindex.find(pr) == _kernels.kindex.end())	{
				_kernels.add_pair(pr);
			}
			_kernels.mvindex.insert((long)mv);
			_kernels.kindex[pr].push_back(u);
		}
		const Value& jys = js["ys"];
		for(SizeType a = 0; a < jys.Size();a++)	{
			val = (long)(0.5+jys[a].GetDouble()*ft);
			pr.second = val;
			if(spairs.find(pr) == spairs.end())	{
				continue;
			}
			if(_kernels.kindex.find(pr) == _kernels.kindex.end())	{
				_kernels.add_pair(pr);
			}
			_kernels.mvindex.insert((long)mv);
			_kernels.kindex[pr].push_back(u);
		}
//		js.SetNull();
//		val_alloc.Clear();
//		parse_alloc.Clear();
	}
	cout << "\n";
//	cout << "skipped: " << skipped << ", hmatches: " << hmatched << ", found: " << found << ", lines: " << lines << endl;
	cout.flush();
	istr.close();
	return true;
}

