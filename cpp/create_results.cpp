/*
#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
# Loads a spectrum file into a vector of spectrum objects
#
*/
#include "pch.h"

#include <fstream>
#include <cstdio>
#include <iostream>
#include <sys/stat.h>
#include <map>
#include <string>
#include <unordered_map>
#include <set>
#include <vector>
#include "parallel_hashmap/phmap.h"
using namespace std;
#include "load_spectra.hpp"
#include "load_kernel.hpp"
#include "create_results.hpp"

create_results::create_results(void) {
}

create_results::~create_results(void) {
}

bool create_results::create(map<string, string>& _p,
	load_spectra& _l,
	kernels& _k,
	map<long, long>& _m) {
	long z = 1;
	double pt = 1.0 / 70.0;
	double ppm = atof(_p["parent tolerance"].c_str()) / 1.0e06;
	vector<long> dvals{ -1,0,1 };
	long d = 0;
	long m = 0;
	long m2 = 0;
	double c13 = 1003.0;
	double rm = 0.0;
	long pm = 0;
	long pz = 0;
	long idi = 0;
	long cpm = 0;
	vector<long> idx;
	vector<vector<long> > ident;
	vector<long> ims;
	bool use2 = false;
	bool use3 = false;
	auto itkend = _k.kindex.end();
	auto itk = _k.kindex.end();
	kPair pv;
	vector<double> mvs;
	vector<double> ms;
	size_t pos = 0;
	for (size_t s = 0; s < _l.spectra.size(); s++) {
		if (s != 0 and s % 500 == 0) {
			cout << '.';
			cout.flush();
		}
		if (s != 0 and s % 10000 == 0) {
			cout << " " << s << endl;
			cout.flush();
		}
		ident.clear();
		ims.clear();
		rm = (double)_l.spectra[s].pm;
		pm = (long)(0.5 + rm * pt);
		pz = (long)_l.spectra[s].pz;
		idx.clear();
		idi = 0;
		use2 = (pm > 1500000 and pz == 2);
		use3 = (pz > 2);
		ms.clear();
		mvs.clear();
		mvs.push_back(rm);
		if (rm > 1500000.0) {
			mvs.push_back(rm - c13);
		}
		for (size_t a = 0; a < mvs.size(); a++) {
			rm = mvs[a];
			cpm = (long)(0.5 + rm * pt);
			for (size_t n = 0; n < dvals.size(); n++) {
				d = dvals[n];
				pv.first = (long)(cpm + d);
				if(_k.mvindex.find(pv.first) == _k.mvindex.end())	{
					continue;
				}
				for(size_t b=0; b < _l.spectra[s].mis.size(); b++)	{
					m = _l.spectra[s].mis[b].first;
					m2 = m * 2;
					pv.second = m;
					itk = _k.kindex.find(pv);
					idx.clear();
					if(itk !=  _k.kindex.end())	{
						idx.insert(idx.end(),_k.kindex[pv].begin(),_k.kindex[pv].end());
						idi = _l.spectra[s].mis[b].second;
					}
					else if(use2 and m2 > 100000)	{
						pv.second = m2;
						itk = _k.kindex.find(pv);
						if(itk !=  _k.kindex.end())	{
							idx.insert(idx.end(),_k.kindex[pv].begin(),_k.kindex[pv].end());
							idi = _l.spectra[s].mis[b].second;
						}
					}
					else if(use3)	{
						pv.second = m2;
						itk = _k.kindex.find(pv);
						if(itk !=  _k.kindex.end())	{
							idx.insert(idx.end(),_k.kindex[pv].begin(),_k.kindex[pv].end());
							idi = _l.spectra[s].mis[b].second;
						}
					}
					else	{
						continue;
					}
					ims.push_back(idi);
					ident.push_back(idx);
					ms.push_back(rm);
				}
			}
		}
		if(ident.empty())	{
			continue;
		}
		long isum = 0;
		for(size_t a = 0; a < _l.spectra[s].mis.size(); a++)	{
			isum += _l.spectra[s].mis[a].second;
		}
		double total = (double)isum/3;
		set<long> aok;
		for(size_t b = 0; b < ident.size(); b++)	{
			for(size_t c = 0; c < ident[b].size(); c++)	{
				long a = ident[b][c];
				if(fabs(ms[b]-_m[a]) < ppm*ms[b])	{
					aok.insert(a);
				}
			}
		}
		map<long,long> ans;
		map<long,long> aint;
		for(size_t b = 0; b < ident.size(); b++)	{
			set<long> sv(ident[b].begin(),ident[b].end());
			set<long>::iterator it = sv.begin();
			while(it != sv.end())	{
				if(aok.find(*it) == aok.end())	{
					ans[*it] = 0;
					aint[*it] = 0;
				}
				else if(ans.find(*it) != ans.end()) {
					ans[*it] += 1;
					aint[*it] += ims[b];
				}
				else	{
					ans[*it] = 1;
					aint[*it] = ims[b];
				}
				it++;
			}
		}
		long mn = 0;
		vector<long> mv;
		vector<long> iv;
		map<long,long>::iterator itans = ans.begin();
		while(itans != ans.end())	{
			if(itans->second > mn)	{
				mn = itans->second;
				mv.clear();
				mv.push_back(itans->first);
				iv.clear();
				iv.push_back(aint[itans->first]);
			}
			else if(itans->second == mn)	{
				mv.push_back(itans->first);
				iv.push_back(aint[itans->first]);
			}
			itans++;
		}
		if(mn > 4)	{
			id r;
			long max_i = *max_element(iv.begin(),iv.end());
			for(size_t b = 0; b < mv.size(); b++)	{
				if(iv[b] < max_i)	{
					continue;
				}
				r.ks.push_back(mv[b]);
			}
			r.sn = z;
			r.peaks = mn;
			r.ri = max_i/total;
			r.pm = _l.spectra[s].pm;
			r.pz = _l.spectra[s].pz;
			r.sc = _l.spectra[s].sc;
			r.ions = (long)_l.spectra[s].mis.size()/3;
			ids.push_back(r);
		}
		z += 1;
	}
	cout << endl;
	cout.flush();
	return true;
}

