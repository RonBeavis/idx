/*
#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
# Identifies kernels corresponding to spectra
#
*/

#include <algorithm>
#include <cmath>

typedef std::pair <unsigned int,unsigned int> sPair;

class spectrum
{
public:
	spectrum(void)	{pm = 0;pi=0;pz=0;sc=0;pt=0.0;}
	virtual ~spectrum(void)	{ clear(); }
	long pm;
	long pi;
	long pz;
	long sc;
	long isum;
	double pt;
	vector<pair<long,long>> mis;
	phmap::flat_hash_set<sPair> spairs;
	string desc;
	bool clear()	{sc = 0;pm = 0; pi = 0; pz = 0; sc=0; desc = ""; mis.clear();spairs.clear();return true;}
	spectrum& operator=(const spectrum &rhs)	{
		mis.clear();
		pm = rhs.pm;
		pi = rhs.pi;
		pz = rhs.pz;
		pt = rhs.pt;
		desc = rhs.desc;
		sc = rhs.sc;
		pair<long,long> p;
		for(size_t i = 0; i < rhs.mis.size(); i++)	{
			p.first = rhs.mis[i].first;
			p.second = rhs.mis[i].second;
			mis.push_back(p);
		}
		spairs.clear();
		spairs.insert(rhs.spairs.begin(),rhs.spairs.end());
		return *this;
	}
	bool condition(long _ires, long _l)	{
		double i_max = 0.0;
		double res = 1.0/(double)_ires;
		long m = 0;
		long i = 0;
		sort(mis.begin(), mis.end(), 
               		[](const auto& x, const auto& y) { return x.first < y.first; } );
		for(size_t a = 0; a < mis.size(); a++)	{
			m = mis[a].first;
			i = mis[a].second;		
			if(m < 150000 or (long)fabs(pm-m) < 45000 or (long)fabs(pm/pz- m) < 2000)	{
				continue;
			}
			if(i > i_max)	{
				i_max = i;
			}
		}
		vector<long> sMs;
		vector<long> sIs;
		i_max = i_max/100.0;
		for(size_t a = 0; a < mis.size(); a++)	{
			m = mis[a].first;
			if(m < 150000 or (long)fabs(pm-m) < 45000 or (long)fabs(pm/pz-m) < 2000)	{
				continue;
			}
			if((float)(mis[a].second)/i_max > 1.00)	{
				i = (long)(mis[a].second);
				if(!sMs.empty())	{
					if((long)fabs(sMs.back() - m) < 500)	{
						if(sIs.back() < i)	{
							sMs.pop_back();
							sIs.pop_back();
							sMs.push_back(m);
							sIs.push_back(i);
						}
					}
					else	{
						sMs.push_back(m);
						sIs.push_back(i);
					}
				}
				else	{
					sMs.push_back(m);
					sIs.push_back(i);
				}
			}
		}
		vector<pair<long,long> > pMs;
		pair<long,long> p(0,0);
		for(size_t a = 0; a < sMs.size();a++)	{
			p.first = sMs[a];
			p.second = sIs[a];
			pMs.push_back(p);
		}

		sort(pMs.begin(), pMs.end(), 
               		[](const auto& x, const auto& y) { return x.second > y.second; } );
		long max_l = 2 * (long)((0.5 + (float)pm)/100000.0);
		if(max_l > _l)	{
			max_l = _l;
		}
		while(pMs.size() > max_l)	{
			pMs.pop_back();
		}
		sort(pMs.begin(), pMs.end(), 
               		[](const auto& x, const auto& y) { return x.first < y.first; } );

//
//		generate a normalized set of spectrum masses
//
		long is = 0;
		mis.clear();
		sPair pr;
		double ptd = 1.0/70.0;
		pr.first = (unsigned int)(0.5+pm*ptd);
		for(size_t a=0; a < pMs.size();a++)	{
			m = pMs[a].first;
			p.first = (long)(0.5+(double)m*res);
			pr.second = p.first;
			spairs.insert(pr);
			p.second = pMs[a].second;
			mis.push_back(p);
			p.first -= 1;
			pr.second = p.first;
			spairs.insert(pr);
			mis.push_back(p);
			p.first += 2;
			pr.second = p.first;
			spairs.insert(pr);
			mis.push_back(p);
			is += 3*p.second;
		}
		isum = is;
		return true;
	}
};

class load_spectra
{
public:
	load_spectra(void);
	virtual ~load_spectra(void);
	bool load(map<string,string>& _p,vector<spectrum>& _s);
};



