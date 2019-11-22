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

typedef std::pair <int64_t,int64_t> sPair;

class spectrum
{
public:
	spectrum(void)	{pm = 0;pi=0;pz=0;sc=0;pt=0.0;}
	virtual ~spectrum(void)	{ clear(); }
	int64_t pm;
	int64_t pi;
	int64_t pz;
	int64_t sc;
	int64_t isum;
	double pt;
	vector<pair<int64_t,int64_t>> mis;
	phmap::parallel_flat_hash_set<sPair> spairs;
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
		pair<int64_t,int64_t> p;
		for(size_t i = 0; i < rhs.mis.size(); i++)	{
			p.first = rhs.mis[i].first;
			p.second = rhs.mis[i].second;
			mis.push_back(p);
		}
		spairs.clear();
		spairs.insert(rhs.spairs.begin(),rhs.spairs.end());
		return *this;
	}
	bool condition(int64_t _ires, int64_t _l)	{
		double i_max = 0.0;
		const double res = 1.0/(double)_ires;
		int64_t m = 0;
		int64_t i = 0;
		sort(mis.begin(), mis.end(), 
               		[](const auto& x, const auto& y) { return x.first < y.first; } );
		for(size_t a = 0; a < mis.size(); a++)	{
			m = mis[a].first;
			i = mis[a].second;		
			if(m < 150000 or (int64_t)fabs(pm-m) < 45000 or (int64_t)fabs(pm/pz- m) < 2000)	{
				continue;
			}
			if(i > i_max)	{
				i_max = (double)i;
			}
		}
		vector<int64_t> sMs;
		vector<int64_t> sIs;
		i_max = i_max/100.0;
		for(size_t a = 0; a < mis.size(); a++)	{
			m = mis[a].first;
			if(m < 150000 or fabs(pm-m) < 45000.0 or fabs(pm/pz-m) < 2000.0)	{
				continue;
			}
			if((double)mis[a].second/i_max >= 1.00)	{
				i = mis[a].second;
				if(!sMs.empty())	{
					if(fabs(sMs.back() - m) < 500.0)	{
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
		vector<pair<int64_t,int64_t> > pMs;
		pair<int64_t,int64_t> p(0,0);
		for(size_t a = 0; a < sMs.size();a++)	{
			p.first = sMs[a];
			p.second = sIs[a];
			pMs.push_back(p);
		}
		sort(pMs.begin(), pMs.end(), 
               		[](const auto& x, const auto& y) { if(x.second > y.second) return true; if(x.second == y.second) return x.first < y.first; return false;} );
		int64_t max_l = 2 * (int64_t)((0.5 + (double)pm)/100000.0);
		if(max_l > _l)	{
			max_l = _l;
		}
		while(pMs.size() > (size_t)max_l)	{
			pMs.pop_back();
		}
		sort(pMs.begin(), pMs.end(), 
               		[](const auto& x, const auto& y) { return x.first < y.first; } );

//
//		generate a normalized set of spectrum masses
//
		int64_t is = 0;
		mis.clear();
		sPair pr;
		const double ptd = 1.0/70.0;
		pr.first = (int64_t)(0.5+(double)pm*ptd);
		for(size_t a=0; a < pMs.size();a++)	{
			m = pMs[a].first;
			p.first = (int64_t)(0.5+(double)m*res);
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
	bool load(map<string,string>& _p);
	vector<spectrum> spectra;
	void set_max(const int64_t _max)	{
		if(spectra.size() <= (size_t)_max)	{
			return;
		}
		spectra.erase(spectra.begin()+_max,spectra.end());
		return;
	}

};



