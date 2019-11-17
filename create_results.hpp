/*
#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
# Identifies kernels corresponding to spectra
#
*/

class id
{
public:
	id(void)	{}
	virtual ~id(void)	{}
	long sn;
	long peaks;
	vector<long> ks;
	double ri;
	long pm;
	long pz;
	long sc;
	long ions;
	id& operator=(const id &rhs)	{
		sn = rhs.sn;
		peaks = rhs.peaks;
		ri = rhs.ri;
		pm = rhs.pm;
		pz = rhs.pz;
		sc = rhs.sc;
		ions = rhs.ions;
		ks.clear();
		for(size_t k = 0; k < rhs.ks.size(); k++)	{
			ks.push_back(rhs.ks[k]);
		}
		return *this;
	}
};

class create_results
{
public:
	create_results(void);
	virtual ~create_results(void);
	bool create(map<string,string>& _p,
			vector<spectrum>& _s,
			kernels& _k,
			map<long,long>& _m);
	long size(void) { return (long)ids.size(); }
	vector<id> ids;
};



