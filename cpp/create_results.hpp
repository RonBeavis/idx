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
	int64_t sn;
	int64_t peaks;
	vector<int64_t> ks;
	double ri;
	int64_t pm;
	int64_t pz;
	int64_t sc;
	int64_t ions;
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
			load_spectra& _l,
			kernels& _k,
			map<int64_t,int64_t>& _m);
	int64_t size(void) { return (int64_t)ids.size(); }
	vector<id> ids;
};



