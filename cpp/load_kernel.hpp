/*
#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
# Identifies kernels corresponding to spectra
#
*/

#include <algorithm>
typedef std::pair <int64_t,int64_t> kPair;

class kernels
{
public:
	kernels(void)	{}
	virtual ~kernels(void)	{}
	phmap::parallel_flat_hash_map<kPair,vector<int64_t> > kindex;
	phmap::parallel_flat_hash_set<int64_t> mvindex;
	void add_pair(kPair _v) {kindex[_v] = vector<int64_t>();}
	int64_t size(void)	{ return (int64_t)kindex.size();}
	void clear(void) { kindex.clear();}
};

class load_kernel
{
public:
	load_kernel(void);
	virtual ~load_kernel(void);
	bool load(map<string,string>& _p,load_spectra& _l,kernels& _k,map<int64_t,int64_t>& _m);
};



