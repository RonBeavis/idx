/*
#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
# Identifies kernels corresponding to spectra
#
*/

#include <algorithm>
typedef std::pair <unsigned int,unsigned int> kPair;

class kernels
{
public:
	kernels(void)	{}
	virtual ~kernels(void)	{}
	phmap::flat_hash_map<kPair,vector<unsigned int> > kindex;
	phmap::flat_hash_set<unsigned int> mvindex;
	void add_pair(kPair _v) {kindex[_v] = vector<unsigned int>();}
	long size(void)	{ return (long)kindex.size();}
	void clear(void) { kindex.clear();}
};

class load_kernel
{
public:
	load_kernel(void);
	virtual ~load_kernel(void);
	bool load(map<string,string>& _p,vector<spectrum>& _s,kernels& _k,map<long,long>& _m);
};



