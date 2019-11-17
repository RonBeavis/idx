/*
#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
# Identifies kernels corresponding to spectra
#
*/

#include <algorithm>

class kernels
{
public:
	kernels(void)	{}
	virtual ~kernels(void)	{}
	unordered_map<unsigned int,unordered_map<unsigned int,vector<unsigned int> > > kindex;
	long size(void)	{ return (long)kindex.size();}
	bool add_map(long _m)	{ 
							kindex.insert(pair<unsigned int,
											unordered_map<unsigned int,vector<unsigned int> > >(_m,
											unordered_map<unsigned int, vector<unsigned int> >())); 
							return true;
	}
	bool add_pair(long _mv, long _val) { kindex[_mv].insert(pair<unsigned int, vector<unsigned int> >(_val, vector<unsigned int>())); return true; }
	bool clear(void) { kindex.clear();}
//	kernel& operator=(const kernel &rhs)	{
//		
//	}
};

class load_kernel
{
public:
	load_kernel(void);
	virtual ~load_kernel(void);
	bool load(map<string,string>& _p,vector<spectrum>& _s,kernels& _k,map<long,long>& _m);
};



