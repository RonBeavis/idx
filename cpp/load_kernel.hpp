/*
#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
# Identifies kernels corresponding to spectra
#
*/

#include <algorithm>

class kernel
{
public:
	kernel(void)	{}
	virtual ~kernel(void)	{}
	bool clear()	{}
//	kernel& operator=(const kernel &rhs)	{
//	}
};

class load_kernel
{
public:
	load_kernel(void);
	virtual ~load_kernel(void);
	bool load(map<string,string>& _p,vector<spectrum>& _s,vector<kernel>& _k,map<long,long>& _m);
};



