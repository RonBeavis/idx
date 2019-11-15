/*
#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
# Identifies kernels corresponding to spectra
#
*/

class mod
{
public:
	mod(void)	{}
	virtual ~mod(void)	{}
	long pos;
	string res;
	double mass;
	mod& operator=(const mod &rhs)	{
		pos = rhs.pos;
		res = rhs.res;
		mass = rhs.mass;
	}
	bool operator<( const mod& rhs ) const { 
		return pos < rhs.pos; 
	}
};

class create_output
{
public:
	create_output(void);
	virtual ~create_output(void);

	bool create(map<string,string>& _p,create_results& _cr);
private:
	bool load_mods(void);
	bool find_window(void);
	bool apply_model(long _r,string& _s,double _pm,long _ions,long _lspectrum,double& pscore,double& p);
	long low;
	long high;
	map<long,id> sv;
	map<long,set<long> > sdict;
	map<long,string> mt;
	map<long,vector<string> > odict;
	map<long,long> ppms;
	const long c13 = 1003;
};



