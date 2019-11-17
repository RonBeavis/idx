/*
#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
# Identifies kernels corresponding to spectra
#
*/

class hypergeom	{
public:
	hypergeom(long _n,long _r,long _N)	{n = _n; r = _r; N = _N;}
	virtual ~hypergeom() {}
	long n;
	long N;
	long r;
	double pdf(long k)	{
		double top = st(n)+st(r)+st(N-n)+st(N-r);
		double bottom = st(N)+st(k)+st(n-k)+st(r-k)+st(N-n-r+k);
		double lp = top - bottom;
		return exp(lp);
	}
	double st(long _n)	{
		if(_n < 50)	{
			double f = 1.0;
			for(long i = 1; i <= _n; i++)    {
        			f *= (double)i;
   			 }
			return log(f);
		}
		double n = (double)_n;
		double pi = 3.1416;
		return (n*log(n) - n + log(sqrt(2*pi*n)));
	}
};

class mod
{
public:
	mod(void)	{}
	virtual ~mod(void)	{}
	long pos;
	string res;
	long mass;
	mod& operator=(const mod &rhs)	{
		pos = rhs.pos;
		res = rhs.res;
		mass = rhs.mass;
		return *this;
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
	long get_cells(double _pm,long _res);
	bool apply_model(long _r,string& _s,double _pm,long _ions,long _lspectrum,double& pscore,double& p);
	long low;
	long high;
	map<long,id> sv;
	map<long,set<long> > sdict;
	map<long,string> mt;
	map<long,vector<string> > odict;
	map<long,long> ppms;
	const long c13 = 1003;
	map<long,vector<double> > distribution;
	bool load_distribution(void)	{
		distribution.clear();
		vector<double> v = {1.3,1.0,0.6};
		distribution[800] = v;
		v = {3.2,2.0,0.8};
		distribution[900] = v;
		v = {4.9,2.7,1.0};
		distribution[1000] = v;
		v = {5.6,3.1,1.0};
		distribution[1100] = v;
		v = {6.2,3.3,1.1};
		distribution[1200] = v;
		v = {6.6,3.5,1.1};
		distribution[1300] = v;
		v = {7.1,3.8,1.2};
		distribution[1400] = v;
		v = {7.5,4.0,1.2};
		distribution[1500] = v;
		v = {7.8,4.1,1.3};
		distribution[1600] = v;
		v = {8.1,4.3,1.3};
		distribution[1700] = v;
		v = {8.2,4.4,1.3};
		distribution[1800] = v;
		v = {8.4,4.5,1.3};
		distribution[1900] = v;
		v = {8.5,4.6,1.3};
		distribution[2000] = v;
		v = {8.4,4.7,1.4};
		distribution[2100] = v;
		v = {8.6,4.7,1.4};
		distribution[2200] = v;
		v = {8.3,4.7,1.4};
		distribution[2300] = v;
		v = {8.2,4.6,1.4};
		distribution[2400] = v;
		v = {8.2,4.7,1.4};
		distribution[2500] = v;
		v = {8.1,4.8,1.4};
		distribution[2600] = v;
		v = {7.8,4.7,1.4};
		distribution[2700] = v;
		v = {7.5,4.7,1.4};
		distribution[2800] = v;
		v = {7.3,4.6,1.4};
		distribution[2900] = v;
		v = {7.3,4.7,1.5};
		distribution[3000] = v;
		v = {6.9,4.5,1.4};
		distribution[3100] = v;
		v = {6.5,4.4,1.4};
		distribution[3200] = v;
		v = {6.4,4.4,1.4};
		distribution[3300] = v;
		v = {5.8,4.1,1.4};
		distribution[3400] = v;
		v = {5.0,3.8,1.4};
		distribution[3500] = v;
		v = {5.5,4.1,1.4};
		distribution[3600] = v;
		v = {4.7,3.6,1.4};
		distribution[3700] = v;
		v = {4.2,3.3,1.4};
		distribution[3800] = v;
		v = {3.8,3.1,1.3};
		distribution[3900] = v;
		v = {3.7,3.0,1.3};
		distribution[4000] = v;
		return true;
	}
};



