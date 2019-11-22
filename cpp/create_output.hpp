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
	hypergeom(int64_t _n,int64_t _r,int64_t _N)	{n = _n; r = _r; N = _N;}
	virtual ~hypergeom() {}
	int64_t n;
	int64_t N;
	int64_t r;
	double pdf(int64_t k)	{
		double top = st(n)+st(r)+st(N-n)+st(N-r);
		double bottom = st(N)+st(k)+st(n-k)+st(r-k)+st(N-n-r+k);
		double lp = top - bottom;
		return exp(lp);
	}
	double st(int64_t _n)	{
		if(_n < 50)	{
			double f = 1.0;
			for(int64_t i = 1; i <= _n; i++)    {
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
	int64_t pos;
	string res;
	int64_t mass;
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
	int64_t get_cells(double _pm,int64_t _res);
	bool apply_model(int64_t _r,string& _s,double _pm,int64_t _ions,int64_t _lspectrum,double& pscore,double& p);
	double low;
	double high;
	map<int64_t,id> sv;
	map<int64_t,set<int64_t> > sdict;
	map<int64_t,string> mt;
	map<int64_t,vector<string> > odict;
	map<int64_t,int64_t> ppms;
	const int64_t c13 = 1003;
	map<int64_t,vector<double> > distribution;
	double get_ppm(string& t)	{
		size_t s = t.find("\t");
		s = t.find("\t",s+1);
		s = t.find("\t",s+1);
		return atof((t.substr(s+1,t.size()-1)).c_str());
	}
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



