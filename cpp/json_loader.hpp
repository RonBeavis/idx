/*
#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
# Identifies kernels corresponding to spectra
#
*/
/*
{"lv":0,"pm":1567745,"lb":"YGR014W","pre":"K","post":"I","beg":1246,"end":1258,"seq":"RMSVQESITQSMR","ns":[0,35,0,0],"bs":[156101,303137,390169,489237,617296,746338,833370,946454,1047502,1175561,1262593,1393633],"ys":[174112,305152,392184,520243,621290,734375,821407,950449,1078508,1177576,1264608,1411644],"mods":[["M",1247,15995]]}
*/
#include <string>
#include <vector>

class mod
{
public:
	mod(void)	{}
	virtual ~mod(void)	{}
	string res;
	long pos;
	long mass;
};

class json_loader
{
public:
	json_loader(void)	{}
	virtual ~json_loader(void)	{}
	long pm;
	long u;
	long h;
	long beg;
	long end;
	string lb;
	string pre;
	string post;
	string seq;
	vector<long> ns;
	vector<long> bs;
	vector<long> ys;
	vector<mod> mods;
	bool load(string& _s)	{
		mod m;
		string t = "\"h\":\"";
		size_t pos = _s.find(t);
		if(pos != _s.npos)	{
			pos += t.size()+1;
			h = atoi(_s.substr(pos,_s.find("\"",pos)).c_str());
		}
		t = "\"u\":\"";
		pos = _s.find(t);
		if(pos != _s.npos)	{
			pos += t.size()+1;
			u = atoi(_s.substr(pos,_s.find("\"",pos)).c_str());
		}
		t = "\"pm\":\"";
		pos = _s.find(t);
		if(pos != _s.npos)	{
			pos += t.size()+1;
			pm = atoi(_s.substr(pos,_s.find("\"",pos)).c_str());
		}
		t = "\"beg\":\"";
		pos = _s.find(t);
		if(pos != _s.npos)	{
			pos += t.size()+1;
			beg = atoi(_s.substr(pos,_s.find("\"",pos)).c_str());
		}
		t = "\"end\":\"";
		pos = _s.find(t);
		if(pos != _s.npos)	{
			pos += t.size()+1;
			end = atoi(_s.substr(pos,_s.find("\"",pos)).c_str());
		}
		t = "\"lb\":\"";
		pos = _s.find(t);
		if(pos != _s.npos)	{
			pos += t.size()+1;
			lb = _s.substr(pos,_s.find("\"",pos));
		}
		t = "\"pre\":\"";
		pos = _s.find(t);
		if(pos != _s.npos)	{
			pos += t.size()+1;
			pre = _s.substr(pos,_s.find("\"",pos));
		}
		t = "\"post\":\"";
		pos = _s.find(t);
		if(pos != _s.npos)	{
			pos += t.size()+1;
			post = _s.substr(pos,_s.find("\"",pos));
		}
		t = "\"seq\":\"";
		pos = _s.find(t);
		if(pos != _s.npos)	{
			pos += t.size()+1;
			seq = _s.substr(pos,_s.find("\"",pos));
		}
	}
};




