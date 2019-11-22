/*
#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
# Loads a spectrum file into a vector of spectrum objects
#
*/
#include "pch.h"

#include <fstream>
#include <cstdio>
#include <iostream>
#include <sys/stat.h>
#include <map>
#include <string>
#include <set>
#include <vector>
#include "parallel_hashmap/phmap.h"
using namespace std;
#include "load_spectra.hpp"

load_spectra::load_spectra(void)	{
}

load_spectra::~load_spectra(void)	{
}

bool load_spectra::load(map<string,string>& _params)	{
	ifstream istr;
	istr.open(_params["spectrum file"]);
	int64_t res = (int64_t)atoi(_params["fragment tolerance"].c_str());
	if(istr.fail())	{
		return false;
	}
	double pt = atoi(_params["parent tolerance"].c_str());
	size_t len = 1024;
	size_t size = 0;
	char *line = new char[len];
	string temp = "";
	string desc = "";
	char *pos = NULL;
	size_t equals = 0;
	double parent = 0.0;
	double charge = 1.0;
	string run_time = "";
	vector<double> masses;
	vector<double> intensities;
	spectrum sp;
	const double proton = 1.007276;
	int64_t scan = 0;
	int64_t s = 1;
	while(istr.good() && !istr.eof())	{
		istr.getline(line,len-1,'\n');
		temp = line;
		size = temp.size();
		equals = temp.find("=");
		if(temp.find("PEPMASS=") != temp.npos)	{
			temp = temp.substr(equals+1,size-equals+1);
			parent = atof(temp.c_str());
		}
		else if(temp.find("BEGIN IONS") != temp.npos)	{
			desc = "";
		}
		else if(temp.find("#") == 0)	{
			desc += line+1;
			desc += " ";
		}
		else if(temp.find("TITLE=") != temp.npos)	{
			desc += temp.substr(equals+1,size-equals+1);
			desc += " ";
		}
		else if(temp.find("RTINSECONDS=") != temp.npos)	{
			run_time = temp.substr(equals+1,size-equals+1);
			desc += "RTINSECONDS=";
			desc += run_time;
			desc += " ";
		}
		else if(temp.find("CHARGE=") != temp.npos)	{
			temp = temp.substr(equals+1,size-equals+1);
			charge = atof(temp.c_str());
		}
		else if(temp.find("SCANS=") != temp.npos)	{
			temp = temp.substr(equals+1,size-equals+1);
			scan = atoi(temp.c_str());
		}
		else if(atof(line) > 0.0)	{
			masses.push_back(atof(line));
			pos = line;
			while(*pos != '\0' && isspace(*pos) != 0)	{
				pos++;
			}
			while(*pos != '\0' && isspace(*pos) == 0)	{
				pos++;
			}
			intensities.push_back(atof(pos));
		}
		else if(temp.find("END IONS") != temp.npos)	{
			sp.clear();
			if(parent*charge > 600.0)	{
				sp.pm = (int64_t)(0.5 + 1000*(parent*charge - proton*charge));
				sp.pz = (int64_t)charge;
				sp.pi = 100;
				sp.pt = pt;
				sp.desc = desc;
				if(scan == 0)	{
					sp.sc = s;
				}
				else	{
					sp.sc = scan;
				}
				size_t i = 0;
				pair<int64_t,int64_t> p;
				while(i < masses.size())	{
					p.first = (int64_t)(0.5+1000.0*(masses[i]-proton));
					p.second = (int64_t)(0.5+intensities[i]);
					sp.mis.push_back(p);
					i++;
				}
				sp.condition(res,50);
				spectra.push_back(sp);
				if(false)	{
					auto it = sp.mis.begin();
					while(it != sp.mis.end())	{
						cout << it->first << "\t" << it->second << endl;
						it++;
					}
				}
			}
			desc = "";
			masses.clear();
			intensities.clear();
			parent = 0.0;
			charge = 0.0;
			s++;
			if(s % 2500 == 0)	{
				cout << '.';
				cout.flush();
			}
			if(s != 0 and s % 50000 == 0)	{
				cout << ' ' << s << endl;
				cout.flush();
			}
			scan = 0;
		}
	}
	cout << "\n";
	cout.flush();
	delete[] line;
	istr.close();
	return true;
}

