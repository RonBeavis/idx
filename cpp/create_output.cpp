/*
#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
# Loads a spectrum file into a vector of spectrum objects
#
*/
#include "pch.h"

#include "rapidjson/document.h"
#include <fstream>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <sys/stat.h>
#include <unordered_map>
#include <map>
#include <string>
#include <set>
#include <vector>
using namespace std;
#include "load_spectra.hpp"
#include "load_kernel.hpp"
#include "create_results.hpp"
#include "create_output.hpp"
//#include <boost/math/distributions/hypergeometric.hpp>
//#include <boost/math/policies/policy.hpp>

create_output::create_output(void)	{
	low = -20.0;
	high = 20.0;
	load_distribution();
}

create_output::~create_output(void)	{
}

long create_output::get_cells(double _pm,long _res)	{
	long pm = 100*(long)(_pm/100000);
	long r = 2;
	if(_res == 50)	{
		r = 1;
	}
	else if(_res == 20)	{
		r = 0;
	}
	if(pm < 800)	{
		return (long)(pm*distribution[800][r]);
	}
	if(pm > 4000)	{
		return (long)(pm*distribution[4000][r]);
	}
	return (long)(pm*distribution[pm][r]);
}

bool create_output::apply_model(long _res,string& _seq,double _pm,long _ions,long _lspectrum,double& pscore,double& p){
	p = 0.0001;
	pscore = 400;
	long sfactor = 20;
	long sadjust = 3;
	if(_res > 100)	{
		sfactor = 40;
	}
	long cells = get_cells(_pm,_res);
	long total_ions = 2*(long)(_seq.size() - 1);
	if(total_ions > sfactor)	{
		total_ions = sfactor;
	}
	if(total_ions < _ions)	{
		total_ions = _ions + 1;
	}
	long sc = _lspectrum * sadjust;
	if(_ions >= sc)	{
		sc = _ions + 2;
	}
//	boost::math::hypergeometric_distribution<double> hyper(sc,total_ions,cells);
//	p = boost::math::pdf<double>(hyper, _ions);
	hypergeom hp(sc,total_ions,cells);
	p = hp.pdf(_ions);
	pscore = -100.0*log(p)/2.3;
	return true;
}

bool create_output::load_mods(void)	{
	ifstream istr;
	istr.open("report_mods.txt");
	if(istr.fail())	{
		return false;
	}
	mt.clear();
	string line;
	while(getline(istr,line))	{
		size_t tab = line.find('\t');
		mt[atoi(line.substr(0,tab).c_str())] = line.substr(tab+1,line.size()-1);
	}
	istr.close();
	return true;
}

bool create_output::find_window(void)	{
	map<long,long> vs;
	long i = 0;
	for(i = -20; i < 21; i++)	{
		vs[i] = 0;
	}
	auto it = ppms.begin();
	long max = -21;
	long center = -1;
	while(it != ppms.end())	{
		i = (long)(0.5 + it->second);
		if(i >= -20 and i <= 20)	{
			vs[i] += 1;
		}
		if(vs[i] > max)	{
			max = vs[i];
			center = i;
		}
		it++;
	}
	double ic = (double)max;
	if(ic < 100.0)	{
		low = -20.0;
		high = 20.0;
		return true;
	}
	long l = -20;
	for(i = -20;i < center; i++)	{
		if(vs[i]/ic >= 0.01)	{
			l = i;
			break;
		}
	}
	long h = 20;
	for(i = 20; i > center; i--)	{
		if(vs[i]/ic >= 0.01)	{
			h = i;
			break;
		}
	}
	low = (double)l;
	high = (double)h;
	return true;
}

bool create_output::create(map<string,string>& _params,create_results& _cr)	{
	ifstream istr;
	istr.open(_params["kernel file"]);
	if(istr.fail())	{
		return false;
	}
	if(!load_mods())	{
		cout << "Warning (idx:1001): annotation file \"report_mods.txt\" was not present" << endl;
	}
	long k = 0;
	for(size_t j = 0; j < _cr.ids.size(); j++)	{
		sv[_cr.ids[j].sn] = _cr.ids[j];
		for(size_t i = 0; i < _cr.ids[j].ks.size(); i++)	{
			k = _cr.ids[j].ks[i];
			if(sdict.find(k) != sdict.end())	{
				sdict[k].insert(_cr.ids[j].sn);
			}
			else	{
				sdict.insert(pair<long,set<long> >(k,set<long>()));
				sdict[k].insert(_cr.ids[j].sn);
			}
		}
	}

	long res = atoi(_params["fragment tolerance"].c_str());
	long inferred = 0;
	long specs = atoi(_params["spectra"].c_str());
	double score_min = 200.0;
	if(specs > 0)	{
		score_min += 100.0 * log(specs)/2.3;
	}
	double total_prob = 0.0;
	long min_c = 8;
	if(res == 50)	{
		min_c = 7;
	}
	else if(res == 20)	{
		min_c = 6;
	}
	string line;
	using namespace rapidjson;
	long c = 0;
	long h = 0;
	long s_count = 0;
	while(getline(istr,line))	{
		if(c != 0 and c % 10000 == 0)	{
			cout << '.';
			cout.flush();
		}
		if(c != 0 and c % 200000 == 0)	{
			cout << " " << c << endl;
			cout.flush();
		}
		c++;
 		Document js;
   		js.Parse(line.c_str(),line.length());
		if(!js.HasMember("pm"))	{
			continue;
		}
		h = (long)js["h"].GetInt();
		if(sdict.find(h) == sdict.end())	{
			continue;
		}
		if(h != (long)js["u"].GetInt())	{
			inferred += 1;
		}
		long last_i = 0;

		double max_prob = 0.0;
		double prob = 0.0;
		double delta = 0.0;
		double pm = 0.0;
		double ppm = 0.0;
		double score = 0.0;
		auto it = sdict[h].begin();
		long s = 0;
		string seq;
		while(it != sdict[h].end())	{
			s = *it;
			it++;
			pm = js["pm"].GetFloat();
			delta = (sv[s].pm-pm)/1000.0;
			ppm = 1.0e6*(sv[s].pm-pm)/pm;
			if(delta > 0.5)	{
				ppm = 1.0e6*(sv[s].pm-c13-pm)/pm;
			}
			seq = js["seq"].GetString();
			apply_model(res,seq,pm,sv[s].peaks,sv[s].ions,score,prob);
			if(score < score_min or sv[s].ri < 0.20 or sv[s].peaks < min_c)	{
				continue;
			}
			if(prob > max_prob)	{
				max_prob = prob;
			}
			ostringstream oline;
			oline << sv[s].sc << "\t" << fixed << setprecision(3) << pm/1000.0 << "\t";
			oline << delta << "\t" << setprecision(1) << ppm << setprecision(3) << "\t";
			oline << sv[s].pz << "\t" << js["lb"].GetString() << "\t";
			oline << js["beg"].GetInt() << "\t" << js["end"].GetInt() << "\t";
			oline << js["pre"].GetString() << "\t" << js["seq"].GetString() << "\t";
			oline << js["post"].GetString() << "\t" << sv[s].peaks << "\t";
			const Value& jbs = js["ns"];
			long lns = 0;
			for(SizeType a = 0; a < jbs.Size();a++)	{
				lns += jbs[a].GetInt();
			}
			oline << fixed << setprecision(2);
			if(lns > 0)	{
				oline << sv[s].ri << "\t" << setprecision(1) << log(lns)/2.3 << "\t" << -0.01*score << "\t";
			}
			else	{
				oline << sv[s].ri << "\t" << setprecision(1) << "-" << "\t" << -0.01*score << "\t";
			}
			if(js.HasMember("mods"))	{
				const Value& jmods = js["mods"];
				vector<mod> mods;
				mod tmod;
				for(SizeType a = 0; a < jmods.Size();a++)	{
					const Value& lmods = jmods[a];
					tmod.pos = lmods[1].GetInt();
					tmod.res = lmods[0].GetString();
					tmod.mass = lmods[2].GetInt();
					mods.push_back(tmod);
				}
				sort(mods.begin(),mods.end());
				for(size_t a = 0; a < mods.size(); a++)	{
					if(mt.find(mods[a].mass) != mt.end())	{
						oline << mods[a].res << mods[a].pos << "~" << mt.find(mods[a].mass)->second << ";";
					}
					else	{
						oline << mods[a].res << mods[a].pos << "#" << mods[a].mass/1000.0 << ";";
					}
				}
			}
			last_i = s+1;
			if(odict.find(s) != odict.end())	{
				odict[s].push_back(oline.str());			
			}
			else	{
				odict.insert(pair<long,vector<string> >(s,vector<string>()));
				odict[s].push_back(oline.str());
			}
			ppms[s] = (long)(0.5+ppm);
		}
		total_prob += max_prob;
		s_count++;

	}
	ofstream ofs;
	ofs.open(_params["output file"]);
	string header = "Id\tSub\tScan\tPeptide mass\tDelta\tppm\tz\tProtein acc\t";
	header += "Start\tEnd\tPre\tSequence\tPost\tIC\tRI\tlog(f)\tlog(p)\tModifications";
	ofs << header << endl;
	find_window();
/*	long ids = 0;
	for(long s = 0; s < odict.size(); s++)	{
		for(long l = 0; l < odict[s].size(); l++)	{
			ofs << ids+1 << "\t" << l+1 << "\t" << odict[s][l] << endl;
		}
		if(!odict[s].empty())	{
			ids++;
		}
	}	*/

	long err = 0;
	long sub = 0;
	long tot = 0;
	string t;
	for(long a = 0; a < (long)odict.size(); a++)	{
		sub = 1;
		for(size_t b = 0; b < odict[a].size(); b++)	{
			t = odict[a][b];
			double ps = get_ppm(t);
			if(ps <= high and ps >= low)	{
				ofs << a << "\t" << sub << "\t" << t << endl;
				sub++;
				tot++;
			}
			else	{
				err++;
			}
		}
	}
	cout << "\n     lines = " << tot << endl;
	if(total_prob > 0)	{
		cout << "     fpr = " << scientific << setprecision(1) << total_prob << endl;
	}
	if(low != -20 and high != 20)	{
		double ble = 100.0*((high-low)/41.0)*(double)err/(double)tot;
		cout << "     baseline error = " << fixed << setprecision(1) << ble << "%" << endl;
	}
	else	{
		cout << "     baseline error = n/a" << endl;
	}
	cout << "     parent ion tolerance = " << low << "," << high << endl;


	ofs.close();
	cout << endl;
	cout.flush();
	istr.close();
	return true;
}

