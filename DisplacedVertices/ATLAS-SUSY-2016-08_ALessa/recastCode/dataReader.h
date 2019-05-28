#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include <ctime>
#include <boost/tokenizer.hpp>
#include <fstream>
#include <stdio.h>
#include <string.h>

using namespace std;

vector< vector<double> > getData(string filename)
{

    ifstream in(filename.c_str());

    typedef boost::tokenizer< boost::escaped_list_separator<char> > Tokenizer;

    vector< string > vec;
    vector< vector<double> > pts;
    string line;

    while (getline(in,line))
    {

    	//Skip comment lines
    	if (strncmp(line.c_str(),"#",1) == 0){continue;}
    	//Skip header lines
    	if (strncmp(line.c_str(),"'",1) == 0){continue;}

    	Tokenizer tok(line);
        vec.assign(tok.begin(),tok.end());
        vector<double> pt;

        if (vec.size() < 2){continue;}
        for (int i=0; i < vec.size(); ++i){
        	double x = atof(vec[i].c_str());
        	pt.push_back(x);
        }
        pts.push_back(pt);
    }

    return pts;
}



//Load all efficiency tables
//Event selection:
vector< vector<double> > effEv_table22 = getData("recastCode/ATLAS_data/HEPData-ins1630632-v2-csv/Table22.csv");
vector< vector<double> > effEv_table23 = getData("recastCode/ATLAS_data/HEPData-ins1630632-v2-csv/Table23.csv");
vector< vector<double> > effEv_table24 = getData("recastCode/ATLAS_data/HEPData-ins1630632-v2-csv/Table24.csv");
//DV reconstruction:
vector< vector<double> > effDV_table25 = getData("recastCode/ATLAS_data/HEPData-ins1630632-v2-csv/Table25.csv");
vector< vector<double> > effDV_table26 = getData("recastCode/ATLAS_data/HEPData-ins1630632-v2-csv/Table26.csv");
vector< vector<double> > effDV_table27 = getData("recastCode/ATLAS_data/HEPData-ins1630632-v2-csv/Table27.csv");
vector< vector<double> > effDV_table28 = getData("recastCode/ATLAS_data/HEPData-ins1630632-v2-csv/Table28.csv");
vector< vector<double> > effDV_table29 = getData("recastCode/ATLAS_data/HEPData-ins1630632-v2-csv/Table29.csv");
vector< vector<double> > effDV_table30 = getData("recastCode/ATLAS_data/HEPData-ins1630632-v2-csv/Table30.csv");
vector< vector<double> > effDV_table31 = getData("recastCode/ATLAS_data/HEPData-ins1630632-v2-csv/Table31.csv");
vector< vector<double> > effDV_table32 = getData("recastCode/ATLAS_data/HEPData-ins1630632-v2-csv/Table32.csv");
vector< vector<double> > effDV_table33 = getData("recastCode/ATLAS_data/HEPData-ins1630632-v2-csv/Table33.csv");
vector< vector<double> > effDV_table34 = getData("recastCode/ATLAS_data/HEPData-ins1630632-v2-csv/Table34.csv");
vector< vector<double> > effDV_table35 = getData("recastCode/ATLAS_data/HEPData-ins1630632-v2-csv/Table35.csv");
vector< vector<double> > effDV_table36 = getData("recastCode/ATLAS_data/HEPData-ins1630632-v2-csv/Table36.csv");


