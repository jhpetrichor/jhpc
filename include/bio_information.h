#ifndef BIO_INFORMATION_H
#define BIO_INFORMATION_H

#include "config.h"
#include <map>
#include <string>
#include <vector>
#include <set>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <algorithm>

using namespace std;

class BioInformation {
public:
    map<string, vector<double>> go_expression;
    map<string, set<string>> go_slim;
    map<string, map<string, double>> subcelluar;
public:
    BioInformation();
    void read_gene_expression();
    void read_subcellular();
    void read_go_slim();
    double go_expression_pearson(const string&, const string&);
    double go_expression_spearman(const string&, const string&);
    double go_slim_slim(const string&, const string&);
    double go_slim_slim_jaccard(const string&, const string&);
    void print_go_slim();
};

#endif 
