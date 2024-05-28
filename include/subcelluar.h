#ifndef __SUBCELLULAR_H__
#define __SUBCELLULAR_H__

#include "graph.h"
#include "ungraph.h"

#include <regex>
#include <vector>

class Subcellular {
public:
    map<string, set<string>> location_protein;
    Subcellular(string& file);
    vector<UnGraph> build_spin(UnGraph& g);
private:
    void read_subcelluar_file_human(string& file);
    void read_subcelluar_file(string& file);
    set<string> extract_go_terms(const std::string& data);
};

#endif