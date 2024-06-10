#include "../../include/subcelluar.h"
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <utility>

Subcellular::Subcellular(string& subcellular_file) {
    // read_subcelluar_file(subcellular_file);
    read_subcelluar_file_human(subcellular_file);
}

void Subcellular::read_subcelluar_file_human(string& file_path) {
    std::fstream file(file_path);
    if (!file.is_open()) {
        cerr << "<ERROR> Failed to open subcellular file! " << file_path
             << endl;
        exit(1);
    }

    string line;
    while (getline(file, line)) {
        auto gos = extract_go_terms(line);
        istringstream is(line);
        string protein;
        is >> protein >> protein;
        for (auto& go : gos) {
            location_protein[go].insert(protein);
        }
    }

    file.close();
}

void Subcellular::read_subcelluar_file(string& subcellular_file) {
    ifstream file(subcellular_file);
    if (!file) {
        cerr << "Failed to open file! " << subcellular_file << endl;
        exit(1);
    }

    string line;
    while (getline(file, line)) {
        istringstream iss(line);
        string protein1, protein2, location;
        iss >> protein1 >> protein2 >> location;
        location_protein[location].insert(protein1);
        location_protein[location].insert(protein2);
    }

    file.close();
}

set<string> Subcellular::extract_go_terms(const std::string& data) {
    set<std::string> goTerms;
    regex pattern("\\bGO:(\\d+)\\b"); // 正则表达式
    smatch matches;

    // 使用迭代器遍历每一行数据
    string::const_iterator start = data.begin();
    string::const_iterator end = data.end();

    while (regex_search(start, end, matches, pattern)) {
        goTerms.insert("GO:" + matches.str(1));
        start = matches.suffix().first;
    }
    return goTerms;
}

vector<UnGraph> Subcellular::build_spin(UnGraph& g) {
    // vector<set<string>> protein_set;
    // vector<vector<string>> edge_list;
    // vector<vector<double>> weight;
    vector<UnGraph> spins;

    for (auto& gl : location_protein) {
        // string location = gl.first;
        set<string> protein_set;
        vector<string> edge_list;
        vector<double> edge_weight;
        // add proteins
        for (auto& p : gl.second) {
            if (g.proteins_name.count(p)) {
                protein_set.insert(p);
            }
        }
        // add edges
        for (auto& e : g.edges) {
            if (protein_set.count(e->node_a->protein_name) &&
                protein_set.count(e->node_b->protein_name)) {
                edge_list.emplace_back(e->node_a->protein_name);
                edge_list.emplace_back(e->node_b->protein_name);
                edge_weight.emplace_back(e->balanced_weight);
            }
        }

        spins.emplace_back(UnGraph(std::move(protein_set), std::move(edge_list),
                                   std::move(edge_weight)));
    }

    return spins;
}
