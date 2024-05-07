/*
 * @Description: 
 * @Author: jh
 * @Date: 2024-04-30 17:11:48
 * @LastEditTime: 2024-05-04 21:27:37
 */
#include "gene_express.h"
#include "ungraph.h"

#include <cmath>

#include <fstream>
#include <istream>
#include <sstream>
#include <algorithm>
#include <utility>

GeneExpress::GeneExpress(string file_path) {
    read_gene_express(file_path);
}

GeneExpress::~GeneExpress() {
    gene_express.clear();
}

void GeneExpress::read_gene_express(string& file_path) {
    fstream file(file_path);
    if(!file.is_open()) {
        cerr << "<ERROR> Failed to open file! " << file_path << endl;
        exit(1); 
    }
    string line;
    while(getline(file, line)) {
        istringstream iss(line);
        string protein_name;
        iss >> protein_name;
        vector<double> temp_express(36, 0.0);
        for(int i = 0; i < 36; ++i) {
            iss >> temp_express[i];
        }
        gene_express.insert(move(make_pair(move(protein_name), move(temp_express))));
    }

    file.close();
}

map<string, double> GeneExpress::activate_by_three_sigma() {
    map<string, double> result;
    // 计算mean and 。。。
    for(auto& it: gene_express) {
        // update average gene expression
        double temp_sum = 0.0;
        for(auto& express: it.second) {
            temp_sum += express;
        }
        double mean = temp_sum / it.second.size();

        // update gene expression variance
        double temp_variance = 0.0;
        for(auto& express: it.second) {
            temp_variance += pow(express - mean, 2);
        }
        double variance = temp_variance / (it.second.size() - 1);
        // cout << "variance: " << variance << endl;
        // update activate threshold
        double active_threshold = active_three_sigma(mean, variance);
        result.insert(make_pair(it.first, active_threshold));
    }
    return move(result);
}

double GeneExpress::active_three_sigma(double mean, double varience) {
    const double f = 1.0 / (1.0 + varience);
    const double s_varience = pow(varience, 0.5);
    return mean + 3  * s_varience * (1.0 - f);
}

// need todo()!
map<string, double> GeneExpress::active_by_top() {
    map<string, double> result;
    for(auto& it: gene_express) {
        vector<double> temp_expression(it.second.begin(), it.second.end());
        sort(temp_expression.begin(), temp_expression.end());
        // return temp_expression[10];
        result.insert(make_pair(it.first, temp_expression[12]));
    }
    return move(result);
}

vector<UnGraph> GeneExpress::build_dynamic_PPI(const UnGraph* g, DPIN_MEHTOD method) {
    switch(method) {
        case DPIN_MEHTOD::THREE_SIGMA: {
            auto active = activate_by_three_sigma();
            return move(build_dynamic_PPI_by_active(g, active, 36));
        }
        case DPIN_MEHTOD::TOP: {
            auto active = active_by_top();
            return move(build_dynamic_PPI_by_active(g, active, 36));
        }
        case DPIN_MEHTOD::TIME: {
            auto active = activate_by_three_sigma();
            return move(build_dynamic_PPI_by_time(g, active, 12));
        }
    }
}

vector<UnGraph> GeneExpress::build_dynamic_PPI_by_active(const UnGraph* g, map<string, double>& active, int count) {
    std:vector<UnGraph> dpins(count);
    // update proteins
    vector<set<string>> set_proteins(count);
    vector<vector<string>> list_edges(count);
    vector<vector<double>> edge_weight(count);
    for(auto& protein: g->ID2Protein) {
        auto it = gene_express.find(protein->protein_name);
        // 没有这个基因的表达信息，则在所有自网络中保留
        if(it == gene_express.end()) {
            for(int i = 0; i < count; ++i) {
                set_proteins[i].insert(protein->protein_name);
            }
            continue;
        }
        // 存在基因的表达信息，则需要筛选，满足阈值则保留
        for(int i = 0; i < count; ++i) {
            double temp_express = (it->second)[i];
            if( temp_express > active[it->first]) {
                set_proteins[i].insert(it->first);
            }
         }
    }

    // update edges 
    for(auto& e: g->edges) {
        for(int i = 0; i < count; ++i) {
            if(set_proteins[i].count(e->node_a->protein_name) && set_proteins[i].count(e->node_b->protein_name)) {
                list_edges[i].emplace_back(e->node_a->protein_name);
                list_edges[i].emplace_back(e->node_b->protein_name);
                edge_weight[i].emplace_back(e->balanced_weight);
            }
        }
    }
    for(int i = 0; i < count; ++i) {
        cout << set_proteins[i].size() << "\t" << list_edges[i].size() << endl;
        UnGraph temp_graph(move(set_proteins[i]), move(list_edges[i]), move(edge_weight[i]));
        dpins[i] = temp_graph;
    }

    return move(dpins);
}

vector<UnGraph> GeneExpress::build_dynamic_PPI_by_time(const UnGraph* g, map<string, double>& active, int count) {
    vector<UnGraph> dpins(12);
    vector<set<string>> set_proteins(count);
    vector<vector<string>> list_edges(count);
    vector<vector<double>> edge_weight(count);
    for(auto& protein: g->ID2Protein) {
        auto it = gene_express.find(protein->protein_name);
        // 没有这个基因的表达信息，则在所有自网络中保留
        if(it == gene_express.end()) {
            for(int i = 0; i < count; ++i) {
                set_proteins[i].insert(protein->protein_name);
            }
            continue;
        }
        // 存在基因的表达信息，则需要筛选，满足阈值则保留
        for(int i = 0; i < count; ++i) {
            vector<double> temp_expression{it->second[i], it->second[i + 1], it->second[i + 2]};
            auto max = std::max(temp_expression.begin(), temp_expression.end());
            if(it->second[i] > active[it->first]) {
                set_proteins[i].insert(it->first);
            }
         }
    }

    // update edges 
    for(auto& e: g->edges) {
        // if(set_proteins.count(e->node_a->protein_name) && set_proteins.count(e->node_b->protein_name)) {

        // }
        for(int i = 0; i < count; ++i) {
            if(set_proteins[i].count(e->node_a->protein_name) && set_proteins[i].count(e->node_b->protein_name)) {
                list_edges[i].emplace_back(e->node_a->protein_name);
                list_edges[i].emplace_back(e->node_b->protein_name);
                edge_weight[i].emplace_back(e->balanced_weight);
            }
        }
    }
    for(int i = 0; i < count; ++i) {
        std::cout << set_proteins.size() << "\t" << list_edges.size() << endl;
        UnGraph temp_graph(move(set_proteins[i]), move(list_edges[i]), move(edge_weight[i]));
        dpins[i] = temp_graph;
    }
    return move(dpins);
}
