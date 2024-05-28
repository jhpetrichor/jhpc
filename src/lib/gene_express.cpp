#include "../../include/gene_express.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <istream>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

GeneExpress::GeneExpress(string file_path) {
    read_gene_express(file_path);
    string essential_protein_file = ESSENTIAL_PROTEIN;
    read_essential_protein(essential_protein_file);
}

GeneExpress::~GeneExpress() { gene_express.clear(); }

void GeneExpress::read_gene_express(string& file_path) {
    fstream file(file_path);
    if (!file.is_open()) {
        cerr << "<ERROR> Failed to open file! " << file_path << endl;
        exit(1);
    }
    string line;
    while (getline(file, line)) {
        istringstream iss(line);
        string protein_name;
        iss >> protein_name;
        vector<double> temp_express(36, 0.0);
        for (int i = 0; i < 36; ++i) {
            iss >> temp_express[i];
        }
        gene_express.insert(
            move(make_pair(move(protein_name), move(temp_express))));
    }

    file.close();
}

void GeneExpress::read_essential_protein(string& file_path) {
    fstream file(file_path);
    if (!file.is_open()) {
        cerr << "<ERROR>: Failed to open the file of essential protein! "
             << file_path << endl;
        exit(1);
    }

    string line;
    while (getline(file, line)) {
        istringstream iss(line);
        string protein;
        iss >> protein;
        essential_proteins.insert(protein);
    }
    file.close();
}

map<string, double> GeneExpress::activate_by_three_sigma() {
    map<string, double> result;
    // 计算mean and 。。。
    for (auto& it : gene_express) {
        // update average gene expression
        double temp_sum = 0.0;
        for (auto& express : it.second) {
            temp_sum += express;
        }
        double mean = temp_sum / it.second.size();

        // update gene expression variance
        double temp_variance = 0.0;
        for (auto& express : it.second) {
            temp_variance += pow(express - mean, 2);
        }
        double variance = temp_variance / (it.second.size() - 1);
        // cout << "variance: " << variance << endl;
        // update activate threshold
        double active_threshold = active_three_sigma(mean, variance);
        result.insert(make_pair(it.first, active_threshold));
    }
    return result;
}

double GeneExpress::active_three_sigma(double mean, double varience) {
    const double f = 1.0 / (1.0 + varience);
    const double s_varience = pow(varience, 0.5);
    return mean + 3 * s_varience * (1.0 - f);
}

// need todo()!
map<string, double> GeneExpress::active_by_top() {
    map<string, double> result;
    for (auto& it : gene_express) {
        vector<double> temp_expression(it.second.begin(), it.second.end());
        sort(temp_expression.begin(), temp_expression.end());
        // return temp_expression[10];
        result.insert(make_pair(it.first, temp_expression[12]));
    }
    return result;
}

vector<UnGraph> GeneExpress::build_dynamic_PPI(const UnGraph* g,
                                               DPIN_MEHTOD method) {
    switch (method) {
    case DPIN_MEHTOD::THREE_SIGMA: {
        auto active = activate_by_three_sigma();
        return build_dynamic_PPI_by_active(g, active, 36);
    }
    case DPIN_MEHTOD::TOP: {
        auto active = active_by_top();
        return build_dynamic_PPI_by_active(g, active, 36);
    }
    case DPIN_MEHTOD::TIME: {
        auto active = activate_by_three_sigma();
        return build_dynamic_PPI_by_time(g, active, 12);
    }
    case DPIN_MEHTOD::ESSENTIAL: {
        auto active = active_by_top();
        return build_dynamic_PPI_by_essential_protein(g, active, 36);
    }
    }
}

vector<UnGraph> GeneExpress::build_dynamic_PPI_by_active(
    const UnGraph* g, map<string, double>& active, int count) {
std:
    vector<UnGraph> dpins(count);
    // update proteins
    vector<set<string>> set_proteins(count);
    vector<vector<string>> list_edges(count);
    vector<vector<double>> edge_weight(count);
    for (auto& protein : g->ID2Protein) {
        auto it = gene_express.find(protein->protein_name);
        // 没有这个基因的表达信息，则在所有自网络中保留
        if (it == gene_express.end()) {
            for (int i = 0; i < count; ++i) {
                set_proteins[i].insert(protein->protein_name);
            }
            continue;
        }
        // 存在基因的表达信息，则需要筛选，满足阈值则保留
        for (int i = 0; i < count; ++i) {
            double temp_express = (it->second)[i];
            if (temp_express > active[it->first]) {
                set_proteins[i].insert(it->first);
            }
        }
    }

    // update edges
    for (auto& e : g->edges) {
        for (int i = 0; i < count; ++i) {
            if (set_proteins[i].count(e->node_a->protein_name) &&
                set_proteins[i].count(e->node_b->protein_name)) {
                list_edges[i].emplace_back(e->node_a->protein_name);
                list_edges[i].emplace_back(e->node_b->protein_name);
                edge_weight[i].emplace_back(e->balanced_weight);
            }
        }
    }
    for (int i = 0; i < count; ++i) {
        cout << set_proteins[i].size() << "\t" << list_edges[i].size() / 2
             << endl;
        UnGraph temp_graph(move(set_proteins[i]), move(list_edges[i]),
                           move(edge_weight[i]));
        dpins[i] = temp_graph;
    }

    return dpins;
}

vector<UnGraph>
GeneExpress::build_dynamic_PPI_by_time(const UnGraph* g,
                                       map<string, double>& active, int count) {
    vector<UnGraph> dpins(count);
    vector<set<string>> set_proteins(count);
    vector<vector<string>> list_edges(count);
    vector<vector<double>> edge_weight(count);
    for (auto& protein : g->ID2Protein) {
        auto it = gene_express.find(protein->protein_name);
        // 没有这个基因的表达信息，则在所有自网络中保留
        if (it == gene_express.end()) {
            for (int i = 0; i < count; ++i) {
                set_proteins[i].insert(protein->protein_name);
            }
            continue;
        }
        // 存在基因的表达信息，则需要筛选，满足阈值则保留
        for (int i = 0; i < count; ++i) {
            vector<double> temp_expression{it->second[i], it->second[i + 1],
                                           it->second[i + 2]};
            auto max = std::max(temp_expression.begin(), temp_expression.end());
            if (it->second[i] > active[it->first]) {
                set_proteins[i].insert(it->first);
            }
        }
    }

    // update edges
    for (auto& e : g->edges) {
        // if(set_proteins.count(e->node_a->protein_name) &&
        // set_proteins.count(e->node_b->protein_name)) {

        // }
        for (int i = 0; i < count; ++i) {
            if (set_proteins[i].count(e->node_a->protein_name) &&
                set_proteins[i].count(e->node_b->protein_name)) {
                list_edges[i].emplace_back(e->node_a->protein_name);
                list_edges[i].emplace_back(e->node_b->protein_name);
                edge_weight[i].emplace_back(e->balanced_weight);
            }
        }
    }
    for (int i = 0; i < count; ++i) {
        std::cout << set_proteins[i].size() << "\t" << list_edges[i].size()
                  << endl;
        UnGraph temp_graph(move(set_proteins[i]), move(list_edges[i]),
                           move(edge_weight[i]));
        dpins[i] = temp_graph;
    }
    return dpins;
}

vector<UnGraph> GeneExpress::build_dynamic_PPI_by_essential_protein(
    const UnGraph* g, map<string, double>& active, int count) {
    vector<UnGraph> dpins(count);
    vector<set<string>> set_proteins(count);
    vector<vector<string>> list_edges(count);
    vector<vector<double>> edge_weight(count);
    for (auto& protein : g->ID2Protein) {
        auto it = gene_express.find(protein->protein_name);
        // 1、 没有这个基因的表达信息
        // 2、 若当前蛋白质为关键蛋白质
        // 则在所有自网络中保留
        if (it == gene_express.end() || essential_proteins.count(it->first)) {
            for (int i = 0; i < count; ++i) {
                set_proteins[i].insert(protein->protein_name);
            }
            continue;
        }
        // 存在基因的表达信息，则需要筛选，满足阈值则保留
        for (int i = 0; i < count; ++i) {
            vector<double> temp_expression{it->second[i], it->second[i + 1],
                                           it->second[i + 2]};
            auto max = std::max(temp_expression.begin(), temp_expression.end());
            if (it->second[i] > active[it->first]) {
                set_proteins[i].insert(it->first);
            }
        }
    }

    // update edges
    for (auto& e : g->edges) {
        // if(set_proteins.count(e->node_a->protein_name) &&
        // set_proteins.count(e->node_b->protein_name)) {

        // }
        for (int i = 0; i < count; ++i) {
            if (set_proteins[i].count(e->node_a->protein_name) &&
                set_proteins[i].count(e->node_b->protein_name)) {
                list_edges[i].emplace_back(e->node_a->protein_name);
                list_edges[i].emplace_back(e->node_b->protein_name);
                edge_weight[i].emplace_back(e->balanced_weight);
            }
        }
    }
    for (int i = 0; i < count; ++i) {
        std::cout << set_proteins[i].size() << "\t" << list_edges[i].size() / 2
                  << endl;
        UnGraph temp_graph(move(set_proteins[i]), move(list_edges[i]),
                           move(edge_weight[i]));
        dpins[i] = temp_graph;
    }
    return dpins;
}

pair<double, double> GeneExpress::calculate_mean_varience(vector<double>& arr) {
    double sum = accumulate(arr.begin(), arr.end(), 0.0);
    const double mean = sum / static_cast<double>(arr.size());
    double sum_varience = 0.0;
    for (auto value : arr) {
        sum_varience += pow(mean - value, 2);
    }
    const double varience = sum_varience / static_cast<double>(arr.size());
    return make_pair(mean, varience);
}

vector<UnGraph> GeneExpress::build_KPIN(const UnGraph* g) {
    int count = 36;
    map<string, double> active;
    for (auto& p : g->proteins) {
        // 查找基因表达
        double active_temp = 0.0;
        auto it = gene_express.find(p->protein_name);
        if (it != gene_express.end()) {
            auto mean_varience = calculate_mean_varience(it->second);
            if (essential_proteins.count(p->protein_name)) {
                // 是关键蛋白质，基因表达活性阈值更低一些，以保证能在生命周期中
                // 的大多数时间内保持活性
                active_temp = mean_varience.first;
            } else { // 非关键蛋白质，具有更高的基因表达活性表达阈值
                active_temp =
                    mean_varience.first - pow(mean_varience.second, 0.5);
            }
        } else {
            active_temp = 0.0;
        }
        active.insert(make_pair(p->protein_name, active_temp));
    }

    // 构建动态网络
    vector<UnGraph> dpins(count);
    vector<set<string>> set_proteins(count);
    vector<vector<string>> list_edges(count);
    vector<vector<double>> edge_weight(count);
    for (auto& protein : g->proteins_name) {
        // 是够含有当前蛋白质的基因表达数据？
        auto it = gene_express.find(protein);
        if (it == gene_express.end()) {
            for (int i = 0; i < count; ++i) {
                set_proteins[i].insert(protein);
            }
            continue;
        }
        for (int i = 0; i < count; ++i) {
            if (it->second[i] > active[protein]) {
                set_proteins[i].insert(protein);
            }
        }
    }

    // update edges
    for (auto& e : g->edges) {
        for (int i = 0; i < count; ++i) {
            if (set_proteins[i].count(e->node_a->protein_name) &&
                set_proteins[i].count(e->node_b->protein_name)) {
                list_edges[i].emplace_back(e->node_a->protein_name);
                list_edges[i].emplace_back(e->node_b->protein_name);
                edge_weight[i].emplace_back(e->balanced_weight);
            }
        }
    }

    for (int i = 0; i < count; ++i) {
        std::cout << set_proteins[i].size() << "\t" << list_edges[i].size() / 2
                  << endl;
        UnGraph temp_graph(move(set_proteins[i]), move(list_edges[i]),
                           move(edge_weight[i]));
        dpins[i] = temp_graph;
    }
    return dpins;
}

/**
 * @brief: 构建动态网络，使用
 * @param {UnGraph*} g
 * @return {*}
 */
vector<UnGraph> GeneExpress::build_KPIN2(const UnGraph* g) {
    int count = 12;
    map<string, double> active;
    for (auto& p : g->proteins) {
        // 查找基因表达
        double active_temp = 0.0;
        auto it = gene_express.find(p->protein_name);
        if (it != gene_express.end()) {
            auto mean_varience = calculate_mean_varience(it->second);
            if (essential_proteins.count(p->protein_name)) {
                // 是关键蛋白质，基因表达活性阈值更低一些，以保证能在生命周期中
                // 的大多数时间内保持活性
                active_temp = mean_varience.first;
            } else { // 非关键蛋白质，具有更高的基因表达活性表达阈值
                active_temp =
                    mean_varience.first - pow(mean_varience.second, 0.5);
            }
        } else {
            active_temp = 0.0;
        }
        active.insert(make_pair(p->protein_name, active_temp));
    }

    // 构建动态网络
    vector<UnGraph> dpins(count);
    vector<set<string>> set_proteins(count);
    vector<vector<string>> list_edges(count);
    vector<vector<double>> edge_weight(count);
    for (auto& protein : g->proteins_name) {
        // 是够含有当前蛋白质的基因表达数据？
        auto it = gene_express.find(protein);
        if (it == gene_express.end()) {
            for (int i = 0; i < count; ++i) {
                set_proteins[i].insert(protein);
            }
            continue;
        }
        for (int i = 0; i < count; ++i) {
            double exp_temp =
                (it->second[i] + it->second[12 + i] + it->second[24 + i]) / 3;
            if (exp_temp > active[protein]) {
                set_proteins[i].insert(protein);
            }
        }
    }

    // update edges
    for (auto& e : g->edges) {
        for (int i = 0; i < count; ++i) {
            if (set_proteins[i].count(e->node_a->protein_name) &&
                set_proteins[i].count(e->node_b->protein_name)) {
                list_edges[i].emplace_back(e->node_a->protein_name);
                list_edges[i].emplace_back(e->node_b->protein_name);
                edge_weight[i].emplace_back(e->balanced_weight);
            }
        }
    }

    for (int i = 0; i < count; ++i) {
        std::cout << set_proteins[i].size() << "\t" << list_edges[i].size() / 2
                  << endl;
        UnGraph temp_graph(move(set_proteins[i]), move(list_edges[i]),
                           move(edge_weight[i]));
        dpins[i] = temp_graph;
    }
    return dpins;
}

vector<UnGraph> GeneExpress::build_KPIN_min(const UnGraph* g) {
    int count = 12;
    map<string, double> active;
    for (auto& p : g->proteins) {
        // 查找基因表达
        double active_temp = 0.0;
        auto it = gene_express.find(p->protein_name);
        if (it != gene_express.end()) {
            auto mean_varience = calculate_mean_varience(it->second);
            if (essential_proteins.count(p->protein_name)) {
                // 是关键蛋白质，基因表达活性阈值更低一些，以保证能在生命周期中
                // 的大多数时间内保持活性
                active_temp = mean_varience.first;
            } else { // 非关键蛋白质，具有更高的基因表达活性表达阈值
                active_temp =
                    mean_varience.first - pow(mean_varience.second, 0.5);
            }
        } else {
            active_temp = 0.0;
        }
        active.insert(make_pair(p->protein_name, active_temp));
    }

    // 构建动态网络
    vector<UnGraph> dpins(count);
    vector<set<string>> set_proteins(count);
    vector<vector<string>> list_edges(count);
    vector<vector<double>> edge_weight(count);
    for (auto& protein : g->proteins_name) {
        // 是够含有当前蛋白质的基因表达数据？
        auto it = gene_express.find(protein);
        if (it == gene_express.end()) {
            for (int i = 0; i < count; ++i) {
                set_proteins[i].insert(protein);
            }
            continue;
        }
        for (int i = 0; i < count; ++i) {
            double exp_temp_min = MAXFLOAT;
            vector<double> arr{it->second[i], it->second[12 + i],
                               it->second[24 + i]};
            for (auto& a : arr) {
                if (a < exp_temp_min) {
                    exp_temp_min = a;
                }
            }

            if (exp_temp_min > active[protein]) {
                set_proteins[i].insert(protein);
            }
        }
    }

    // update edges
    for (auto& e : g->edges) {
        for (int i = 0; i < count; ++i) {
            if (set_proteins[i].count(e->node_a->protein_name) &&
                set_proteins[i].count(e->node_b->protein_name)) {
                list_edges[i].emplace_back(e->node_a->protein_name);
                list_edges[i].emplace_back(e->node_b->protein_name);
                edge_weight[i].emplace_back(e->balanced_weight);
            }
        }
    }

    for (int i = 0; i < count; ++i) {
        std::cout << set_proteins[i].size() << "\t" << list_edges[i].size() / 2
                  << endl;
        UnGraph temp_graph(move(set_proteins[i]), move(list_edges[i]),
                           move(edge_weight[i]));
        dpins[i] = temp_graph;
    }
    return dpins;
}

// 关键蛋白质又更高的阈值
vector<UnGraph> GeneExpress::build_KPINN(const UnGraph* g) {
    int count = 36;
    map<string, double> active;
    for (auto& p : g->proteins) {
        // 查找基因表达
        double active_temp = 0.0;
        auto it = gene_express.find(p->protein_name);
        if (it != gene_express.end()) {
            auto mean_varience = calculate_mean_varience(it->second);
            if (essential_proteins.count(p->protein_name)) {
                // 是关键蛋白质，基因表达活性阈值更低一些，以保证能在生命周期中
                // 的大多数时间内保持活性
                active_temp = mean_varience.first - pow(mean_varience.second, 0.5);

            } else { // 非关键蛋白质，具有更高的基因表达活性表达阈值
                active_temp = mean_varience.first;
            }
        } else {
            active_temp = 0.0;
        }
        active.insert(make_pair(p->protein_name, active_temp));
    }

    // 构建动态网络
    vector<UnGraph> dpins(count);
    vector<set<string>> set_proteins(count);
    vector<vector<string>> list_edges(count);
    vector<vector<double>> edge_weight(count);
    for (auto& protein : g->proteins_name) {
        // 是够含有当前蛋白质的基因表达数据？
        auto it = gene_express.find(protein);
        if (it == gene_express.end()) {
            for (int i = 0; i < count; ++i) {
                set_proteins[i].insert(protein);
            }
            continue;
        }
        for (int i = 0; i < count; ++i) {
            if (it->second[i] > active[protein]) {
                set_proteins[i].insert(protein);
            }
        }
    }

    // update edges
    for (auto& e : g->edges) {
        for (int i = 0; i < count; ++i) {
            if (set_proteins[i].count(e->node_a->protein_name) &&
                set_proteins[i].count(e->node_b->protein_name)) {
                list_edges[i].emplace_back(e->node_a->protein_name);
                list_edges[i].emplace_back(e->node_b->protein_name);
                edge_weight[i].emplace_back(e->balanced_weight);
            }
        }
    }

    for (int i = 0; i < count; ++i) {
        std::cout << set_proteins[i].size() << "\t" << list_edges[i].size() / 2
                  << endl;
        UnGraph temp_graph(move(set_proteins[i]), move(list_edges[i]),
                           move(edge_weight[i]));
        dpins[i] = temp_graph;
    }
    return dpins;
}
