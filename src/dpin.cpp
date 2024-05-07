/*
 * @brief: 
 * @Author: jh
 * @Date: 2024-05-06 16:00:21
 * @LastEditTime: 2024-05-07 16:03:48
 */
#include "bio_information.h"
#include "config.h"
#include "dag.h"
#include "gene_express.h"
#include "ungraph.h"
#include <algorithm>
#include <cstdio>
#include <queue>
#include <set>
#include <vector>

void display(string& ppi_file, string& result_file) {
    printf("ppi_file: %s\n", ppi_file.c_str());
    printf("result_file: %s\n", ppi_file.c_str());
}

int main() {
    string ppi_file = COLLINS_PPI;
    string result_file = "./dpin_result.txt";
    display(ppi_file, result_file);
    BioInformation bio;
    DAG dag;
    UnGraph g(ppi_file, false);
    g.weight_by_go_term(bio, dag);
    g.calculate_balanced_weight();

    GeneExpress gene_express(GENE_EXPRESSION);
    vector<UnGraph> dpins = gene_express.build_dynamic_PPI(&g);
    std::cout << "dpins.size " << dpins.size() << endl;

    queue<UnGraph> queue_ppi;
    for(auto& dg: dpins) {
        queue_ppi.push(dg);
    }

    vector<UnGraph> splitted_ppi;
    while(!queue_ppi.empty()) {
        UnGraph::split_graph(queue_ppi, splitted_ppi, bio, dag);
    }
    std::cout << splitted_ppi.size() << std::endl;

    vector<set<string>> complexes;
    for(auto& sg: splitted_ppi) {
        vector<set<string>> temp_complexes;
        UnGraph::get_complexes(sg, temp_complexes, 0.45);
        complexes.insert(complexes.end(), temp_complexes.begin(), temp_complexes.end());
        std::cout << complexes.size() << endl;
    }
    sort(complexes.begin(), complexes.end(), Complex::CompareSetBySize);
    
    // cout << "last" << endl;
    // vector<set<string>> result;
    // for(auto& c: complexes) {
    //     Complex::update_complexes(result, c);
    // }

    Complex::write_complex_to_file(complexes, result_file);

    return 0;
}