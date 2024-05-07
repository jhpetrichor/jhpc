/*
 * @brief: 
 * @Author: jh
 * @Date: 2024-05-06 16:00:21
 * @LastEditTime: 2024-05-07 19:48:04
 */
#include "bio_information.h"
#include "config.h"
#include "dag.h"
#include "ungraph.h"
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <queue>
#include <set>
#include <sstream>
#include <vector>

void display(string& ppi_file, string& result_file) {
    printf("ppi_file: %s\n", ppi_file.c_str());
    printf("result_file: %s\n", ppi_file.c_str());
}



// double cacluate_p_value(UnGraph& g, BioInformation& bio, set<string>& complex) {
//     string file_path = "/home/jh/code/JHPC/dataset/Yeast/DAG/protein-go.txt";
//     read_go_protein(file_path);
    
// }

int main() {
    string ppi_file = COLLINS_PPI;
    string result_file = "./bops_result.txt";
    string go_protein_file = "/home/jh/code/JHPC/dataset/Yeast/DAG/protein-go.txt";
    display(ppi_file, result_file);

    UnGraph g(ppi_file, false);
    BioInformation bio;
    DAG dag;
    g.weight_by_go_term(bio, dag);
    g.calculate_balanced_weight();

    queue<UnGraph> queue_ppi;
    queue_ppi.push(std::move(g));
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

    Complex::write_complex_to_file(complexes, result_file);

    return 0;
}