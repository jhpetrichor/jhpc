/*
 * @brief: 
 * @Author: jh
 * @Date: 2024-05-06 16:00:21
 * @LastEditTime: 2024-05-09 16:10:30
 */
/*
 * @brief: 
 * @Author: jh
 * @Date: 2024-05-06 16:00:21
 * @LastEditTime: 2024-05-09 15:17:55
 */
#include "bio_information.h"
#include "config.h"
#include "dag.h"
#include "ungraph.h"
#include <algorithm>
#include <bits/getopt_core.h>
#include <cstdio>
#include <fstream>
#include <queue>
#include <set>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <unistd.h>

// 只加权一次
int main(int argc, char** argv) {
    int opt;
    string ppi_file;
    string out_file;

    while ((opt = getopt(argc, argv, "i:o:")) != -1) {
        string temp = optarg;
        switch (opt) {
            case 'i':
                ppi_file = "/home/jh/code/JHPC/dataset/Yeast/PPI/" + temp +".txt";
                break;
            case 'o':
                out_file = "/home/jh/code/JHPC/result" + temp + ".txt";
                break;
            default:
                std::cerr << "Usage: " << argv[0] << " -n <value> -s <value>" << std::endl;
                return 1;
        }
    }

    if (ppi_file.empty() || out_file.empty()) {
        std::cerr << "Error: -i and -n options require values." << std::endl;
        return 1;
    }
    cout << "ppi_file: " << ppi_file << endl;
    cout << "out_file: " << out_file << endl;

    UnGraph g(ppi_file);
    std::cout << "Read ppi finished!\t";
    std::cout << "proteins: " <<  g.proteins.size() << "\t edges: " << g.edges.size() << endl;
    BioInformation bio;
    DAG dag;
    
    std::cout << "Weightting ppi ... \t";
    // 使用go-DAG加权
    g.weight_by_go_term(bio, dag);
    // 同质性指数加权
    for(auto& e: g.edges) {
        e->balanced_weight += e->node_a->jaccard_similarity(e->node_b);
    }
    // 对便权重计算平衡技术
    g.calculate_balanced_weight();
    std::cout << "weight ppi finished!" << endl;

    std::cout << "Splitting ppi ... \t";
    queue<UnGraph> queue_ppi;
    queue_ppi.push(std::move(g));
    vector<UnGraph> splitted_ppi;
    while(!queue_ppi.empty()) {
        UnGraph::split_graph(queue_ppi, splitted_ppi, bio, dag);
    }
    std::cout << "split ppi finished! \t split ppi count: " << splitted_ppi.size() << endl;

    std::cout << "detecting complex ... \t";
    vector<set<string>> complexes;
    for(auto& sg: splitted_ppi) {
        vector<set<string>> temp_complexes;
        UnGraph::get_complexes(g,  sg, temp_complexes, 0.45);
        complexes.insert(complexes.end(), temp_complexes.begin(), temp_complexes.end());
        std::cout << complexes.size() << endl;
    }
    sort(complexes.begin(), complexes.end(), Complex::CompareSetBySize);
    std::cout << "detected complex finished!";

    Complex::write_complex_to_file(complexes, out_file);
    std::cout << "Write result to file: " << out_file << endl;
    return 0;
}