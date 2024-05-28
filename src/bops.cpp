/*
 * @brief:
 * @Author: jh
 * @Date: 2024-05-06 16:00:21
 * @LastEditTime: 2024-05-28 17:21:32
 */
#include "../include/config.h"
#include "../include/dag.h"
#include "../include/ungraph.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>

using namespace std;

DAG UnGraph::ancestor_child(IS_A_A2C, PART_OF_A2C);
DAG UnGraph::child_ancestor(IS_A_C2A, PART_OF_C2A);
map<string, set<string>> UnGraph::protein_go = UnGraph::read_protein_go(GO_SLIM);

// 只加权一次
int main(int argc, char** argv) {
    int opt;
    string ppi_file;
    string out_file;

    while ((opt = getopt(argc, argv, "i:o:")) != -1) {
        string temp = optarg;
        switch (opt) {
        case 'i':
            ppi_file = "/home/jh/code/JHPC/dataset/Yeast/PPI/" + temp + ".txt";
            break;
        case 'o':
            out_file = "/home/jh/code/JHPC/result" + temp + ".txt";
            break;
        default:
            std::cerr << "Usage: " << argv[0] << " -n <value> -s <value>"
                      << std::endl;
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
    std::cout << "proteins: " << g.proteins.size()
              << "\t edges: " << g.edges.size() << endl;

    std::cout << "Weightting ppi ... \t";
    // 使用go-DAG加权
    g.weight_by_go_term();
    // 同质性指数加权
    for (auto& e : g.edges) {
        e->balanced_weight += e->node_a->jaccard_similarity(e->node_b);
    }
    // 对便权重计算平衡系数
    // g.calculate_balanced_weight();
    std::cout << "weight ppi finished!" << endl;

    std::cout << "Splitting ppi ... \t";
    queue<UnGraph> queue_ppi;
    queue_ppi.push(std::move(g));
    vector<UnGraph> splitted_ppi;
    while (!queue_ppi.empty()) {
        UnGraph::split_graph(queue_ppi, splitted_ppi);
    }
    std::cout << "split ppi finished! \t split ppi count: "
              << splitted_ppi.size() << endl;

    std::cout << "detecting complex ... \t";
    vector<set<string>> complexes;
    for (auto& sg : splitted_ppi) {
        vector<set<string>> temp_complexes;
        UnGraph::get_complexes(g, sg, temp_complexes, 0.45);
        complexes.insert(complexes.end(), temp_complexes.begin(),
                         temp_complexes.end());
        std::cout << complexes.size() << endl;
    }
    sort(complexes.begin(), complexes.end(), Complex::CompareSetBySize);
    std::cout << "detected complex finished!";

    Complex::write_complex_to_file(complexes, out_file);
    std::cout << "Write result to file: " << out_file << endl;
    return 0;
}
