#include "../include/ungraph.h"

#include <iostream>

auto get_seed_nodes(UnGraph& g) {
    set<ProteinPtr> seed_nodes;
    // 计算每个节点的影响力
    double max = 0, min = 999999;
    double d = 0;
    for(auto& e: g.edges) {
        // d =  1 - log10(0.1 / (1 + e->balanced_weight));
        d =  -log10(e->balanced_weight);
        std::cout << e->balanced_weight << "\t" << d << std::endl;
        if(d > max) max = d;
        if(d < min) min = d;
    }
    std::cout << "max: " << max << std::endl;
    std::cout << "min: " << min << std::endl;
}


int main() {
    string ppi_file = COLLINS_PPI;
    string result_file = "./cns0.txt";

    UnGraph g(ppi_file, false);
    DAG dag;
    BioInformation bio;
    g.weight_by_go_term(bio, dag);
    get_seed_nodes(g);

    return 0;
}
