/*
 * @Description: 
 * @Author: jh
 * @Date: 2024-04-30 17:11:48
 * @LastEditTime: 2024-05-07 15:44:10
 */

#include "bio_information.h"
#include "dag.h"
#include "graph.h"
#include "ungraph.h"
#include <algorithm>
#include <iterator>
#include <utility>
#include <vector>

void display_massage() {

}

int main() {
    string ppi_file = COLLINS_PPI;

    UnGraph g(ppi_file, false);
    DAG dag;
    BioInformation bio;
    g.weight_by_go_term(bio, dag);
    // g.write_to_file("/home/jh/code/complex_predict/bin/collins.txt");
    g.calculate_balanced_weight();
    for(auto& e: g.edges) {
        e->balanced_weight += e->node_a->jaccard_similarity(e->node_b);
    }
    std::queue<UnGraph> queue_ppi;
    queue_ppi.push(std::move(g));
    std::vector<UnGraph> splitted_ppi;
    while (!queue_ppi.empty()) {
        UnGraph::split_graph(queue_ppi, splitted_ppi, bio, dag);
    }
    std::cout << "splitted_ppi size: " << splitted_ppi.size() << std::endl;
    vector<std::set<std::string>> complexes;
    for (auto &g: splitted_ppi) {
        g.get_complexes1(g, complexes, 0.45);
        std::cout << complexes.size() << std::endl;
        std::cout << g.proteins.size() << "\t " << g.edges.size() << endl;
    }
    std::cout << "End get_complex! " << complexes.size() << endl;
    //
    vector<set<string>> result;
    for(auto& c: complexes) {
        bool matched = false;
        for(auto& cc: result) {

            vector<string> common;
            set_intersection(c.begin(), c.end(),
                            cc.begin(), cc.end(),
                            std::inserter(common, common.begin()));
            double score = (double)common.size() / (double)std::max(c.size(), cc.size());
            if(score >= 0.65) {
                matched = true;
                break;
            }
        }
        if(!matched) {
            result.emplace_back(c);
        }
    }


    std::string rest_file = "./test0.txt";
    Complex::write_complex_to_file(result, rest_file);

    return 0;
}
