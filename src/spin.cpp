/*
 * @brief: 
 * @Author: jh
 * @Date: 2024-05-19 00:02:34
 * @LastEditTime: 2024-05-28 17:10:38
 */
#include "../include/ungraph.h"
#include "../include/subcelluar.h"
#include "../include/bio_information.h"
#include "../include/dag.h"

DAG UnGraph::ancestor_child(IS_A_A2C, PART_OF_A2C);
DAG UnGraph::child_ancestor(IS_A_C2A, PART_OF_C2A);
map<string, set<string>> UnGraph::protein_go = UnGraph::read_protein_go(GO_SLIM);

int main() {
    string result_file = "./spin_collins.txt";
    string subcelluar_file = "/home/jh/code/JHPC/dataset/Yeast/DAG/subcellular.txt";
    Subcellular s(subcelluar_file);
    cout << s.location_protein.size() << endl;
    
    string ppi_file = COLLINS_PPI;
    UnGraph g(ppi_file);


    g.weight_by_go_term();
    g.calculate_balanced_weight();
    
    auto spin = s.build_spin(g);
    std::cout << spin.size() << endl;

    int count = 0;
    for(auto& sg: spin) {
        if(sg.proteins.size() >= 3 && sg.proteins.size() - 2 <= sg.edges.size()) {
            count += 1;
        }
    }
    std::cout << count << endl;
    // build queue<UnGraph>
    queue<UnGraph> q;
    for(auto& sg: spin) {
        q.push(sg);
    }

    vector<UnGraph> splitted_ppi;
    while(!q.empty()) {
        UnGraph::split_graph(q, splitted_ppi);
    }
    std::cout << "split ppi finish!" << endl;
    vector<set<string>> complex;
    for(auto& sg: splitted_ppi) {
        vector<set<string>> temp_complexes;
        UnGraph::get_complexes(g, sg, temp_complexes, 0.45);
        complex.insert(complex.end(), temp_complexes.begin(), temp_complexes.end());
    }

    sort(complex.begin(), complex.end(), Complex::CompareSetBySize);
    Complex::write_complex_to_file(complex, result_file);
    return 0;
}
