/*
 * @brief:
 * @Author: jh
 * @Date: 2024-05-06 13:23:07
 * @LastEditTime: 2024-05-28 17:10:47
 */
#include "../include/bio_information.h"
#include "../include/dag.h"
#include "../include/gene_express.h"
#include "../include/ungraph.h"

#include <algorithm>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <vector>

DAG UnGraph::ancestor_child(IS_A_A2C, PART_OF_A2C);
DAG UnGraph::child_ancestor(IS_A_C2A, PART_OF_C2A);
map<string, set<string>> UnGraph::protein_go = UnGraph::read_protein_go(GO_SLIM);

class CNS {
  public:
    UnGraph* graph;
    map<ProteinPtr, double> protein_weight;
    vector<vector<double>> attraction; // 记录节点之间的相互吸引力
    map<ProteinPtr, double> protein_attraction; // 节点的影响力之和
    explicit CNS(UnGraph* g) { graph = g; }
    /**
     * @brief: 计算每个节点的权重，节点的权重为关联边权重之和
     */
    auto calculateod_node_weight() {
        for (const auto& e : graph->edges) {
            protein_weight[e->node_a] += e->balanced_weight;
            protein_weight[e->node_b] += e->balanced_weight;
        }
    }

    /**
     * @brief: 计算每两个节点直接的吸引力，以及节点的影响力
     */
    auto calculate_attraction() {
        attraction.resize(graph->proteins.size(),
                          vector<double>(graph->proteins.size(), 0.0));
        calculateod_node_weight();
        for (const auto& e : graph->edges) {
            const double d = 1.0 - log(1.2 - e->balanced_weight);
            const double score = protein_weight[e->node_a] *
                                 protein_weight[e->node_b] / pow(d, 2);
            // e->display();
            // std::cout << d << "\t attraction: " << score << endl;
            protein_attraction[e->node_a] += score;
            protein_attraction[e->node_b] += score;
            attraction[e->node_a->id][e->node_b->id] = score;
        }
    }

    auto find_compelx() -> set<set<string>> {
        calculate_attraction();
        set<set<string>> complexes;
        for (auto& p : graph->proteins) {
            set<string> temp_complex{p->protein_name};
            for (auto& nei : p->neighbor) {
                temp_complex.insert(nei->protein_name);
            }

            // expand
            set<string> candidate;
            // temp_complex 成员
            for (auto& m : temp_complex) {
                // 邻居， 如果已经包含在temp_complex 或者candidate 中跳过
                for (auto& nei :
                     graph->ID2Protein[graph->protein_name_id[m]]->neighbor) {
                    if (candidate.count(nei->protein_name) ||
                        temp_complex.count(nei->protein_name)) {
                        continue;
                    }
                    double sum_attraction = 0.0;
                    for (auto& nnei : nei->neighbor) {
                        if (temp_complex.count(nnei->protein_name)) {
                            sum_attraction += attraction[nnei->id][nei->id];
                        }
                    }
                    if (sum_attraction / protein_attraction[nei] >= 0.40) {
                        candidate.insert(nei->protein_name);
                    }
                }
            }
            temp_complex.insert(candidate.begin(), candidate.end());
            if (temp_complex.size() <= 2 || temp_complex.size() >= 30)
                continue;
            if (!Complex::evaluate_by_weight(*graph, temp_complex))
                continue;
            complexes.insert(temp_complex);
        }
        return complexes;
    }
};

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

    vector<set<string>> complexes;

    BioInformation bio;
    UnGraph g(ppi_file);
    GeneExpress exp;
    vector<UnGraph> dpins = exp.build_KPIN(&g);
    for (auto& dg : dpins) {
        // 加权
        dg.weight_by_go_term();
        // 同质性指数
        for (auto& e : dg.edges) {
            const double s = e->node_a->jaccard_similarity(e->node_b);
            e->weight += s;
            e->balanced_weight += s;
        }
        // 平衡
        dg.calculate_balanced_weight();
        CNS cns(&dg);
        auto temp_complexes = cns.find_compelx();
        complexes.insert(complexes.end(), temp_complexes.begin(),
                         temp_complexes.end());
    }
    std::sort(complexes.begin(), complexes.end(), Complex::CompareSetBySize);
    vector<set<string>> result;
    for (auto& c : complexes) {
        // if(Complex::evaluate_by_weight(UnGraph &g, set<string> &complex))
        Complex::update_complexes(result, c);
    }

    Complex::write_complex_to_file(result, out_file);
    return 0;
}
