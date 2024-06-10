#include "graph.h"
#include "config.h"
#include "gene_express.h"

#include <cmath>
#include <sstream>

class KPIN_CNS {
  public:
    Graph* g;
    double ct;
    double at;
    
    vector<double> node_weight;        // 节点的权重
    vector<vector<double>> attraction; // 节点之间的相互吸引力
    vector<double> node_attraction;    // 节点的影响力之和
    vector<vector<double>> uni_sim;

    KPIN_CNS(Graph* _g, double _ct, double _at) {
        g = _g;
        ct = _ct;
        at = _at;
    }
    set<set<string>> find_complex() {
        calculate_node_weight();
        calculate_attraction();
        update_uni_sim();
        set<set<string>> complexes = _m_find_complex();
        return complexes;
    }

  private:
    void calculate_node_weight() {
        node_weight.resize(g->node_count, 0.0);
        for (auto& e: g->edges) {
            node_weight[e.smllar] += e.weight;
            node_weight[e.bigger] += e.weight;
        }
    }

    void calculate_attraction() {
        attraction.resize(g->node_count, vector<double>(g->node_count, 0));
        node_attraction.resize(g->node_count, 0.0);
        for (auto& e: g->edges) {
            const double dis = 1.0 - log(1.0 / (1 + e.weight));
            const double f = node_weight[e.smllar] * node_weight[e.bigger] / pow(dis, 2);
            node_attraction[e.smllar] += f;
            node_attraction[e.bigger] += f;
            attraction[e.smllar][e.bigger] = f;
            attraction[e.bigger][e.smllar] = f;
        }
    }

    void update_uni_sim() {
        uni_sim.resize(g->node_count, vector<double>(g->node_count, 0.0));
        for(auto& e: g->edges) {
            pair<double, double> it = g->unidirectional_similarity(e.smllar, e.bigger);
            uni_sim[e.smllar][e.bigger] = it.first;
            uni_sim[e.bigger][e.smllar] = it.second;
        }
    }

    set<set<string>> _m_find_complex() {
        set<set<string>> complexes;    // core？
        for(int i = 0; i < g->node_count; ++i) {
            set<int> temp_complex{i};
            set<int> neighbor = g->get_neighbor(i);
            for(auto& j: neighbor) {
                if(uni_sim[i][j] > ct) {
                    temp_complex.insert(j);
                }
            }

            // candidate
            set<int> candidate;  // attrachment
            for(auto& n: temp_complex) {
                auto neighbor = g->get_neighbor(n);
                for(auto& nei: neighbor) {
                    if(!temp_complex.count(nei) && !candidate.count(nei)) {
                        // 计算和社区内的节点的吸引力之和
                        double sum = 0.0;
                        for(auto& c: temp_complex) {
                            sum += attraction[c][nei];
                        }
                        if(sum / node_attraction[nei]  >= at) {
                            candidate.insert(nei);
                        }
                    }
                }
            }
            temp_complex.insert(candidate.begin(), candidate.end());
            set<string> complex;
            for(auto& n: temp_complex) {
                complex.insert(g->id_protein[n]);
            }
            if(complex.size() <= 2 || complex.size() >= 20) continue;
            complexes.insert(complex);
        }
        return complexes;
    }
};

DAG Graph::ancestor_child(IS_A_A2C, PART_OF_A2C);
DAG Graph::child_ancestor(IS_A_C2A, PART_OF_C2A);
int main(int argc, char const **argv) {
    string file_path = argv[1];
    string result = argv[2];
    double ct = atof(argv[3]);   // core threshold
    double at = atof(argv[4]);   // uni similarity threshold
    double ot = atof(argv[5]);   // overlap threshold

    // fstream file("config.txt");  // 配置文件
    // string line;
    // while(getline(file, line)) {
    //     istringstream iss(line);
    //     iss >> file_path;
    //     iss >> result;
    // }

    set<set<string>>  complexes;
    
    Graph g(file_path);
    GeneExpress ges;
    vector<Graph> graphs = ges.build_KPIN(&g);
    for(auto& pg: graphs) {
        pg.weighted_by_go_term();
        for(auto& e: pg.edges) {
            e.weight += pg.jaccard_similarity_more(e.smllar, e.bigger);
        }
        // 归一化
        pg.calculate_balanced_weight();
        pg.normalize_edge_weight_min_max();
        KPIN_CNS kpin(&pg, ct, at);
        auto temp_complexes = kpin.find_complex();
        complexes.insert(temp_complexes.begin(), temp_complexes.end());
    }
    vector<set<string>> temp_complexes(complexes.begin(), complexes.end());
    sort(temp_complexes.begin(), temp_complexes.end(), Complex::CompareSetBySize<string>);
    std::cout << complexes.size() << endl;
    vector<set<string>> result_complexes;
    for(auto& c: temp_complexes) { 
        Complex::update_complexes(result_complexes, c, ot);
    }

    Complex::write_complex(result, temp_complexes);
    // std::cout << result_complexes.size() << endl;
    return 0;
}