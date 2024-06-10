#include "config.h"
#include "gene_express.h"
#include "graph.h"

#include <cmath>
#include <sstream>

class KPIN_CNS {
  public:
    Graph* g;

    vector<double> node_weight;        // 节点的权重
    vector<vector<double>> attraction; // 节点之间的相互吸引力
    vector<double> node_attraction;    // 节点的影响力之和
    vector<vector<double>> uni_sim;

    KPIN_CNS(Graph* _g) : g(_g) {}
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
        for (auto& e : g->edges) {
            node_weight[e.smllar] += e.weight;
            node_weight[e.bigger] += e.weight;
        }
    }

    void calculate_attraction() {
        attraction.resize(g->node_count, vector<double>(g->node_count, 0));
        node_attraction.resize(g->node_count, 0.0);
        // 一介邻居
        for (auto& e : g->edges) {
            const double dis = 1.0 - log(1.0 / (1 + e.weight));
            const double f =
                node_weight[e.smllar] * node_weight[e.bigger] / pow(dis, 2);
            node_attraction[e.smllar] += f;
            node_attraction[e.bigger] += f;
            attraction[e.smllar][e.bigger] = f;
            attraction[e.bigger][e.smllar] = f;
        }
        // 找出所有二阶邻居"边"
        set<Edge> edges;
        for (int i = 0; i < g->node_count; ++i) {
            auto neighbor = g->get_neighbor(i); // 一阶邻居
            for (auto j : neighbor) {
                auto neighbor2 = g->get_neighbor(j); // 二阶邻居
                for (auto k : neighbor2) {
                    if (i == k || g->connected[i][k])
                        continue;
                    Edge e(i, k);
                    edges.insert(e);
                }
            }
        }
        vector<Edge> edges2(edges.begin(), edges.end());
        g->weighted_by_go_term(edges2);
        // 结构相似性
        for (auto& e : edges2) {
            e.weight += g->jaccard_similarity(e.smllar, e.bigger);
        }

        // 归一化
        double max_weight = 0.0;
        double min_weight = 999999;
        for (auto& e : edges2) {
            if (e.weight > max_weight) {
                max_weight = e.weight;
            }
            if (e.weight < min_weight) {
                min_weight = e.weight;
            }
        }
        for (auto& e : edges2) {
            e.weight = (e.weight - min_weight) / (max_weight - min_weight);
        }

        // 计算吸引力
        for (auto& e :edges2) {
            const double dis = 1.0 - log(1.0 / (1 + e.weight));
            const double f =
                node_weight[e.smllar] * node_weight[e.bigger] / pow(dis, 2);
            node_attraction[e.smllar] += f;
            node_attraction[e.bigger] += f;
            attraction[e.smllar][e.bigger] = f;
            attraction[e.bigger][e.smllar] = f;
        }
    }

    void update_uni_sim() {
        uni_sim.resize(g->node_count, vector<double>(g->node_count, 0.0));
        for (auto& e : g->edges) {
            pair<double, double> it =
                g->unidirectional_similarity(e.smllar, e.bigger);
            uni_sim[e.smllar][e.bigger] = it.first;
            uni_sim[e.bigger][e.smllar] = it.second;
        }
    }

    set<set<string>> _m_find_complex(double threshold = 0.4) {
        set<set<string>> complexes; // core？
        for (int i = 0; i < g->node_count; ++i) {
            set<int> temp_complex{i};
            set<int> neighbor = g->get_neighbor(i);
            for (auto& j : neighbor) {
                if (uni_sim[i][j] > 0.5) {
                    temp_complex.insert(j);
                }
            }

            // candidate
            set<int> candidate; // attrachment
            for (auto& n : temp_complex) {
                auto neighbor = g->get_neighbor(n);
                for (auto& nei : neighbor) {
                    if (!temp_complex.count(nei) && !candidate.count(nei)) {
                        // 计算和社区内的节点的吸引力之和
                        double sum = 0.0;
                        for (auto& c : temp_complex) {
                            sum += attraction[c][nei];
                        }
                        if (sum / node_attraction[nei] >= threshold) {
                            candidate.insert(nei);
                        }
                    }
                }
            }
            temp_complex.insert(candidate.begin(), candidate.end());
            set<string> complex;
            for (auto& n : temp_complex) {
                complex.insert(g->id_protein[n]);
            }
            if(complex.size() <= 2 || complexes.size() >= 20) continue;
            complexes.insert(complex);
        }
        return complexes;
    }
};

DAG Graph::ancestor_child(IS_A_A2C, PART_OF_A2C);
DAG Graph::child_ancestor(IS_A_C2A, PART_OF_C2A);
int main(int argc, char const** argv) {
    string file_path;
    string result;
    fstream file("config.txt");  // 配置文件
    string line;
    while(getline(file, line)) {
        istringstream iss(line);
        iss >> file_path;
        iss >> result;
    }

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
        KPIN_CNS kpin(&pg);
        auto temp_complexes = kpin.find_complex();
        complexes.insert(temp_complexes.begin(), temp_complexes.end());
    }
    vector<set<string>> temp_complexes(complexes.begin(), complexes.end());
    sort(temp_complexes.begin(), temp_complexes.end(), Complex::CompareSetBySize<string>);
    cout << temp_complexes.size() << endl;
    // vector<set<string>> result_complexes;
    // for(auto& c: temp_complexes) { 
    //     Complex::update_complexes(result_complexes, c, 0.65);
    // }

    Complex::write_complex(result, temp_complexes);
    // cout << result_complexes.size() << endl;

    return 0;
}