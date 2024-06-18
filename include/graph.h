#ifndef __GRAPH_H__
#define __GRAPH_H__

#include "dag.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>

using namespace std;

class Edge {
  public:
    int smllar;
    int bigger;
    double weight;

    Edge(int a = 0, int b = 0, double __weight = 0.0);
    bool operator==(const Edge&) const;
    bool operator<(const Edge&) const;
    static bool CompareByWeight(const Edge& e1, const Edge& e2);
    void debug() const;
};

// 简单无向图
class Graph {
  private:
    struct PairHash {
        size_t operator()(const std::pair<int, int>& p) const {
            return std::hash<int>()(p.first) ^ std::hash<int>()(p.second);
        }
    };
    void read_edge(string& file_path, vector<int>& edge_list,
                   map<string, int>& __protein_id,
                   map<int, string>& __id_protein) const;

  public:
    static DAG child_ancestor;
    static DAG ancestor_child;

    int node_count;
    int edge_count;
    vector<vector<bool>> connected;
    map<int, set<int>> adjancy_list; // 邻接表
    vector<Edge> edges;
    std::unordered_map<std::pair<int, int>, int, PairHash> edge_id;
    map<string, int> protein_id;
    map<int, string> id_protein;

    Graph() = default;
    explicit Graph(string& ppi_file);
    explicit Graph(vector<string>& edge_list, vector<double>* = nullptr);
    ~Graph() = default;
    set<int> get_neighbor(int n) const;
    void add_edge(int u, int v, double weight = 0.0);
    // 谨慎使用不安全  
    void remove_edge(int u, int v);
    int get_edge_id(int u, int v) const;
    int degree(int id) const;

    // set<int> get_common_neighbor(int u, int v) const;
    double jaccard_similarity(int u, int v) const;
    double jaccard_similarity_more(int u, int v) const;
    void weighted_by_go_term();
    void weighted_by_go_term(vector<Edge>& edges) const;
    // void weighted_by_go_term_Lin();
    void normalize_edge_weight_min_max();
    void update_edge_id();
    void calculate_balanced_weight(double balanced_index = 1.5);
    set<int> get_common_neighbor(int u, int v) const;
    pair<double, double> unidirectional_similarity(int u, int v) const;
    double density(set<int>& nodes) const;
    double density(set<int>& nodes1, set<int>& nodes2) const;
    double evaluate(set<int>& complex) const;

};

namespace Complex {

template<typename T>
bool CompareSetBySize(const set<T>& a, const set<T>& b) {
    return a.size() > b.size();
}

template <typename T> double overlapping_score(set<T>& a, set<T>& b) {
    set<T> common;
    set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                     inserter(common, common.end()));
    return static_cast<double>(common.size()) /
           static_cast<double>(a.size(), b.size());
}
/**
 * @brief Updates the complexes vector by adding a new complex if it does not
 * have a significant overlap with existing complexes.
 *
 * This function iterates over the existing complexes in the complexes vector
 * and checks if the overlap score (OS) between the new complex and any existing
 * complex is greater than or equal to the given threshold. If no significant
 * overlap is found, the new complex is added to the complexes vector.
 *
 * @param complexes A vector of sets, where each set represents a complex.
 * @param complex The new complex to be added to the complexes vector.
 * @param threshold The minimum overlap score required for a new complex to be
 * considered a significant overlap with an existing complex.
 *
 * @return void
 */
void update_complexes(vector<set<string>>& complexes, set<string>& complex,
                      double threshold = 0.90);
/**
 * @brief Writes the complexes to a file.
 *
 * This function writes the complexes to a file, where each complex is
 * represented as a line with tab-separated protein identifiers.
 *
 * @param file_path The path to the output file.
 * @param id_protein A map that maps protein IDs to their identifiers.
 * @param complexes A vector of sets, where each set represents a complex.
 *
 * @return void
 *
 * @throws std::runtime_error If the file cannot be opened.
 */
void write_complex(const string& file_path, const map<int, string>& id_protein,
                   const vector<set<int>>& comeplxes);

void write_complex(const string& file_path, const set<set<string>>& comeplxes);
void write_complex(const string& file_path,
                   const vector<set<string>>& comeplxes);
}; // namespace Complex

#endif
