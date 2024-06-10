/*
 * @Description: 
 * @Author: jh
 * @Date: 2024-04-30 17:11:48
 * @LastEditTime: 2024-06-01 22:09:48
 */

#ifndef COMPLEX_PREDICT_UNGRAPH_H
#define COMPLEX_PREDICT_UNGRAPH_H

#include "dag.h"

#include <iostream>
#include <set>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <memory>
#include <queue>

class Complex;
class Protein;
class Edge;

// 是否为加权图
typedef bool Weight;
typedef shared_ptr<Protein> ProteinPtr;
typedef shared_ptr<Edge> EdgePtr;


class Protein {
public:
    int id;
    string protein_name;
    set<ProteinPtr> neighbor;
    double weight;

public:
    Protein(int _id, string _protein_name, double _weight = 0.0);
    void display()const;
    void add_neighbor(const ProteinPtr& protein);

    void remove_neighbor(const ProteinPtr& protein);
    double jaccard_similarity(const ProteinPtr& other) const;
    double more_jaccard_similarity(const ProteinPtr& other) const;
    int degree() const;
    /**
     * @brief: 计算两个节点之间的单向相似性
     * @param {ProteinPtr&} other
     * @return pair{self -> other sim, other -> other sim}
     */    
    pair<double, double> unidirectional_similarity(const ProteinPtr& other) const;
    static bool ProteinCompareByWeight(const ProteinPtr& p1, const ProteinPtr& p2);
};


class Edge {
public:
    explicit Edge(const Edge *pEdge);

    ProteinPtr node_a;
    ProteinPtr node_b;
    double weight;
    double balanced_weight;

    Edge(const ProteinPtr& node_a, const ProteinPtr& node_b,
        double _weight = 0.0, double _balanced_weight = 0.0, int _visited_count = 0);
    bool operator<(const Edge&)const;
    // static bool CompareByBalancedWeight(const Edge* e1, const Edge* e2);
    // static bool CompareByBalancedWeight(const Edge& e1, const Edge& e2);
    static bool CompareByBalancedWeight(const EdgePtr& e1, const EdgePtr& e2);
    void display() const;
};

class UnGraph {
public:
    static DAG child_ancestor;
    static DAG ancestor_child;
    // static map<string, set<string>> protein_go;
    
    vector<ProteinPtr> ID2Protein;
    map<ProteinPtr, int> Protein2ID;
    map<string, int> protein_name_id;
    set<ProteinPtr> proteins;
    set<string> proteins_name;
    vector<EdgePtr> edges;
    vector<vector<bool>> connected;    // 存储个节点是否直接相连
    map<set<string>, int> Edge2ID;
    
    /**
     * @description: 从文件中读取文件
     * @param {string} ppi_file PPI连边文件
     * @param {weight} 是否为加权图（边加权）
     */    
    explicit UnGraph(string& ppi_file, Weight weight = false);
    /**
     * @description: 从节点和连边列表创建图
     * @param {&&} set_proteins: 蛋白质列表
     * @param {&&}  list_edges: PPI连边列表 每相邻两个蛋白质为一组边 [0, 1] [2, 3]，
     *              list_edges默认为空，即为无权重网络
     */    
    UnGraph(set<string>&& set_proteins, vector<string>&& list_edges, vector<double>&& edge_weight = vector<double>());
    UnGraph() =  default;
    ~UnGraph();
    int protein_count() const;
    int edge_count() const;
    void display() const;
    void sort_edges_balanced_weight();
    EdgePtr getEdge(const ProteinPtr& protein1, const ProteinPtr& protein2);
    double agglomeration_coefficient(const vector<ProteinPtr>& nodes);
    void weight_by_go_term();
    void clean_protein_weight();

    // ewca 相关算法
    __attribute__((unused)) void calculate_structure_similarty(vector<vector<double>>&) const;
    void get_common_neighbor_size(vector<vector<int>>&) const;
    void get_JCS(vector<vector<double>>&, const vector<vector<int>>&) const;
    void get_CNS(vector<vector<double>>&, const vector<vector<double>>&) const;

    /**
     * @brief: 通过jaccard计算两个节点之间的同质性指数
     * @return 返回同质性指数
     */    
    static double homogeneity_index_jaccard(const ProteinPtr&, const ProteinPtr&);

    // cns相关算法
    void calculate_attraction(vector<double>& attractions) const;
    void calculate_average_attraction(vector<double>& av_attractions) const;

    // 计算转移概率
    void calculate_walk_probability(vector<vector<double>>& probability) const;


    // BOPS相关算法
//    static int get_fa(int fa[], int x);
    void calculate_balanced_weight();
    static int find_parent(int protein, map<int, int>& parent);
    static void split_graph(queue<UnGraph>& ppi_queue, vector<UnGraph>& splited_ppi);
    static void get_complexes(UnGraph& origin_ppi, UnGraph& g, vector<set<string>>& complexes, double similarity_threshold);
    static void get_complexes1(UnGraph& g, vector<set<string>>& complexes, double similarity_threshold);
    static map<string, set<string>> read_protein_go(string file_path);
    /**
     * @brief: 为节点添加权重。节点的权重是连边的权重之和
     * @return 返回节点的权重。节点 i 的权重存储在对应的下标位置
     */
    vector<double> calculate_protein_weight() const;
    void write_to_file(const string& file_path) const;
    void write_xgmml(string& file_path) ;
    void normalize_edge_weight_min_max();

    static void split_graph1(queue<UnGraph>& ppi_queue, vector<UnGraph>& splited_ppi);
    static void get_complexes1(UnGraph& origin_ppi, UnGraph& g, vector<set<string>>& complexes, double similarity_threshold);
private:
    static void read_edge_list(string&, set<string>&, vector<string>&);
    static void read_edge_list_with_weight(string&, set<string>&, vector<string>&, vector<double>&);
    void add_edge(const ProteinPtr&, const ProteinPtr&, double weight = 0.0);
    void update_edge_id();

    static bool compare_pairs(const pair<EdgePtr, int>& pair1, const pair<EdgePtr, int>& pair);
};

class Complex {
public:
    double cohesion;
    vector<string> proteins;

    Complex(const UnGraph& g, vector<string>& _proteins);
    double complex_match_score(Complex& other);
    static bool CompareBySize(const Complex& a, const Complex& b);
    static bool CompareSetBySize(const set<string>& a, const set<string>& b);
    static set<set<string>> read_complex_from_file(string& file_path);
    static void write_complex_to_file(vector<Complex>&, string& file_path);
    static void write_complex_to_file(vector<set<string>>&, string& file_path);
    static void write_complex_to_file(set<set<string>>&, string& file_path);
    static double calculate_cohesion(const UnGraph& g, set<string>& _proteins);
    static void update_complexes(vector<Complex>& complexes, Complex& complex);
    static void update_complexes(vector<set<string>>& complexes, set<string>& complex, double threshold = 0.45);
    static bool evaluate_by_weight(UnGraph& g, set<string>& complex);
    static bool evaluate_by_density(UnGraph& g, set<string>& complex, double threshold = 0.50);
    void display() const;
};

#endif
