/*
 * @Description: 
 * @Author: jh
 * @Date: 2024-04-30 17:11:48
 * @LastEditTime: 2024-05-09 15:33:39
 */

#ifndef COMPLEX_PREDICT_UNGRAPH_H
#define COMPLEX_PREDICT_UNGRAPH_H

#include "bio_information.h"
#include "dag.h"

#include <set>
#include <map>
#include <string>
#include <vector>
#include <memory>

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

    void add_neighbor(const ProteinPtr& protein);

    void remove_neighbor(const ProteinPtr& protein);
    double jaccard_similarity(const ProteinPtr& other) const;
    int degree() const;

    static bool ProteinCompareByWeight(const ProteinPtr& p1, const ProteinPtr& p2);
};


class Edge {
public:
    explicit Edge(const Edge *pEdge);

    ProteinPtr node_a;
    ProteinPtr node_b;
    double weight;
    double balanced_weight;
    double visited_count;   // 访问次数

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
    void display() const;
    /**
     * @brief: 计算蛋白质复合物的p-value， 一般而言p-value越小，则越具有生物学意义，通常取value <= 0.01
     * @param {UnGraph&} g
     * @return p-value
     */    
    static double calculate_p_value(UnGraph& g, set<string>& complex);
    EdgePtr getEdge(const ProteinPtr& protein1, const ProteinPtr& protein2);
    double agglomeration_coefficient(const vector<ProteinPtr>& nodes);
    void weight_by_go_term(BioInformation& bio, DAG& dag);

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
    static void split_graph(queue<UnGraph>& ppi_queue, vector<UnGraph>& splited_ppi, BioInformation& bio, DAG& dag);
    static void get_complexes(UnGraph& origin_ppi, UnGraph& g, vector<set<string>>& complexes, double similarity_threshold);
    static void get_complexes1(UnGraph& g, vector<set<string>>& complexes, double similarity_threshold);

    /**
     * @brief: 为节点添加权重。节点的权重是连边的权重之和
     * @return 返回节点的权重。节点 i 的权重存储在对应的下标位置
     */
    vector<double> calculate_protein_weight() const;
    void write_to_file(const string& file_path) const;


    static void split_graph1(queue<UnGraph>& ppi_queue, vector<UnGraph>& splited_ppi, BioInformation& bio, DAG& dag);
    static void get_complexes1(UnGraph& origin_ppi, UnGraph& g, vector<set<string>>& complexes, double similarity_threshold);
private:
    static void read_edge_list(string&, set<string>&, vector<string>&);
    static void read_edge_list_with_weight(string&, set<string>&, vector<string>&, vector<double>&);
    void add_edge(const ProteinPtr&, const ProteinPtr&, double weight = 0.0);

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
    static void update_complexes(vector<set<string>>& complexes, set<string>& complex);
    static bool evaluate_by_weight(UnGraph& g, set<string>& complex);
    void display() const;
};

#endif
