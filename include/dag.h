/*
 * @Description: 
 * @Author: jh
 * @Date: 2024-04-30 17:11:48
 * @LastEditTime: 2024-05-07 19:38:53
 */
#ifndef __DAG_H__
#define __DAG_H__

#include "config.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <map>
#include <queue>
#include <vector>
#include <stack>
#include <utility>
#include <algorithm>

using namespace std;

enum Relation {
    IS_A,
    PART_OF,
    NONE,
};

// 每一个几点指向自己的祖先
class Node {
public:
    int id;
    string go;
    set<Node*> ancestor;

    explicit Node(int id) : id(id) {}
    Node(int id, string go): id(id), go(move(go)) {}
    void add(Node*);
    void print()const;
};

class DAG {
public:
    map<pair<string, string>, double> Similarity;
    vector<Node*> nodes;
    vector<vector<Relation>> relation;
    map<string, int> GO2ID;
public:
    // 需要传入两个参数，第一个是用来初始化所有的GOterm
    // 并且固定所有的GOTerm对应的参数 
    DAG();
    ~DAG();
    Node* addNode(int id, string&);
    void addEdge(int, int, Relation);
    void add_edges(vector<vector<string>>&, Relation);
    void print() const;
    double similarity(const string&, const string&);
//    set<Node*> get_all_ancestor(Node*);
    static void get_all_ancestors(Node*, set<Node*>&);
    static void read_ancestor_child(const string&, vector<vector<string>>&);
    static void read_all_go_terms(const string&, set<string>&);
    void calculate_SValue(Node*, map<Node*, double>&);
    // 查找a到b的全部路径
    void find_all_paths(Node* a, Node* b, vector<Node*>& path, vector<vector<Node*>>& all_path);
    static void print_path(const vector<vector<Node*>>&);
    double get_similarity_go_gos_by_max(const string&, const set<string>&);
    double get_similarity_protein(const set<string>& gos1, const set<string>& gos2);
private:
    double similarity(int, int);
};

#endif // __DAG_H__
