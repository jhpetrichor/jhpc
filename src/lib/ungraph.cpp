/*
 * @brief: 
 * @Author: jh
 * @Date: 2024-05-06 13:23:07
 * @LastEditTime: 2024-05-09 15:50:00
 */
#include "ungraph.h"
#include "config.h"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <memory>
#include <ostream>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

Protein::Protein(int idx, string _protein_name, double _weight) {
    id = idx;
    protein_name = move(_protein_name);
    weight = _weight;
}

void Protein::add_neighbor(const ProteinPtr& _protein) {
    if(!_protein) {
        cerr << "Protein not founded!" << endl;
        return;
    }

    if(_protein->protein_name == this->protein_name) return;
    neighbor.insert(_protein);
}

void Protein::remove_neighbor(const ProteinPtr& _protein) {
    if(!_protein) {
        cerr << "Protein not founded!" << endl;
        return;
    }
    neighbor.erase(_protein);
}

double Protein::jaccard_similarity(const ProteinPtr& other) const{
    set<ProteinPtr> common;
    set_intersection(neighbor.begin(), neighbor.end(),
        other->neighbor.begin(), other->neighbor.end(),
        inserter(common, common.begin()));
    return static_cast<double>(common.size()) / static_cast<double>(degree() + other->degree() - common.size());
}

int Protein::degree() const {
    return static_cast<int>(neighbor.size());
}

bool Protein::ProteinCompareByWeight(const ProteinPtr& p1, const ProteinPtr& p2){
    return p1->weight < p2->weight;
}

Edge::Edge(const ProteinPtr& _node_a, const ProteinPtr& _node_b, double _weight,
    double _balanced_weight, int _visited_count) {
    if(!_node_a || !_node_b) {
        cerr << "Node not found!" << endl;
        return;
    }
    node_a = _node_a;
    node_b = _node_b;
    weight = _weight;
    balanced_weight = _balanced_weight;
    visited_count = _visited_count;
}

bool Edge::operator<(const Edge& other) const {
    return balanced_weight < other.balanced_weight;
}

bool Edge::CompareByBalancedWeight(const EdgePtr& e1, const EdgePtr& e2) {
    return e1->balanced_weight < e2->balanced_weight;
}

Edge::Edge(const Edge *pEdge) {
    node_a = pEdge->node_a;
    node_b = pEdge->node_b;
    weight = pEdge->weight;
    balanced_weight = pEdge->balanced_weight;
    visited_count = pEdge->visited_count;
}

void Edge::display() const {
        // printf("%s(%d) -- %s(%d)\t %f\t %f\n", node_a->protein_name.c_str(), node_a->id, \
        //         node_b->protein_name.c_str(), node_b->id, weight, balanced_weight);
    cout << node_a->protein_name << "(" << node_a->id << ") -- " \
        << node_b->protein_name << "(" << node_b->id << ")\t" \
        << weight << "\t" << balanced_weight << std::endl;
}

UnGraph::UnGraph(string& ppi_file, Weight weight) {
    // update Node
    if(weight) {
        vector<string> edge_list;
        vector<double> edge_weight;
        read_edge_list_with_weight(ppi_file, proteins_name, edge_list, edge_weight);
        assert(edge_list.size() % 2 == 0);
        ID2Protein.resize(proteins_name.size());
        connected.resize(proteins_name.size(), vector<bool>(proteins_name.size(), false));
        // update node
        for(auto& protein_name: proteins_name) {
            ProteinPtr protein(new Protein(static_cast<int>(proteins.size()), protein_name));
            ID2Protein[proteins.size()] = protein;  // ID2Protein
            Protein2ID[protein] = static_cast<int>(proteins.size());  // Protein2ID
            protein_name_id[protein_name] = static_cast<int>(proteins.size());  // protein_name_id
            proteins.insert(protein);            //  Proteins
        }
        // update edge
        assert(edge_list.size() == edge_weight.size() * 2);
        for(int i = 0; i < edge_list.size(); i += 2) {
            auto protein_a = ID2Protein[protein_name_id[edge_list[i]]];
            auto protein_b = ID2Protein[protein_name_id[edge_list[i+1]]];
            add_edge(protein_a, protein_b, edge_weight[i / 2]); // 默认权重为0
        }

    } else {
        // set<string> protein_name;
        vector<string> edge_list;
        read_edge_list(ppi_file, proteins_name, edge_list);
        assert(edge_list.size() % 2 == 0);
        ID2Protein.resize(proteins_name.size());
        connected.resize(proteins_name.size(), vector<bool>(proteins_name.size(), false));
        // update node
        for(auto& protein_name: proteins_name) {
            ProteinPtr protein(new Protein(static_cast<int>(proteins.size()), protein_name));
            ID2Protein[proteins.size()] = protein;  // ID2Protein
            Protein2ID[protein] = static_cast<int>(proteins.size());  // Protein2ID
            protein_name_id[protein_name] = static_cast<int>(proteins.size());  // protein_name_id
            proteins.insert(protein);            //  Proteins
        }

        // update edges
        for(int i = 0; i < edge_list.size(); i += 2) {
            auto protein_a = ID2Protein[protein_name_id[edge_list[i]]];
            auto protein_b = ID2Protein[protein_name_id[edge_list[i+1]]];
            add_edge(protein_a, protein_b); // 默认权重为0
        }
    }
}

UnGraph::UnGraph(set<string>&& set_proteins, vector<string>&& list_edges, vector<double>&& edge_weight) {
    assert(list_edges.size() % 2 == 0);
    ID2Protein.resize(set_proteins.size());
    connected.resize(set_proteins.size(), vector<bool>(set_proteins.size(), false));
    // update node
    for(auto& protein_name: set_proteins) {
        ProteinPtr protein(new Protein(static_cast<int>(proteins.size()), protein_name));
        ID2Protein[proteins.size()] = protein;  // ID2Protein
        Protein2ID[protein] = static_cast<int>(proteins.size());  // Protein2ID
        protein_name_id[protein_name] = static_cast<int>(proteins.size());  // protein_name_id
        proteins.insert(protein);            //  Proteins
    }

    // update edge
    if(edge_weight.empty()) {
        for(int i = 0; i < list_edges.size(); i += 2) {
            auto protein_a = ID2Protein[protein_name_id[list_edges[i]]];
            auto protein_b = ID2Protein[protein_name_id[list_edges[i+1]]];
            add_edge(protein_a, protein_b); // 默认权重为0
        } 
    } else {
        assert(edge_weight.size() * 2 == list_edges.size());
        for(int i = 0; i < list_edges.size(); i += 2) {
            auto protein_a = ID2Protein[protein_name_id[list_edges[i]]];
            auto protein_b = ID2Protein[protein_name_id[list_edges[i+1]]];
            add_edge(protein_a, protein_b, edge_weight[i / 2]);
        }
    }

}

UnGraph::~UnGraph() {
    ID2Protein.clear();
    Protein2ID.clear();
    protein_name_id.clear();
    proteins.clear();
    edges.clear();
    connected.clear();
    Edge2ID.clear();
}

void UnGraph::display() const {
    cout << "protein_name_id: " << endl;
    for(auto& it: protein_name_id) {
        cout << it.second << "\t" << it.first << endl;
    }

    cout << "ID2Protein: " << endl;
    for(int i = 0; i < ID2Protein.size(); ++i) {
        cout << i << "\t" << ID2Protein[i]->protein_name << endl;
    }

    cout << "edges: " << endl;
    for(auto& e: edges) {
        cout << e->node_a->protein_name << "\t" << e->node_b->protein_name << endl;
    }

    cout << "edges_id" << endl;
    for(auto& it: Edge2ID) {
        for(auto& p: it.first) {
            cout << p << "\t";
        }
        cout << it.second << endl;
    }
    cout << "connected: " << endl;
    for(auto& row: connected) {
        for(auto v: row) {
            cout << v << "  ";
        }
        cout << endl;
    }
}


void UnGraph::read_edge_list(string& file_path, set<string>& proeins_list, vector<string>& edge_list) {
    fstream file(file_path);
    if(!file.is_open()) {
        cout << "Failed to open file! " << file_path << endl;
        exit(1);
    }
    string line;
    while(getline(file, line)) {
        istringstream iss(line);
        string protein_name;
        while(iss >> protein_name) {
            proeins_list.insert(protein_name);
            edge_list.emplace_back(protein_name);
        }
    }
    file.close();
}

void UnGraph::read_edge_list_with_weight(string& file_path, set<string>& proeins_list, vector<string>& edge_list, vector<double>& edge_weight) {
    fstream file(file_path);
    if(!file.is_open()) {
        cout << "Failed to open file! " << file_path << endl;
        exit(1);
    }
    string line;
    while(getline(file, line)) {
        istringstream iss(line);
        string protein_name;
        // first protein
        iss >> protein_name;
        proeins_list.insert(protein_name);
        edge_list.emplace_back(protein_name);
        // second protein
        iss >> protein_name;
        proeins_list.insert(protein_name);
        edge_list.emplace_back(protein_name);
        // weight
        double weight;
        iss >> weight;
        edge_weight.emplace_back(weight);
    }
    file.close();
}

void UnGraph::add_edge(const ProteinPtr& protein_a, const ProteinPtr& protein_b, double weight) {
    if(!protein_a || !protein_b) {
        cerr << "Protein not founded!" << endl;
        return;
    }
    if(protein_a == protein_b) {
        return;
    }
    protein_a->add_neighbor(protein_b);
    protein_b->add_neighbor(protein_a);
    connected[protein_a->id][protein_b->id] = true;
    connected[protein_b->id][protein_a->id] = true;

    auto e = new Edge(protein_a, protein_b, weight, weight);
    edges.emplace_back(e);
    set<string> e_set{protein_a->protein_name, protein_b->protein_name};
    Edge2ID.insert(make_pair(move(e_set),  Edge2ID.size()));
}

EdgePtr UnGraph::getEdge(const ProteinPtr& protein_a, const ProteinPtr& protein_b) {
    set<string> e = {protein_a->protein_name, protein_b->protein_name};
    auto edge_it = Edge2ID.find(e);
    if(edge_it == Edge2ID.end()) {
        return nullptr;
    } else {
        return edges[edge_it->second];
    }
}

// 使用go term为连边添加权重
void UnGraph::weight_by_go_term(BioInformation& bio, DAG& dag) {
    for(auto & edge : edges) {
        string protein_1 = edge->node_a->protein_name;
        string protein_2 = edge->node_b->protein_name;
        set<string> go_term_1 = bio.go_slim[protein_1];
        set<string> go_term_2 = bio.go_slim[protein_2];
        double weight = dag.get_similarity_protein(go_term_1, go_term_2);
        edge->weight = weight;
        edge->balanced_weight = weight;
    }
}

// 直接算两个或者更多的节点
double UnGraph::agglomeration_coefficient(const vector<ProteinPtr>& nodes){
    int edge_count = 0;
    for(auto node1 =  nodes.begin(); node1 != nodes.end(); ++node1) {
        for(auto node2 = next(node1); node2 != nodes.end(); ++node2) {
            if(getEdge(*node1, *node2)) {
                edge_count +=1 ;
            }
        }
    }
    return static_cast<double>(edge_count * 2) / static_cast<double>((nodes.size() * (nodes.size() - 1)));
}

// 计算平衡系数重新计算边权重
void UnGraph::calculate_balanced_weight(){
    map <int,double> Sum;
    for(const auto & edge : edges) {
        Sum[edge->node_a->id] += edge->weight;
        Sum[edge->node_b->id] += edge->weight;
    }
    for(const auto & edge : edges) {
        if(edge->weight - 0 <= 0.0000001) {
            edge->balanced_weight = 0.0;
            continue;
        }
        edge->balanced_weight = pow(edge->weight, BALANCED_INDEX) \
            / pow(Sum[edge->node_a->id], BALANCED_INDEX - 1 )
                                    + pow(edge->weight, BALANCED_INDEX) \
            / pow(Sum[edge->node_b->id],BALANCED_INDEX - 1);
    }
    // for(auto& e: edges) {
    //     cout << e->node_a->protein_name << "\t" << e->node_b->protein_name << "\t" << e->balanced_weight << endl;
    // }
}

void UnGraph::calculate_structure_similarty(vector<vector<double>>& ss_weight) const {
    ss_weight.resize(ID2Protein.size(), vector<double>(ID2Protein.size(), 0.0));
    vector<vector<int>> common_neighbor_size;
    get_common_neighbor_size(common_neighbor_size);
    // 计算每条边的jcs
    vector<vector<double>> JCS;
    get_JCS(JCS, common_neighbor_size);
    // 计算cns
    vector<vector<double>> CNS;
    get_CNS(CNS, JCS);
    for(auto it1 = proteins.begin(); it1 != proteins.end(); ++it1) {
        for(auto it2 = next(it1); it2 != proteins.end(); ++it2) {
            double ss_w = (JCS[(*it1)->id][(*it2)->id] + CNS[(*it1)->id][(*it2)->id]) / (common_neighbor_size[(*it1)->id][(*it2)->id] + 1);
            ss_weight[(*it1)->id][(*it2)->id] = ss_w;
            ss_weight[(*it2)->id][(*it1)->id] = ss_w;
        }
    }
}

void UnGraph::get_common_neighbor_size(vector<vector<int>>& common_size) const {
    common_size.resize(ID2Protein.size(), vector<int>(ID2Protein.size(), 0));
    for(auto it1 = proteins.begin(); it1 != proteins.end(); ++it1) {
        for(auto it2 = next(it1); it2 != proteins.end(); ++it2) {
            vector<ProteinPtr> common;
            set_intersection((*it1)->neighbor.begin(), (*it1)->neighbor.end(),
                             (*it2)->neighbor.begin(), (*it2)->neighbor.end(),
                             inserter(common, common.begin()));
            common_size[(*it1)->id][(*it2)->id] = static_cast<int>(common.size());
            common_size[(*it2)->id][(*it1)->id] = static_cast<int>(common.size());
        }
    }
}

void UnGraph::get_JCS(vector<vector<double>>& JCS, const vector<vector<int>>& common_size) const {
    JCS.resize(ID2Protein.size(), vector<double>(ID2Protein.size(), 0));
    for(auto it1 = proteins.begin(); it1 != proteins.end(); ++it1) {
        for(auto it2 = next(it1); it2 != proteins.end(); ++it2) {
            double jcs = static_cast<double>(common_size[(*it1)->id][(*it2)->id]) / static_cast<double>(((*it1)->neighbor.size() + (*it2)->neighbor.size() - common_size[(*it1)->id][(*it2)->id]));
            JCS[(*it1)->id][(*it2)->id] = jcs;
            JCS[(*it2)->id][(*it1)->id] = jcs;
        }
    }
}

void UnGraph::get_CNS(vector<vector<double>>& CNS, const vector<vector<double>>& JCS) const {
    CNS.resize(ID2Protein.size(), vector<double>(ID2Protein.size(), 0));
    for(auto it1 = proteins.begin(); it1 != proteins.end(); ++it1) {
        for(auto it2 = next(it1); it2 != proteins.end(); ++it2) {
            vector<ProteinPtr> common;
            set_intersection((*it1)->neighbor.begin(), (*it1)->neighbor.end(),
                             (*it2)->neighbor.begin(), (*it2)->neighbor.end(),
                             inserter(common, common.begin()));
            double cns = 0.0;
            for(const auto& neighbor: common) {
                cns += JCS[(*it1)->id][neighbor->id] + JCS[(*it2)->id][neighbor->id];
            }
            CNS[(*it1)->id][(*it2)->id] = cns;
            CNS[(*it2)->id][(*it1)->id] = cns;
        }
    }
}

double UnGraph::homogeneity_index_jaccard(const ProteinPtr& a, const ProteinPtr& b) {
    return a->jaccard_similarity(b);
}

// 高度节点必然大
void UnGraph::calculate_attraction(vector<double>& attractions) const {
    attractions.resize( ID2Protein.size(), 0.0);

    vector<double> weight_node(ID2Protein.size(), 0.0);   // 每个节点对应的质量
    for(auto& e: edges) {
        weight_node[e->node_b->id] += e->balanced_weight;
        weight_node[e->node_a->id] += e->balanced_weight;
    }
    // 计算每一条边上面的吸引力GMm / (1 - sim)
    for(auto& e: edges) {
        double d = 2.0 - e->balanced_weight;  // 相异度   避免出现除以零的情况
        double attraction = weight_node[e->node_a->id] * weight_node[e->node_b->id] / pow(d, 2);
        attractions[e->node_a->id] += attraction;
        attractions[e->node_b->id] += attraction;
    }
}

void UnGraph::calculate_average_attraction(vector<double>& attractions) const {
    attractions.resize( ID2Protein.size(), 0.0);

    vector<double> weight_node(ID2Protein.size(), 0.0);   // 每个节点对应的质量
    for(auto& e: edges) {
        weight_node[e->node_b->id] += e->balanced_weight;
        weight_node[e->node_a->id] += e->balanced_weight;
    }
    // 计算每一条边上面的吸引力GMm / (1 - sim)
    for(auto& e: edges) {
        double d = 1.1 - e->balanced_weight;  // 相异度   避免出现除以零的情况
        double attraction = weight_node[e->node_a->id] * weight_node[e->node_b->id] / pow(d, 2);
        attractions[e->node_a->id] += attraction;
        attractions[e->node_b->id] += attraction;
    }

    for(int i = 0; i < attractions.size(); ++i) {
        attractions[i] /= static_cast<double>(ID2Protein[i]->neighbor.size());
    }
}

void UnGraph::calculate_walk_probability(vector<vector<double>>& probability) const {
    probability.resize(static_cast<int>(ID2Protein.size()), vector<double>(static_cast<int>(ID2Protein.size()), 0.0));

    // 计算每个节点的权重
    vector<double> node_weight(ID2Protein.size(), 0.0);
    for(auto & edge : edges) {
        node_weight[edge->node_a->id] += edge->balanced_weight;
        node_weight[edge->node_b->id] += edge->balanced_weight;
    }
    for(auto& e: edges) {
        // a ---> b 的概率
        probability[e->node_a->id][e->node_b->id] = e->balanced_weight / node_weight[e->node_a->id];
        probability[e->node_b->id][e->node_a->id] = e->balanced_weight / node_weight[e->node_b->id];
    }

}

int UnGraph::find_parent(int protein, map<int, int>& parent) {
    int root = protein;
    // find root
    while (root != parent[root]) {
        root = parent[root];
    }
    // path compression
    while (protein != root) {
        int next = parent[protein];
        parent[protein] = root;
        protein = next;
    }
    return root;
}

int getfa1(int x,int fa[]) {
    int a = x;
    while (x != fa[x]) 
	{
        x = fa[x];
    }
    while (a != fa[a]) 
	{
        int z = a;
        a = fa[a];
        fa[z] = x;
    }
    return x;
}

void UnGraph::split_graph(queue<UnGraph>& ppi_queue, vector<UnGraph>& splited_ppi, BioInformation& bio, DAG& dag) {
    UnGraph current_ppi = ppi_queue.front();
    ppi_queue.pop();

    if(current_ppi.proteins.size() <= COMPLEX_MAX_SIZE) {
        if(current_ppi.proteins.size() >= 3) {
            splited_ppi.emplace_back(move(current_ppi));
        }
        return;
    }
    
    current_ppi.weight_by_go_term(bio, dag);
    for(auto& e: current_ppi.edges) {
        e->balanced_weight += e->node_a->jaccard_similarity(e->node_b);
    }


    int count = 0, location = -1;
    // map<int, int> parent;
    // for(auto& p: current_ppi.ID2Protein) {
    //     parent.insert(make_pair(p->id, p->id));
    // }
    int fa[current_ppi.proteins.size()];
    for(int i = 0; i < current_ppi.proteins.size(); ++i) {
        fa[i] = i;
    }

    // 从小大到排列
    sort(current_ppi.edges.begin(), current_ppi.edges.end(), Edge::CompareByBalancedWeight);
    for(int i = current_ppi.edges.size() - 1; i >= 0; i--) {
        if(getfa1(current_ppi.edges[i]->node_a->id, fa) == getfa1(current_ppi.edges[i]->node_b->id, fa)) {
            continue;
        }
        fa[getfa1(current_ppi.edges[i]->node_a->id, fa)] = getfa1(current_ppi.edges[i]->node_b->id, fa);
        count += 1;
        if(count == current_ppi.proteins.size() - 2) {
            break;
        }
    }

    while(true) {
        bool success = false;
        set<string> new_protein_set;
        vector<string> new_edge_list;
        vector<double> edge_weight;
        for(int i = location + 1; i < current_ppi.ID2Protein.size(); ++i) {
            assert(i == current_ppi.ID2Protein[i]->id);
            if(getfa1(current_ppi.ID2Protein[i]->id, fa) == current_ppi.ID2Protein[i]->id) {
                location = i; 
                success = true;
                break;
            }
        }
        if(!success) {
            break;
        }
        // add new protein_name
        for(int i = 0; i < current_ppi.ID2Protein.size(); ++i) {
            if(getfa1(current_ppi.ID2Protein[i]->id, fa) == current_ppi.ID2Protein[location]->id) {
                new_protein_set.insert(current_ppi.ID2Protein[i]->protein_name);
            }
        }
        // add new edge_list;
        for(auto& e: current_ppi.edges) {
            if(new_protein_set.count(e->node_a->protein_name) && new_protein_set.count(e->node_b->protein_name)) {
                new_edge_list.emplace_back(e->node_a->protein_name);
                new_edge_list.emplace_back(e->node_b->protein_name);
                edge_weight.emplace_back(e->balanced_weight);
            }
        }
        UnGraph new_ppi(move(new_protein_set), move(new_edge_list));
        // 添加同质性指数
        // for(auto& e: new_ppi.edges) {
        //     e->weight += e->node_a->jaccard_similarity(e->node_b);
        // }
        // new_ppi.weight_by_go_term(bio, dag);
        // new_ppi.calculate_balanced_weight();
        ppi_queue.push(move(new_ppi));
    }

}

void UnGraph::get_complexes(UnGraph& origin_ppi, UnGraph& g, vector<set<string>>& complexes, double similarity_threshold) {
    if(g.proteins.empty()) return;
    map<int, bool> complex_record;
    queue<int> complex_queue;

    for(auto& protein: g.ID2Protein) {
        complex_queue.push(1 << protein->id);
        complex_record[1 << protein->id] = true;
    }

    while(!complex_queue.empty()) {
        int complex_state = complex_queue.front();
        complex_queue.pop();
        for(int i = 0; i < g.ID2Protein.size(); ++i) {
            if((complex_state & (1 << i)) == 0) {
                continue;
            }
            for(int j = 0; j < g.ID2Protein.size(); ++j) {
                if(!g.connected[i][j]) {
                    continue;
                }
                if((complex_state & (1 << j)) != 0) {
                    continue;
                }
                if(complex_record.count(complex_state | (1 << j)) != 0) {
                    continue;
                }
                complex_record[complex_state | (1 << j)] = true;
                complex_queue.push(complex_state | (1 << j));
            }
        }
    }
    set<string> complex;
    auto it = complex_record.begin();
    vector<set<string>> complex_result;
    while(it != complex_record.end()) {
        complex.clear();
        for(int i = 0; i < g.ID2Protein.size(); ++i) {
            if((it->first & (1 << i)) != 0) {
                complex.insert(g.ID2Protein[i]->protein_name);
            }
        }
        if(complex.size() >= 3) {
            // Complex c(g, move(complex));
            // complex_result.insert(complex_result.end(), move(c));
            if(Complex::evaluate_by_weight(origin_ppi, complex)) {
                complex_result.emplace_back(move(complex));
            }
        }
        ++it;
    }
    if(complex_result.empty()) {
        return;
    }
    sort(complexes.begin(), complexes.end(), Complex::CompareSetBySize);
    for(auto& c: complex_result) {
        Complex::update_complexes(complexes, c);
    }
}

void UnGraph::get_complexes1(UnGraph& g, vector<set<string>>& complexes, double similarity_threshold) {
    set<string> complex;
    for(auto& p: g.proteins) {
        complex.clear();
        complex.insert(p->protein_name);
        for(auto& nei: p->neighbor) {
            complex.insert(nei->protein_name);
        }
        complexes.emplace_back(complex);
    }
}

vector<double> UnGraph::calculate_protein_weight() const {
    vector<double> protein_weight(ID2Protein.size(), 0.0);

    for(auto& e: edges) {
        protein_weight[e->node_a->id] += e->balanced_weight;
        protein_weight[e->node_b->id] += e->balanced_weight;
    }

    return move(protein_weight);
}

void UnGraph::write_to_file(const string& file_path) const {
    ofstream file(file_path);
    if(!file.is_open()) {
        cerr << "<ERROR>: Failed to open file! " << file_path << endl;
        exit(1);
    }
    for(auto& e: edges) {
        file << e->node_a->protein_name << "\t" << e->node_b->protein_name << "\t" << e->balanced_weight << endl;
    }
    file.close();
    cout << "Write to file finished!" << endl;
}


bool UnGraph::compare_pairs(const pair<EdgePtr, int>& pair1, const pair<EdgePtr, int>& pair2) {
    return pair1.second > pair2.second;
}

Complex::Complex(const UnGraph& g, vector<string>& _proteins) {
    set<string> protein_set(_proteins.begin(), _proteins.end());
    proteins = move(_proteins);
    cohesion = calculate_cohesion(g, protein_set);
}

double Complex::complex_match_score(Complex& other) {
    set<string> common;
    set_intersection(proteins.begin(), proteins.end(),
                          other.proteins.begin(), other.proteins.end(),
                          inserter(common, common.begin()));
    double score = static_cast<double>(common.size()) / static_cast<double>(max(other.proteins.size(), proteins.size()));
    display();
    other.display();
    cout << score << endl;
    return score;
}

void Complex::display() const {
    cout << cohesion << "\t";
    for(auto& p: proteins) {
        cout << p << "\t";
    }
    cout << endl;
}

bool Complex::CompareBySize(const Complex& a, const Complex& b) {
    return a.proteins.size() > b.proteins.size();
}

bool Complex::CompareSetBySize(const set<string>& a, const set<string>& b)  {
    return a.size() > b.size();        
}

set<set<string>> Complex::read_complex_from_file(string& file_path) {
    set<set<string>> complexes;
    fstream file(file_path);
    if(!file.is_open()) {
        cout << "<ERROR>: Failed to open file! " << file_path << endl;
        exit(1);
    }

    string line;
    while(getline(file, line)) {
        set<string> complex;
        istringstream iss(line);
        string protein;
        while(iss >> protein) {
            complex.insert(protein);
        }
        complexes.insert(move(complex));
    }
    file.close();
    return move(complexes);
}

void Complex::write_complex_to_file(vector<Complex>& complexes, string& file_path) {
    ofstream file(file_path);
    if(!file.is_open()) {
        cerr << "<ERROR>: Failed to open file! " << file_path << endl;
        exit(1);
    }
    for(auto& c: complexes) {
        for(auto& p: c.proteins) {
            file << p << "\t";
        }
        file << endl;
    } 
    file.close();
}


void Complex::write_complex_to_file(vector<set<string>>& complexes, string& file_path) {
    ofstream file(file_path);
    if(!file.is_open()) {
        cerr << "<ERROR>: Failed to open file! " << file_path << endl;
        exit(1);
    }
    for(auto& c: complexes) {
        for(auto& p: c) {
            file << p << "\t";
        }
        file << endl;
    } 
    file.close();
}

void Complex::write_complex_to_file(set<set<string>>& complexes, string& file_path) {
    ofstream file(file_path);
    if(!file.is_open()) {
        cerr << "<ERROR>: Failed to open file! " << file_path << endl;
        exit(1);
    }
    for(auto& c: complexes) {
        for(auto& p: c) {
            file << p << "\t";
        }
        file << endl;
    }
    file.close();
}

// need to modify
double Complex::calculate_cohesion(const UnGraph& g, set<string>& _proteins) {
    map<string, double> sum;
    map<string, int> count;

    for(auto& e: g.edges) {
        if(!_proteins.count(e->node_a->protein_name) || !_proteins.count(e->node_b->protein_name)) {
            continue;
        }
        count[e->node_a->protein_name] += 1;
        count[e->node_b->protein_name] += 1;
        sum[e->node_a->protein_name] += e->balanced_weight;
        sum[e->node_b->protein_name] += e->balanced_weight;
    }
    double temp_cohesion = 0.0;
    for(auto& p: _proteins) {
        temp_cohesion += sum[p] * (count[p] + 1) / static_cast<int>(_proteins.size());
    }
    double cohesion = temp_cohesion / static_cast<double>(_proteins.size());
    cout << "complex cohesion: " << cohesion << endl;
    return cohesion;
}

void Complex::update_complexes(vector<Complex>& complexes, Complex& temp_complex) {
    sort(complexes.begin(), complexes.end(), CompareBySize);
    for(auto& complex: complexes) {
        if(complex.complex_match_score(temp_complex) >=MAX_MATCH_INDEX ) {
            return;
        }       
    }
    complexes.emplace_back(move(temp_complex));
}

bool CompareSet(const set<string>& s1, const set<string>& s2) {
    return s1.size() > s2.size();
}

double complex_match_score1(const set<string>& s1, const set<string>& s2) {
    set<string> common;
    set_intersection(s1.begin(), s1.end(),
                     s2.begin(), s2.end(),
                     inserter(common, common.begin()));

    return static_cast<double>(common.size()) / static_cast<double>(max(s1.size(), s2.size()));
}

void Complex::update_complexes(vector<set<string>>& complexes, set<string>& complex) {
    // std::cout << complex.size() << endl;
    sort(complexes.begin(), complexes.end(), CompareSet);
    for(auto& c: complexes){
        if(complex_match_score1(c, complex) >= MAX_MATCH_INDEX) {
            return;
        }
    }
    complexes.emplace_back(move(complex));
}

bool Complex::evaluate_by_weight(UnGraph& g, set<string>& complex) {
    double weight_in = 0.0;
    double weight_out = 0.0;
    for(auto it_p = complex.begin(); it_p != complex.end(); ++it_p) {
        auto p_ptr = g.ID2Protein[g.protein_name_id[*it_p]];
        for(auto& nei: p_ptr->neighbor) {
            auto e = g.getEdge(p_ptr, nei);
            if(e == nullptr) {
                continue;
            }
            if(complex.count(nei->protein_name)) {
                weight_in += e->balanced_weight;
            } else {
                weight_out += e->balanced_weight;
            }
        }
    }
    return weight_in * 2.2 >= weight_out;
}


void UnGraph::split_graph1(queue<UnGraph>& ppi_queue, vector<UnGraph>& splited_ppi, BioInformation& bio, DAG& dag) {
    UnGraph current_ppi = ppi_queue.front();
    ppi_queue.pop();

    if(current_ppi.proteins.size() <= COMPLEX_MAX_SIZE) {
        if(current_ppi.proteins.size() >= 3) {
            splited_ppi.emplace_back(move(current_ppi));
        }
        return;
    }

    int count = 0;
    int fa[current_ppi.proteins.size()];
    for(auto& p: current_ppi.proteins) {
        fa[p->id] = p->id;
    }

    // 从小大到排列
    sort(current_ppi.edges.begin(), current_ppi.edges.end(), Edge::CompareByBalancedWeight);
    std::reverse(current_ppi.edges.begin(), current_ppi.edges.end());
    for(auto& e: current_ppi.edges) {
        if(getfa1(e->node_a->id, fa) == getfa1(e->node_b->id, fa)) {
            continue;
        }
        fa[getfa1(e->node_a->id, fa)] = getfa1(e->node_b->id, fa);
        count += 1;
        if(count + 2 == current_ppi.proteins.size()) {
            break;
        }
    }

    set<ProteinPtr>::iterator location;
    while(true) {
        bool success = false;
        set<string> new_protein_set;
        vector<string> new_edge_list;
        vector<double> edge_weight;
        for(auto it = current_ppi.proteins.begin(); it != current_ppi.proteins.end(); ++it) {
            if(getfa1((*it)->id, fa) == (*it)->id) {
                location = it; 
                success = true;
                break;
            }
        }
        if(!success) {
            break;
        }
        // add new protein_name
        for(auto it = current_ppi.proteins.begin(); it != current_ppi.proteins.end(); ++it) {
            if(getfa1((*it)->id, fa) == (*location)->id) {
                new_protein_set.insert((*it)->protein_name);
            }
        }
        // add new edge_list;
        for(auto& e: current_ppi.edges) {
            if(new_protein_set.count(e->node_a->protein_name) && 
                    new_protein_set.count(e->node_b->protein_name)) {
                
                new_edge_list.emplace_back(e->node_a->protein_name);
                new_edge_list.emplace_back(e->node_b->protein_name);
                edge_weight.emplace_back(e->balanced_weight);
            }
        }
        UnGraph new_ppi(move(new_protein_set), move(new_edge_list), move(edge_weight));
        ppi_queue.push(move(new_ppi));
    }
}

void UnGraph::get_complexes1(UnGraph& origin_ppi, UnGraph& g, vector<set<string>>& complexes, double similarity_threshold) {
    if(g.proteins.empty()) return;
    map<int, bool> complex_record;
    queue<int> complex_queue;

    for(auto& protein: g.ID2Protein) {
        complex_queue.push(1 << protein->id);
        complex_record[1 << protein->id] = true;
    }

    while(!complex_queue.empty()) {
        int complex_state = complex_queue.front();
        complex_queue.pop();
        for(int i = 0; i < g.ID2Protein.size(); ++i) {
            if((complex_state & (1 << i)) == 0) {
                continue;
            }
            for(int j = 0; j < g.ID2Protein.size(); ++j) {
                if(!g.connected[i][j]) {
                    continue;
                }
                if((complex_state & (1 << j)) != 0) {
                    continue;
                }
                if(complex_record.count(complex_state | (1 << j)) != 0) {
                    continue;
                }
                complex_record[complex_state | (1 << j)] = true;
                complex_queue.push(complex_state | (1 << j));
            }
        }
    }
    set<string> complex;
    auto it = complex_record.begin();
    vector<set<string>> complex_result;
    while(it != complex_record.end()) {
        complex.clear();
        for(int i = 0; i < g.ID2Protein.size(); ++i) {
            if((it->first & (1 << i)) != 0) {
                complex.insert(g.ID2Protein[i]->protein_name);
            }
        }
        if(complex.size() >= 3) {
            if(Complex::evaluate_by_weight(origin_ppi, complex)) {
                complex_result.emplace_back(move(complex));
            }
        }
        ++it;
    }
    if(complex_result.empty()) {
        return;
    }
    sort(complexes.begin(), complexes.end(), Complex::CompareSetBySize);
    for(auto& c: complex_result) {
        Complex::update_complexes(complexes, c);
    }
}
