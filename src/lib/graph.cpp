#include "graph.h"

#include <cmath>
#include <fstream>
#include <ostream>
#include <sstream>

Edge::Edge(int a, int b, double __weight) {
    if (a < b) {
        smllar = a;
        bigger = b;
    } else {
        smllar = b;
        bigger = a;
    }
    weight = __weight;
}

bool Edge::operator==(const Edge& other) const {
    return smllar == other.smllar && bigger == other.bigger;
}

bool Edge::operator<(const Edge& other) const { return weight < other.weight; }

bool Edge::CompareByWeight(const Edge& e1, const Edge& e2) {
    return e1.weight > e2.weight;
}

void Edge::debug() const {
    printf("smarr: %d\t bigger: %d\tweight: %lf\t\n", smllar, bigger, weight);
}

//=======================Graph=============================

Graph::Graph(string& ppi_file) {
    vector<int> edge_list;
    map<string, int> __protein_id;
    map<int, string> __id_protein;
    read_edge(ppi_file, edge_list, __protein_id, __id_protein);
    id_protein = std::move(__id_protein);
    protein_id = std::move(protein_id);
    node_count = __protein_id.size();
    connected.resize(node_count, vector<bool>(node_count, false));
    for (int i = 0; i < edge_list.size(); i += 2) {
        add_edge(edge_list[i], edge_list[i + 1], 0.0);
    }
}

Graph::Graph(vector<string>& edge_list, vector<double>* _edge_weight) {
    // update protein_id and id_protein
    vector<int> __edge_list;
    for (auto& p : edge_list) {
        int id = 0;
        if (protein_id.count(p)) {
            id = protein_id.find(p)->second;
        } else {
            id = protein_id.size();
            protein_id.insert(make_pair(p, id));
            id_protein.insert(make_pair(id, p));
        }
        __edge_list.emplace_back(id);
    }
    node_count = protein_id.size();
    connected.resize(node_count, vector<bool>(node_count, false));

    // add edges
    if (_edge_weight == nullptr) {
        for (int i = 0; i < __edge_list.size(); i += 2) {
            add_edge(__edge_list[i], __edge_list[i + 1]);
        }
    } else if (_edge_weight->size() << 1 == edge_list.size()) {
        for (int i = 0; i < __edge_list.size(); i += 2) {
            add_edge(__edge_list[i], __edge_list[i + 1], i / 2);
        }
    } else {
        std::cerr << "边和权重数量不匹配! " << endl;
        exit(1);
    }
}

set<int> Graph::get_neighbor(int n) const {
    auto it = adjancy_list.find(n);
    return it == adjancy_list.end() ? set<int>() : it->second;
}

void Graph::add_edge(int u, int v, double weight) {
    // Resize the connected array if necessary
    if (u >= connected.size() || v >= connected.size()) {
        connected.resize(std::max(u, v) + 1,
                         std::vector<bool>(connected.size(), false));
    }
    connected[u][v] = true;
    connected[v][u] = true;

    // Update the adjacency list
    auto it_u = adjancy_list.find(u);
    if (it_u == adjancy_list.end()) {
        adjancy_list[u] = {v};
    } else {
        it_u->second.insert(v);
    }

    auto it_v = adjancy_list.find(v);
    if (it_v == adjancy_list.end()) {
        adjancy_list[v] = {u};
    } else {
        it_v->second.insert(u);
    }

    // Add the new edge
    Edge e(u, v, weight);
    edges.emplace_back(e);

    // Update the edge_id map
    auto edge_key = std::make_pair(u, v);
    auto it = edge_id.find(edge_key);
    if (it == edge_id.end()) {
        edge_id[edge_key] = edge_count;
        edge_count += 1;
    }
}

void Graph::remove_edge(int u, int v) {
    // 邻接矩阵删除
    connected[u][v] = false;
    connected[v][u] = false;

    // 邻接表删除
    auto it = adjancy_list.find(u);
    it->second.erase(v);
    it = adjancy_list.find(v);
    it->second.erase(u);

    // edge_id
    int offset = edge_id.find({u, v})->second;
    edges.erase(edges.begin() + offset);
    edge_count -= 1;

    update_edge_id();
    // 不更新节点个数了。如果节点度为1，鼓励节点，注意忽略计算
}

int Graph::get_edge_id(int u, int v) const {
    auto it = edge_id.find({u, v});
    if (it != edge_id.end()) {
        return it->second;
    } else {
        return -1;
    }
}

int Graph::degree(int n) const {
    if (n >= node_count) {
        cerr << "节点编号不合法！" << endl;
        exit(1);
    }
    auto it = adjancy_list.find(n);
    if (it != adjancy_list.end()) {
        return it->second.size();
    } else {
        return 0;
    }
}

void Graph::read_edge(string& file_path, vector<int>& edge_list,
                      map<string, int>& __protein_id,
                      map<int, string>& __id_protein) const {
    fstream file(file_path);
    if (!file.is_open()) {
        cerr << "Failed to open file! " << file_path << endl;
        exit(1);
    }

    string line;
    while (getline(file, line)) {
        istringstream iss(line);
        string protein;
        // first one
        iss >> protein;
        if (!__protein_id.count(protein)) {
            __protein_id.insert(make_pair(protein, __protein_id.size()));
            __id_protein.insert(make_pair(__id_protein.size(), protein));
        }
        edge_list.emplace_back(__protein_id[protein]);

        // second one
        iss >> protein;
        if (!__protein_id.count(protein)) {
            __protein_id.insert(make_pair(protein, __protein_id.size()));
            __id_protein.insert(make_pair(__id_protein.size(), protein));
        }
        edge_list.emplace_back(__protein_id[protein]);
    }
    file.close();
}

set<int> Graph::get_common_neighbor(int u, int v) const {
    auto it_u = adjancy_list.find(u);
    auto it_v = adjancy_list.find(v);
    if (it_u == adjancy_list.end() || it_v == adjancy_list.end()) {
        return set<int>();
    }
    set<int> common_neighbors;
    set_intersection(it_u->second.begin(), it_u->second.end(),
                     it_v->second.begin(), it_v->second.end(),
                     inserter(common_neighbors, common_neighbors.end()));
    return common_neighbors;
}

void Graph::weighted_by_go_term() {
    for (auto& e : edges) {
        string smaller = id_protein[e.smllar];
        string bigger = id_protein[e.bigger];
        e.weight =
            0.5 * (child_ancestor.get_similarity_protein(smaller, bigger) +
                   ancestor_child.get_similarity_protein(smaller, bigger));
    }
}

void Graph::weighted_by_go_term(vector<Edge>& __edges) const {
    for (auto& e : __edges) {
        string smaller = id_protein.find(e.smllar)->second;
        string bigger = id_protein.find(e.bigger)->second;
        e.weight =
            0.5 * (child_ancestor.get_similarity_protein(smaller, bigger) +
                   ancestor_child.get_similarity_protein(smaller, bigger));
    }
}

// void Graph::weighted_by_go_term_Lin(){
//         for (auto& e : edges) {
//         string smaller = id_protein.find(e.smllar)->second;
//         string bigger = id_protein.find(e.bigger)->second;
//         e.weight = child_ancestor.get_similarity_protein_Lin(smaller,
//         bigger);
//             // 0.5 * (child_ancestor.get_similarity_protein_Lin(smaller,
//             bigger) +
//             //        ancestor_child.get_similarity_protein_Lin(smaller,
//             bigger));
//     }
// }

void Graph::normalize_edge_weight_min_max() {
    double max_weight = 0.0;
    double min_weight = 999999;
    for (auto& e : edges) {
        if (e.weight > max_weight) {
            max_weight = e.weight;
        }
        if (e.weight < min_weight) {
            min_weight = e.weight;
        }
    }

    for (auto& e : edges) {
        e.weight = (e.weight - min_weight) / (max_weight - min_weight);
    }
}

void Graph::update_edge_id() {
    edge_id.clear();
    for (int i = 0; i < edges.size(); i++) {
        edge_id[std::make_pair(edges[i].smllar, edges[i].bigger)] = i;
    }
}

void Graph::calculate_balanced_weight(double balanced_index) {
    vector<double> sum(node_count, 0.0);
    for (auto& e : edges) {
        sum[e.smllar] += e.weight;
        sum[e.bigger] += e.weight;
    }

    for (auto& e : edges) {
        if (e.weight - 0.0 <= -0.001) {
            e.weight = 0.0;
            continue;
        }
        double temp = pow(e.weight, balanced_index);
        e.weight = temp / pow(sum[e.smllar], balanced_index - 1.0) +
                   temp / pow(sum[e.bigger], balanced_index - 1.0);
    }
}

double Graph::jaccard_similarity(int u, int v) const {
    auto it_u = adjancy_list.find(u);
    auto it_v = adjancy_list.find(v);
    set<int> common_neighbors;
    set_intersection(it_u->second.begin(), it_u->second.end(),
                     it_v->second.begin(), it_v->second.end(),
                     inserter(common_neighbors, common_neighbors.end()));
    return static_cast<double>(common_neighbors.size()) /
           static_cast<double>(it_u->second.size() + it_v->second.size() -
                               common_neighbors.size());
}

double Graph::jaccard_similarity_more(int u, int v) const {
    auto it_u = adjancy_list.find(u);
    auto it_v = adjancy_list.find(v);
    set<int> common_neighbors;
    set_intersection(it_u->second.begin(), it_u->second.end(),
                     it_v->second.begin(), it_v->second.end(),
                     inserter(common_neighbors, common_neighbors.end()));
    if (connected[u][v]) {
        return static_cast<double>(common_neighbors.size() + 2) /
               static_cast<double>(it_u->second.size() + it_v->second.size() -
                                   common_neighbors.size());
    } else {
        return static_cast<double>(common_neighbors.size()) /
               static_cast<double>(it_u->second.size() + it_v->second.size() -
                                   common_neighbors.size());
    }
}

pair<double, double> Graph::unidirectional_similarity(int u, int v) const {
    set<int> common = get_common_neighbor(u, v);
    const double u_v = static_cast<double>(common.size()) / degree(v);
    const double v_u = static_cast<double>(common.size()) / degree(u);
    return pair<double, double>(u_v, v_u);
}

double Graph::density(set<int>& nodes) const {
    int count = 0;
    for(auto i: nodes) {
        for(auto j: nodes) {
            if (i == j) {
                continue;
            }
            if (connected[i][j]) {
                count++;
            }
        }
    }
    return static_cast<double>(count) /
           static_cast<double>(nodes.size() * (nodes.size() - 1));
}

double Graph::density(set<int>& nodes1, set<int>& nodes2) const {
    int count = 0;
    for(auto i: nodes1) {
        for(auto j: nodes2) {
            if(connected[i][j]) {
                count += 1;
            }
        }
    }
    return static_cast<double>(count) /
           static_cast<double>(nodes1.size() * nodes2.size());
}
double Graph::evaluate(set<int>& complex) const {
    set<int> comunity_neighbor;
    for(auto i: complex) {
        auto i_neighbor = get_neighbor(i);
        for(auto j: i_neighbor) {
            if(complex.count(j)) {
                continue;
            }
            comunity_neighbor.insert(j);
        }
    }

    double in_density = density(complex);
    double out_density = density(complex, comunity_neighbor);
    return in_density / out_density;
}

/*======================== namespace Complex ================================ */

void Complex::update_complexes(vector<set<string>>& complexes,
                               set<string>& complex, double threshold) {
    for (auto& c : complexes) {
        if (overlapping_score(c, complex) >= threshold) {
            return;
        }
    }
    complexes.push_back(complex);
}

void Complex::write_complex(const string& file_path,
                            const map<int, string>& id_protein,
                            const vector<set<int>>& complexes) {
    ofstream file(file_path);
    if (!file.is_open()) {
        cerr << "<ERROR>: Failed to open file! " << file_path << endl;
        exit(1);
    }
    for (auto& c : complexes) {
        for (auto& p : c) {
            file << id_protein.find(p)->second << "\t";
        }
        file << endl;
    }
    file.close();
}

void Complex::write_complex(const string& file_path,
                            const set<set<string>>& complexes) {
    ofstream file(file_path);
    if (!file.is_open()) {
        cerr << "<ERROR>: Failed to open file! " << file_path << endl;
        exit(1);
    }
    for (const auto& c : complexes) {
        for (const auto& p : c) {
            file << p << "\t";
        }
        file << endl;
    }
    file.close();
}

void Complex::write_complex(const string& file_path,
                            const vector<set<string>>& complexes) {
    ofstream file(file_path);
    if (!file.is_open()) {
        cerr << "<ERROR>: Failed to open file! " << file_path << endl;
        exit(1);
    }
    for (const auto& c : complexes) {
        for (const auto& p : c) {
            file << p << "\t";
        }
        file << endl;
    }
    file.close();
}
