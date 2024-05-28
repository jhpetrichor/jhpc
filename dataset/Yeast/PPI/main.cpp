/*
 * @brief: 
 * @Author: jh
 * @Date: 2024-05-14 21:16:08
 * @LastEditTime: 2024-05-14 22:16:55
 */
// 整理protein complex 数据集
// 将CYC2008 和 sgd数据集合二为一
// 删除大小为1的蛋白质复合物
#include <cstdio>
#include <iostream>
#include <iterator>
#include <set>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

void read_protein(string file_path, set<string>& proteins) {
    std::fstream file(file_path);
    if(!file.is_open()) {
        cerr << "Faile to open file! " << file_path << endl;
        exit(1);
    }

    string line;
    while(getline(file, line)) {
        istringstream iss(line);
        string protein;
        while(iss >> protein) {
            proteins.insert(protein);
        }
    }
    file.close();
}

set<string> read_go_slim(string file_path) {
    std::fstream file(file_path);
    if(!file.is_open()) {
        cerr << "Faile to open file! " << file_path << endl;
        exit(1);
    }

    set<string> proteins;

    string line;
    while(getline(file, line)) {
        istringstream iss(line);
        string protein;
        iss >> protein;
        proteins.insert(protein);
    }
    file.close();
    return proteins;
}

void writet_to_file(string file_path, set<string>& proteins) {
    std::ofstream file(file_path);
    for(auto& p: proteins) {
        file << p << "\n";
    }
    file.close();
}

int main() {
    set<string> proteins;
    read_protein("collins.txt", proteins);
    read_protein("gavin.txt", proteins);
    read_protein("krogan_core.txt", proteins);
    read_protein("krogan_extended.txt", proteins);
    read_protein("biogrid.txt", proteins);
    read_protein("yeast_all.txt", proteins);
    read_protein("original_dip.txt", proteins);
    read_protein("original_krogan.txt", proteins);
    read_protein("original_dip.txt", proteins);
    read_protein("original_Gavin.txt", proteins);
    read_protein("original_MIPS.txt", proteins);



    std::cout << "proteins size: " << proteins.size() << endl;
    auto exit_protein = read_go_slim("../DAG/go-slim.txt");
    auto exit_protein_expression = read_go_slim("../DAG/gene-expression.txt");
    set<string> both_exit;
    set_intersection(exit_protein.begin(), exit_protein.end(),
                    exit_protein_expression.begin(), exit_protein_expression.end(),
                    inserter(both_exit, both_exit.begin()));

    std::cout << "exit protein size: " << both_exit.size() << endl;

    set<string> no_exit;
    std::set_difference(proteins.begin(), proteins.end(),
                        both_exit.begin(), both_exit.end(), 
                        std::inserter(no_exit, no_exit.begin()));
    std::cout << "no_exit: " << no_exit.size() << endl;
    writet_to_file("./no_exi1t.txt", no_exit);
    return 0;
}