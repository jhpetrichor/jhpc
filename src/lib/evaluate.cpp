/*
 * @brief: 
 * @Author: jh
 * @Date: 2024-05-06 15:18:33
 * @LastEditTime: 2024-06-16 20:59:29
 */
#include <iostream>
#include <vector>
#include <set>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>

using namespace std;

#define MIN_SIZE  3
#define MAX_SIZE  100

set<string> read_protein(string& ppi_file) {
    set<string> protein_set;
    fstream file(ppi_file);
    if(!file.is_open()) {
        cout << "open file error" << ppi_file << endl;
        exit(1);
    }
    string line;
    while(getline(file, line)) {
        istringstream iss(line);
        string protein1, protein2;
        iss >> protein1 >> protein2;
        protein_set.insert(protein1);
        protein_set.insert(protein2);
    }
    cout << "protein num: " << protein_set.size() << endl;
    return protein_set;
}

vector<vector<string>> read_complex(string& file_path, set<string>& protein_set) {
    vector<vector<string>> complexes;
    fstream file(file_path);
    if(!file.is_open()) {
        cout << "open file error" << file_path << endl;
        exit(1);
    }
    string line;
    while(getline(file, line)) {
        vector<string> temp_complex;
        istringstream iss(line);
        string protein;
        while(iss >> protein) {
            if(!protein_set.count(protein)) break;

            temp_complex.push_back(protein);
        }
        if(temp_complex.size() >= MIN_SIZE && temp_complex.size() <= MAX_SIZE )
            complexes.push_back(temp_complex);
    }
    return complexes;
}

double OS(vector<string>& a, vector<string>& b) {
    vector<string> common;
    set_intersection(a.begin(), a.end(), b.begin(), b.end(), inserter(common, common.end()));
    return static_cast<double>(pow(common.size(), 2)) / static_cast<double>(a.size() * b.size());
}

int calculate_NPC(vector<vector<string>>& pre_complexes, vector<vector<string>>& ref_complexes, double threshold = 0.25) {
    int count = 0;
    for(auto& pre : pre_complexes) {
        for(auto& ref : ref_complexes) {
            if(OS(pre, ref) >= threshold) {
                count += 1;
                break;
            }
        }
    }
    return count;
}

int calculate_NKC(vector<vector<string>>& pre_complexes, vector<vector<string>>& ref_complexes, double threshold = 0.25) {
    int count = 0;
    for(auto& ref: ref_complexes) {
        for(auto& pre: pre_complexes) {
            if(OS(ref, pre) >= threshold) {
                count += 1;
                break;
            }
        }
    }

    return count;
}

// int main() {
//     string ppi_file = "/home/jh/code/JHPC/dataset/Yeast/PPI/Collins.txt";
//     string ref_complex_file = "/home/jh/code/JHPC/dataset/Yeast/complex/CYC2008.txt";
//     string pre_complex_file = "/home/jh/code/JHPC/bin/collins.txt";

//     set<string> protein_set = read_protein(ppi_file);
//     vector<vector<string>> pre_complexes = read_complex(pre_complex_file, protein_set);
//     vector<vector<string>> ref_complexes = read_complex(ref_complex_file, protein_set);
//     cout << "pre complex num: " << pre_complexes.size() << endl;
//     cout << "ref complex num: " << ref_complexes.size() << endl;

//     int NKC = calculate_NKC(pre_complexes, ref_complexes, 0.25);
//     int NPC = calculate_NPC(pre_complexes, ref_complexes, 0.25);
//     cout << "NPC: " << NPC << endl;
//     cout << "NKC: " << NKC << endl;

//     double precision = static_cast<double>(NPC) / static_cast<double>(pre_complexes.size());
//     double recall = static_cast<double>(NPC) / static_cast<double>(ref_complexes.size());
//     cout << "precision: " << precision << endl;
//     cout << "recall: " << recall << endl;
//     cout << "F1: " << 2 * precision * recall / (precision + recall) << endl;
//     return 0;
// }