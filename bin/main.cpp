/*
 * @brief: 
 * @Author: jh
 * @Date: 2024-05-20 21:38:38
 * @LastEditTime: 2024-05-22 16:26:12
 */
#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace std;


int main(int argc, char** argv) {
    // read
    std::fstream nfile("/home/jh/code/JHPC/dataset/Yeast/PPI/no_exi1t.txt");
    if(!nfile.is_open()) {
        cerr << "open no_exit failed!" << endl;
    }
    string line;
    set<string> no_exit;
    while(getline(nfile, line)) {
        istringstream iss(line);
        string protein;
        while(iss >> protein) {
            no_exit.insert(protein);
        }
    }

    set<set<string>> complexes;
    std::fstream file(argv[1]);
    if (!file.is_open()) {
        cerr << "Failed to open complexes file! " << endl;
        exit(1);
    }

    vector<string> edge_list;
    while(getline(file, line)) {
        string a, b;
        istringstream iss(line);
        iss >> a >> b;
        if(no_exit.count(a) || no_exit.count(b)) {
            continue;
        }
        edge_list.emplace_back(a);
        edge_list.emplace_back(b);
    }

    std::ofstream ofile(argv[2]);
    for(int i = 0; i < edge_list.size(); i+=2) {
        ofile << edge_list[i] << "\t" << edge_list[i+1] << endl;
    }

    return 0;
}
