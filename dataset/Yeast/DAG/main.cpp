/*
 * @brief: 
 * @Author: jh
 * @Date: 2024-05-06 13:23:07
 * @LastEditTime: 2024-06-01 16:02:47
 */
#include <cstdio>
#include <iostream>
#include <set>
#include <fstream>
#include <ostream>
#include <sstream>
#include <vector>
#include <regex>

using namespace std;

int main() {
    std::fstream file("./is_a.txt");
    string line;
    map<string, set<string>> child_ancestor;
    
    while(getline(file, line)) {
        istringstream iss(line);
        string ancestor;
        iss >> ancestor;
        string child;
        while(iss >> child) {
            child_ancestor[child].insert(ancestor);
        }
    }
    file.close();

    std::ofstream ofile("./is_a_child_ancestor.txt");
    for(auto& it: child_ancestor) {
        ofile << it.first << "\t";
        for(auto& a: it.second) {
            ofile << a << "\t";
        }
        ofile << endl;
    }
    ofile.close();

    return 0;
}