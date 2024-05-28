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

bool CompareBySize(const set<string>& a, const set<string>& b) {
    return a.size() > b.size();
}

double matched_score(const set<string>& a, const set<string>& b) {
    set<string> comon;
    set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                     inserter(comon, comon.begin()));
    return static_cast<double>(comon.size()) /
           static_cast<double>(max(a.size(), b.size()));
}

int main() {
    // read
    set<set<string>> complexes;
    std::fstream file("./spin_collins.txt");
    if (!file.is_open()) {
        cerr << "Failed to open complexes file! " << endl;
        exit(1);
    }

    string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        set<string> complex;
        string protein;
        while (iss >> protein) {
            complex.insert(protein);
        }
        complexes.insert(complex);
    }
    std::cout << complexes.size() << endl;
    vector<set<string>> complex_list(complexes.begin(), complexes.end());
    sort(complex_list.begin(), complex_list.end(), CompareBySize);
    vector<set<string>> result;
    for (auto& c : complex_list) {
        for (auto& d : result) {
            bool insert = true;
            double score = matched_score(d, c);
            if (score >= 0.45) {
                insert = false;
                break;
            }
            if (insert) {
                result.push_back(d);
            }
        }
    }

    std::ofstream ofile("./last.txt");

    return 0;
}
