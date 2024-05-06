/*
 * @Description: 
 * @Author: jh
 * @Date: 2024-04-30 17:11:48
 * @LastEditTime: 2024-05-01 22:07:55
 */
#include "bio_information.h"

BioInformation::BioInformation() {
    read_gene_expression();
    read_go_slim();
    read_subcellular();   
}

void BioInformation::read_gene_expression() {
    fstream file(GENE_EXPRESSION);
    if(!file.is_open()) {
        cout << "--- [ERROR]: Failed to open file: " << GENE_EXPRESSION << " ---" << endl;
        exit(1);
    }

    string line;
    while(getline(file, line)) {
        istringstream iss(line);
        
        string protein;
        iss >> protein;
        
        vector<double> arr(36, 0);
        for(int i = 0; i < 36; ++i) {
            iss >> arr[i];
        }
        
        go_expression.insert(make_pair(protein, arr));
    }    

    file.close();
}

void BioInformation::print_go_slim() {
    for(auto const& it: go_slim) {
        cout << it.first << " ";
        for(auto const& go: it.second) {
            cout << go << " ";
        }
        cout << endl;
    }
}

void BioInformation::read_go_slim() {
    fstream file(GO_SLIM);
    if(!file.is_open()) {
        cout << "--- [ERROR]: Failed to open file: " << GO_SLIM << " ---" << endl;
        exit(1);
    }
    string line;
    while(getline(file, line)) {
        istringstream iss(line);
        string protein;
        iss >> protein;
        set<string> gos;
        string s;
        while(iss >> s) {
            gos.insert(s);
        }
        go_slim.insert(make_pair(protein, gos));
    }
    file.close();
}

void BioInformation::read_subcellular() {
    fstream file(SUBCELLULAR);
    if(!file.is_open()) {
        cout << "--- [ERROR]: Failed to open file: " << SUBCELLULAR << " ---" << endl;
        exit(1);
    }
    string line;
    while(getline(file, line)) {
        istringstream iss(line);
        string str;
        vector<string> splitstring;
        while(iss >> str) {
            splitstring.emplace_back(str);
        }
        string protein = splitstring[0];
        string go = splitstring[2];
        double score = stod(splitstring[splitstring.size() - 1]);

        
    }
    file.close();
}

// 计算go_expression相似性, pearson correlation coefficient
double BioInformation::go_expression_pearson(const string& protein1, const string& protein2) {
    auto it1 = go_expression.find(protein1);
    auto it2 = go_expression.find(protein2);
    if(it1 == go_expression.end()) {
        // cout << "--- [ERROR]: Does not contain this protein in Gene Expression! --- " << protein1 << endl;
        return 0.0;
    }

    if(it2 == go_expression.end()) {
        // cout << "--- [ERROR]: Does not contain this protein in Gene Expression! --- " << protein2 << endl;
        return 0.0;
    }
    
    size_t n = it1->second.size();
    double sum_x = 0.0;
    double sum_y = 0.0;
    for(int i = 0; i < n; ++i) {
        sum_x += it1->second[i];
        sum_y += it2->second[i];
    }
    double aver_x = sum_x / n;
    double aver_y = sum_y / n;
    double cov = 0.0;
    double s1 = 0.0, s2 = 0.0;
    for(int i = 0; i < n; ++i) {
        cov += (it1->second[i] - aver_x) * (it2->second[i] - aver_y);        
        s1 += (it1->second[i] - aver_x) * (it1->second[i] - aver_x);        
        s2 += (it2->second[i] - aver_y) * (it2->second[i] - aver_y);        
    }
    double denominator = sqrt(s1) * sqrt(s2);
    if(denominator == 0.0) {
        cout << "--- [ERROR]: denominator is zero ---" << endl;
        exit(1);
    }
    return cov / denominator;
}

double BioInformation::go_expression_spearman(const string& protein1, const string& protein2) {
    auto it1 = go_expression.find(protein1);
    auto it2 = go_expression.find(protein2);
    if(it1 == go_expression.end()) {
        // cout << "--- [ERROR]: Does not contain this protein in Gene Expression! --- " << protein1 << endl;
        return 0.0;
    }

    if(it2 == go_expression.end()) {
        // cout << "--- [ERROR]: Does not contain this protein in Gene Expression! --- " << protein2 << endl;
        return 0.0;
    }
    double sum_x = 0.0, sum_y = 0.0;
    int n = it1->second.size();
    for(int i = 0; i < n; ++i) {
        sum_x += it1->second[i];
        sum_y += it2->second[2];
    }
    double average_x = sum_x / n;
    double average_y = sum_y / n;

    double molecular = 0.0;
    double denominator_x = 0.0, denominator_y = 0.0;
    for(int i = 0; i < n; ++i) {
        double xi = it1->second[i] - average_x;
        double yi = it2->second[i] - average_y;
        molecular += xi * yi;
        denominator_x += pow(xi, 2);
        denominator_y += pow(yi, 2);
    }

    return molecular / sqrt(denominator_x * denominator_y);
}

// 
double BioInformation::go_slim_slim(const string& protein1, const string& protein2) {
    auto it1 = go_slim.find(protein1); 
    auto it2 = go_slim.find(protein2);

    if(it1 == go_slim.end()) {
        // cout << "--- [WARNNING]: Does not contain this protein go_slim in GO SLIM ---" << protein1 << endl;
        return 0.0;
    }
    if(it2 == go_slim.end()) {
        // cout << "--- [WARNNING]: Does not contain this protein go_slim in GO SLIM --- " << protein2 << endl;   
        return 0.0;
    }

    set<string> common_go_slim;
    set_intersection(it1->second.begin(), it1->second.end(),
        it2->second.begin(), it2->second.end(), 
        inserter(common_go_slim, common_go_slim.begin()));

    return (double)pow(common_go_slim.size(), 2) / double(it1->second.size() * it2->second.size());
}

double BioInformation::go_slim_slim_jaccard(const string& protein1, const string& protein2) {
    auto it1 = go_slim.find(protein1);
    auto it2 = go_slim.find(protein2);
    if(it1 == go_slim.end() || it2 == go_slim.end()) {
        return 0.5;
    }

    set<string> common_go_slim;
    set<string> all_go_slim;
    set_union(it1->second.begin(), it1->second.end(),
        it2->second.begin(), it2->second.end(), 
        inserter(all_go_slim, all_go_slim.begin()));
    set_intersection(it1->second.begin(), it1->second.end(),
        it2->second.begin(), it2->second.end(), 
        inserter(common_go_slim, common_go_slim.begin()));
    
    if(all_go_slim.size() == 0) {
        return 0.0;
    }
    return (double)common_go_slim.size() / (double) all_go_slim.size();
}
