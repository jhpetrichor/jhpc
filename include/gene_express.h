/*
 * @Description: 
 * @Author: jh
 * @Date: 2024-04-30 17:11:48
 * @LastEditTime: 2024-06-12 07:57:23
 */

#ifndef __GENE_EXPRESS_H__
#define __GENE_EXPRESS_H__

// #include "ungraph.h"
#include "config.h"
#include "graph.h"

#include <iostream>
#include <map>
#include <utility>
#include <vector>
#include <unordered_map>

using namespace std;

enum DPIN_MEHTOD{
    THREE_SIGMA,   // three-sigma
    TOP,           // 前。。。
    TIME,          // 时间周期
    ESSENTIAL
};

enum Level {
    HIGH,
    MEDIUM,
    LOW,
};

class GeneExpress {
public:
    map<string, vector<double>> gene_express;
    set<string> essential_proteins;

    unordered_map<string, vector<Level>> expression_level;
    unordered_map<string, pair<double, double>> level_probability;

    explicit GeneExpress(string file_path = GENE_EXPRESSION);
    ~GeneExpress();
    void read_gene_express(string& file_path);
    void read_essential_protein(string& file_path);
    map<string, double> activate_by_three_sigma();
    map<string, double> active_by_top();
    // vector<UnGraph> build_KPIN(const UnGraph* g);
    vector<Graph> build_KPIN(const Graph* g);
    vector<Graph> build_KPIN1(const Graph* g);
    // vector<UnGraph> build_KPIN2(const UnGraph* g);
    // vector<UnGraph> build_KPIN_min(const UnGraph* g);
    // vector<UnGraph> build_KPINN(const UnGraph* g);

    // vector<UnGraph> build_dynamic_PPI(const UnGraph* g, DPIN_MEHTOD method = DPIN_MEHTOD::THREE_SIGMA);


    // 相关性
    double mutual_information(string& p1, string& p2);

private:
    void calculate_expression_level();

    /**
     * calculate gene active threshold by three-sigma principle
     * @param mean 基因表达的均值
     * @param varience 基因表达的方差
     * @return 基因的活性阈值
     */
    static double active_three_sigma(double mean, double varience);
    /**
     * @description: 通过基因表达活性阈值构建动态网络，构建的动态网络数量和基因表达采样点一致
     * @param {UnGraph*} g 原始PPI网络
     * @param {int} count 构建的动态网路数量
     * @return {*} 动态网络（原始网络的子网络）
     */    
    // vector<UnGraph> build_dynamic_PPI_by_active(const UnGraph* g, map<string, double>& active, int count);
    // /**
    //  * @description: 通过时间周期构建动态网络
    //  * @param {UnGraph*} g 原始PPI网络
    //  * @param {int} count 构建的动态网路数量
    //  * @return {*} 动态网络（原始网络的子网络）
    //  */    
    // vector<UnGraph> build_dynamic_PPI_by_time(const UnGraph* g, map<string, double>& active, int count);
    // vector<UnGraph> build_dynamic_PPI_by_essential_protein(const UnGraph* g, map<string, double>& active, int count);
    /**
     * @brief: 为g构建基于关键蛋白质的动态网络
     * @param {UnGraph*} g
     * @return 动态网络的结果
     */    
    pair<double, double> 
    static inline calculate_mean_varience(vector<double>& arr);
};


#endif
