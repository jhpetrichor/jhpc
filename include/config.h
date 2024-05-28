/*
 * @brief: 
 * @Author: jh
 * @Date: 2024-05-06 13:24:02
 * @LastEditTime: 2024-05-28 14:49:39
 */
#ifndef __CONFIG_H__
#define __CONFIG_H__
#define ESSENTIAL_PROTEIN "/home/jh/code/JHPC/dataset/Yeast/essential_protein.txt"
// Yeast
// dag.h
#define IS_A_A2C        "/home/jh/code/JHPC/dataset/Yeast/DAG/is_a.txt"
#define PART_OF_A2C     "/home/jh/code/JHPC/dataset/Yeast/DAG/part_of.txt"
#define PATH_GO_TERMS    "/home/jh/code/JHPC/dataset/Yeast/DAG/go_term.txt"

#define IS_A_C2A        "/home/jh/code/JHPC/dataset/Yeast/DAG/is_a_child_ancestor.txt"
#define PART_OF_C2A     "/home/jh/code/JHPC/dataset/Yeast/DAG/part_of_child_ancestor.txt"
// #define PROTEIN_GO_FILE  "/home/jh/code/JHPC/dataset/Yeast/DAG/protein-go.txt"

// go_information
#define GO_SLIM          "/home/jh/code/JHPC/dataset/Yeast/DAG/go-slim.txt"
#define SUBCELLULAR      "/home/jh/code/JHPC/dataset/Yeast/DAG/subcellular.txt"
#define GENE_EXPRESSION  "/home/jh/code/JHPC/dataset/Yeast/DAG/gene-expression.txt"

// collins
#define COLLINS_PPI         "/home/jh/code/JHPC/dataset/Yeast/PPI/collins.txt"
#define COLLINS_RESULT      "/home/jh/code/JHPC/temp/collins"
// gavin
#define GAVIN_PPI         "/home/jh/code/JHPC/dataset/Yeast/PPI/gavin.txt"
#define GAVIN_RESULT      "/home/jh/code/JHPC/temp/gavin"
// krogan
#define KROGAN_CORE_PPI         "/home/jh/code/JHPC/dataset/Yeast/PPI/krogan_core.txt"
#define KROGAN_CORE_RESULT      "/home/jh/code/JHPC/temp/krogan_core"
// krogan
#define KROGAN_EXTENDED_PPI         "/home/jh/code/JHPC/dataset/Yeast/PPI/krogan_extended.txt"
#define KROGAN_EXTENDED_RESULT      "/home/jh/code/JHPC/temp/krogan_extended"
// biogrid
#define BIOGRID_PPI         "/home/jh/code/JHPC/dataset/Yeast/PPI/biogrid.txt"
#define BIOGRID_RESULT     "/home/jh/code/JHPC/temp/biogrid"

// 标准蛋白质复合物
#define COMPLEX_FILE     "/home/jh/code/JHPC/dataset/Yeast/complex/yeast_complex.txt"

// BOPS算法参数
#define BALANCED_INDEX    1.5
#define MATCH_INDEX       0.2
#define MAX_MATCH_INDEX   0.45     // 允许最大的匹配阈值
// Human

#define COMPLEX_MAX_SIZE   20

// p-value


#endif  // __CONFIG_H__
