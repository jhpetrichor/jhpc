#ifndef __CONFIG_H__
#define __CONFIG_H__
// Yeast
// dag.h
#define PATH_IS_A        "/home/jh/code/JHPC/dataset/Yeast/DAG/is_a.txt"
#define PATH_PART_OF     "/home/jh/code/JHPC/dataset/Yeast/DAG/part_of.txt"
#define PATH_GO_TERMS    "/home/jh/code/JHPC/dataset/Yeast/DAG/go_term.txt"

// go_information
#define GO_SLIM          "/home/jh/code/JHPC/dataset/Yeast/DAG/protein-go.txt"
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

#define MAX_MATCHED        4
#define COMPLEX_MAX_SIZE   20

// p-value


#endif  // __CONFIG_H__
