/*
 * @brief:
 * @Author: jh
 * @Date: 2024-05-06 16:00:21
 * @LastEditTime: 2024-05-28 17:19:56
 */
#include "../include/config.h"
#include "../include/dag.h"
#include "../include/gene_express.h"
#include "../include/ungraph.h"
#include <algorithm>
#include <queue>
#include <set>
#include <unistd.h>
#include <vector>

DAG UnGraph::ancestor_child(IS_A_A2C, PART_OF_A2C);
DAG UnGraph::child_ancestor(IS_A_C2A, PART_OF_C2A);
map<string, set<string>> UnGraph::protein_go = UnGraph::read_protein_go(GO_SLIM);
int main(int argc, char **argv) {
  int opt;
  string ppi_file;
  string out_file;

  while ((opt = getopt(argc, argv, "i:o:")) != -1) {
    string temp = optarg;
    switch (opt) {
      case 'i':
        ppi_file = "/home/jh/code/JHPC/dataset/Yeast/PPI/" + temp + ".txt";
        break;
      case 'o':
        out_file = "/home/jh/code/JHPC/result/" + temp + ".txt";
        break;
      default:
        std::cerr << "Usage: " << argv[0] << " -n <value> -s <value>"
            << std::endl;
        return 1;
    }
  }

  if (ppi_file.empty() || out_file.empty()) {
    std::cerr << "Error: -i and -n options require values." << std::endl;
    return 1;
  }
  cout << "ppi_file: " << ppi_file << endl;
  cout << "out_file: " << out_file << endl;

  UnGraph g(ppi_file);
  g.weight_by_go_term();
  // g.calculate_balanced_weight();

  GeneExpress gene_express(GENE_EXPRESSION);
  vector<UnGraph> dpins = gene_express.build_dynamic_PPI(&g, DPIN_MEHTOD::TOP);
  std::cout << "dpins.size " << dpins.size() << endl;

  std::cout << "splitting ppi ... \n";
  queue<UnGraph> queue_ppi;
  for (auto &dg: dpins) {
    queue_ppi.push(dg);
  }

  vector<UnGraph> splitted_ppi;
  while (!queue_ppi.empty()) {
    UnGraph::split_graph(queue_ppi, splitted_ppi);
  }
  std::cout << splitted_ppi.size() << std::endl;

  cout << "get_complex ...\n";
  vector<set<string> > complexes;
  for (auto &sg: splitted_ppi) {
    vector<set<string> > temp_complexes;
    UnGraph::get_complexes(g, sg, temp_complexes, 0.45);
    complexes.insert(complexes.end(), temp_complexes.begin(),
                     temp_complexes.end());
    std::cout << complexes.size() << endl;
  }
  sort(complexes.begin(), complexes.end(), Complex::CompareSetBySize);

  Complex::write_complex_to_file(complexes, out_file);

  return 0;
}
