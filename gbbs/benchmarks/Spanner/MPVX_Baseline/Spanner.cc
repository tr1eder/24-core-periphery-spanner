//MIT License

// Copyright (c) 2022 Maulein
// This code is part of our work titled “Scalable algorithms for compact spanners on real-world graphs"

// This work is adapted from the work done under project "Theoretically Efficient Parallel Graph
// Algorithms Can Be Fast and Scalable", presented at Symposium on Parallelism
// in Algorithms and Architectures, 2018. 
// GBBS: Graph Based Benchmark Suite - https://github.com/ParAlg/gbbs
// Copyright (c) 2018 Laxman Dhulipala, Guy Blelloch, and Julian Shu

// This code is part of the project "Theoretically Efficient Parallel Graph
// Algorithms Can Be Fast and Scalable", presented at Symposium on Parallelism
// in Algorithms and Architectures, 2018.
// Copyright (c) 2018 Laxman Dhulipala, Guy Blelloch, and Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all  copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// Usage:
// numactl -i all ./Spanner -rounds 3 -s -m twitter_SJ
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -rounds : the number of times to run the algorithm
//     -stats : print the #ccs, and the #vertices in the largest cc

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <filesystem>
#include "Spanner.h"

namespace gbbs {
// Beta should be set to log n/2k. See Corollary 3.1 and Lemma 3.2 in MPVX'15.
template <class Graph>
double Spanner_runner(Graph& G, commandLine P) {
  size_t n = G.n;
  size_t k = P.getOptionLongValue("-k", 4);
  double beta = log(n) / (2 * k);
  std::cout << "### Application: Spanner (O(k)-spanner from MPXV)" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### Params: -k = " << k << " => \\beta = \\log n/2k = " << beta
            << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  assert(P.getOption("-s"));
  timer t;
  t.start();
  auto spanner = spanner::Spanner(G, beta);
  double tt = t.stop();
  std::cout << "### Running Time: " << tt << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  
  // ! ADDED
  std::string filename = std::string(P.getArgument(0));
  std::filesystem::path filepath(filename);

  std::string graphname = filepath.stem().string();
  std::string pathname  = filepath.parent_path().parent_path().string();
  // create the resulting spanner at the subfolder graphs-results
  std::string output_filename = pathname + "/graphs-results/" + graphname + "_spanner_MPVXbase.txt";

  std::ofstream output_file(output_filename);
  if (!output_file.is_open()) {
    std::cerr << "Error: could not open file " << output_filename << std::endl;
    return -1;
  }

  // write the number of nodes and edges to the output file
  output_file << "# Spanner" << std::endl;
  output_file << "# Nodes: " << G.n << " Edges: " << spanner.size() << std::endl;
  output_file << "# FromNodeId	ToNodeId" << std::endl;

  // first write spanner edges into a string stream -- maybe split into multiple batches for large graphs
  std::stringstream ss;
  for (auto edge : spanner) {
    ss << edge.first << " " << edge.second << std::endl;
  }
  output_file << ss.str();

  // close file
  output_file.close();
  // ! END ADDED




  return tt;
}
}  // namespace gbbs

generate_main(gbbs::Spanner_runner, false);