//MIT License

// Copyright (c) 2022 Maulein
// This code is part of our work titled “Scalable algorithms for compact spanners on real-world graphs"

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


#pragma once

#include "benchmarks/LowDiameterDecomposition/MPX13/LowDiameterDecomposition.h"
#include "gbbs/gbbs.h"
#include "gbbs/graph.h"
#include "gbbs/graph_io.h"
#include "gbbs/helpers/dyn_arr.h"
#include "gbbs/helpers/sparse_table.h"
#include "math.h"
namespace gbbs {
namespace spanner {

using edge = std::pair<uintE, uintE>;

struct cluster_and_parent {
  uintE cluster;
  uintE parent;
  cluster_and_parent(uintE _cluster, uintE _parent)
      : cluster(_cluster), parent(_parent) {}
};

sequence<size_t> generate_shifts_geomcap(size_t n, size_t k) {
  // Create k levels
  uintE last_round = k-1;
  uintE r = k-1;
  double p = (1 - (1/pow(n,(double)1/k))) ; 
 
  auto shifts = sequence<size_t>(last_round + 2);
  parallel_for(2, last_round+2, kDefaultGranularity, [&](size_t i) { 
 	shifts[i] = ceil(n*p*pow((1-p),last_round - (i-1)));
   });
  shifts[0] = 0;
  shifts[1] = ceil(n*pow((1-p),r));

  for(uintE j = 1;j<last_round+2;j++)
   shifts[j] += shifts[j-1];
  return shifts;
}


template <class Graph>
dyn_arr<uintE> HighCov_for_centers(Graph& G,uintE k_param)
{
	size_t n = G.n;
	auto twoHop = sequence<uintE>(n + 1);
	auto reduced_twoHop = sequence<uintE>(n + 1);
	auto reduced_degree = sequence<uintE>(n + 1);
	auto reduced_twoHop_2 = sequence<uintE>(n + 1);

	parallel_for(0, n, 1, [&](size_t i) 
	{
  		uintE deg_u = G.get_vertex(i).out_degree();	
  		reduced_degree[i] = deg_u;
		twoHop[i]=0;

		if(k_param==2)
			reduced_twoHop[i]=reduced_degree[i];
		else
		{
      			auto v_iter = G.get_vertex(i).out_neighbors().get_iter();

			uintE ct_u = 0;
			while (ct_u < deg_u) 
			{
    			  	auto vertex = std::get<0>(v_iter.cur());
          		  	if (v_iter.has_next()) v_iter.next();
				ct_u++;
    				twoHop[i] += G.get_vertex(vertex).out_degree();
			}
			reduced_twoHop[i]=twoHop[i];
		}
	});

  auto inMIS = sequence<uintE>(n + 1);
  auto active = sequence<uintE>(n + 1);
  auto frontier = sequence<uintE>(n + 1);
  auto propagate = sequence<uintE>(n+1);
  auto parent = sequence<uintE>(n + 1);
  auto propagate_next = sequence<uintE>(n+1);
  auto parent_next = sequence<uintE>(n + 1);
  auto more_next = sequence<uintE>(n + 1);
  auto more = sequence<uintE>(n + 1);
  uintE frontierSize=n+1, old_frontierSize =-1;

  parallel_for(0, n, 1, [&](size_t i) {
    inMIS[i] = 0;
    active[i] = 1;
    frontier[i] = 1;
    propagate[i] = twoHop[i];
    propagate_next[i] = twoHop[i];
    });
  
  uintE remove_hops = k_param-2;  //till remove_hops + 1 distanced neighbours removed
  auto highcov_centers = gbbs::dyn_arr<uintE>(20);
  uintE iter = 0;
  uintE alpha=10; 
  while(frontierSize!=0 && iter < alpha)
  {
		iter++;
                parallel_for(0, n, 1, [&](size_t i) {
			if(frontier[i] && reduced_twoHop[i] == 0)
				active[i] = 0;
			if(active[i])
				inMIS[i] = 1;
		
			more[i] = 0;
			parent[i] = parent_next[i] = 0;
			propagate[i]=reduced_twoHop[i];
		});
  		parallel_for(0, n, 1, [&](size_t i) {
		if(frontier[i])
		{
  			uintE deg_u = G.get_vertex(i).out_degree();	
      			auto v_iter = G.get_vertex(i).out_neighbors().get_iter();
			uintE ct_u = 0;
			while (ct_u < deg_u) {
    			   auto vertex = std::get<0>(v_iter.cur());
          		   if (v_iter.has_next()) v_iter.next();
			        ct_u++;          		  
			   if(frontier[vertex])
			   {
				if((propagate[i]<reduced_twoHop[vertex]) ||((propagate[i] == reduced_twoHop[vertex]) && (i>vertex )))
				{
					inMIS[i] = 0;
					propagate[i]= reduced_twoHop[vertex];		
					parent[i]= vertex;
					more[i] = remove_hops;
				}
			   }	
  			}
 		}
		});
		uintE k =remove_hops;
		while(k>=1)
		{
  			parallel_for(0, n, 1, [&](size_t i) {
  				uintE deg_u = G.get_vertex(i).out_degree();	
				uintE propagate_final = propagate[i];
	 			uintE parent_final = parent[i];
	 			uintE more_final = more[i];
				uintE ct_u = 0;
				while (ct_u < deg_u) {
      					auto v_iter = G.get_vertex(i).out_neighbors().get_iter();
    				  	auto vertex = std::get<0>(v_iter.cur());
          			 	if (v_iter.has_next()) v_iter.next();
					ct_u++;
				  	if((more[vertex] == k) && (parent[vertex] !=i))
					{
					
						if((propagate_final<propagate[vertex]) || ((propagate_final == propagate[vertex]) && parent_final>parent[vertex]))
						{
							inMIS[i] = 0;			
							propagate_final=propagate[vertex];
							parent_final = parent[vertex];
							more_final = k-1;
						}
					}
				}
				propagate_next[i] = propagate_final;
				parent_next[i] = parent_final;
				more_next[i] = more_final;
			});
			k--;
			
  			parallel_for(0, n, 1, [&](size_t i) {
				propagate[i] = propagate_next[i];
				parent[i] = parent_next[i];
				more[i] = more_next[i];
			});
		}
                parallel_for(0, n, 1, [&](size_t i) {
			more[i] = 0;
		});

  		parallel_for(0, n, 1, [&](size_t i) {
		if(frontier[i])
		{
  			uintE deg_u = G.get_vertex(i).out_degree();
			uintE max =0;	
      			auto v_iter = G.get_vertex(i).out_neighbors().get_iter();
			uintE ct_u = 0;
			while (ct_u < deg_u) {
    			   auto vertex = std::get<0>(v_iter.cur());
          		   if (v_iter.has_next()) v_iter.next();
			     ct_u++;          		  
			   if(frontier[vertex] && inMIS[vertex] && max<reduced_twoHop[vertex])
			   {
				active[i] = 0;
				inMIS[i] = 0;
				more[i] = remove_hops;
			   }	
  			}
 		}
		});
		k =remove_hops;
		while(k>=1)
		{

  			parallel_for(0, n, 1, [&](size_t i) {
				{
					if(frontier[i])
					{
  						uintE deg_u = G.get_vertex(i).out_degree();	
		      				auto v_iter = G.get_vertex(i).out_neighbors().get_iter();
						uintE ct_u = 0;
						while (ct_u < deg_u) {
    						   auto vertex = std::get<0>(v_iter.cur());
         			 		   if (v_iter.has_next()) v_iter.next();
						   ct_u++;          		  
					  	   if(more[vertex] == k)
						   {
							if(parent[vertex] != i)
							{
								inMIS[i] = 0;			
								active[i] = 0;
								more[i] = k-1;
							}
						  }  
			
						}
					}
				}	
			});
			k--;
		}

      auto candidates_f = [&](size_t i) {
     			return static_cast<uintE>(i);
      };

      auto candidates = parlay::delayed_seq<uintE>(n+1, candidates_f);
      auto pred = [&](uintE v) { return (frontier[v] && inMIS[v]); };
      auto centers = parlay::filter(candidates, pred);
      

      highcov_centers.copyIn(centers, centers.size());
                parallel_for(0, n, 1, [&](size_t i) {
			frontier[i] = 0;
		});
      int flag=0;
      if(frontierSize == old_frontierSize)
		flag=1;
      old_frontierSize = frontierSize;
      frontierSize=0;
		
      parallel_for(0, n, 1, [&](size_t i) {
		if(active[i] && !inMIS[i])
				frontier[i] = 1;
     });
     frontierSize =
    		  parlay::reduce(frontier, parlay::addm<size_t>());

     parallel_for(0, n, 1, [&](size_t i) {
		uintE deg_u = G.get_vertex(i).out_degree();	
      		auto v_iter = G.get_vertex(i).out_neighbors().get_iter();
		uintE ct_u = 0;
		while (ct_u < deg_u) {
    			   auto vertex = std::get<0>(v_iter.cur());
        		   if (v_iter.has_next()) v_iter.next();
				   ct_u++;          		  
			   if(frontier[vertex]==1 && !active[vertex])
			   	reduced_degree[i]=reduced_degree[i] - 1;
		}
     });
		
     parallel_for(0, n, 1, [&](size_t i) {
		if(k_param==2)
			reduced_twoHop[i]=reduced_degree[i];
		else
		{
			reduced_twoHop[i]=0;
			uintE deg_u = G.get_vertex(i).out_degree();	
	      		auto v_iter = G.get_vertex(i).out_neighbors().get_iter();
			uintE ct_u = 0;
			while (ct_u < deg_u) {
				   auto vertex = std::get<0>(v_iter.cur());
       				   if (v_iter.has_next()) v_iter.next();
					   ct_u++;
					   if(frontier[vertex])          		  
					   	reduced_twoHop[i]+=reduced_degree[vertex];
			}
		}
	});
	if(flag==1)
		break;	
  }
  return highcov_centers;
}

template <class Graph, class C>
sequence<edge> fetch_intercluster_te(Graph& G, C& clusters,
                                     size_t num_clusters, sequence<uintE> &diam,  sequence<uintE>& level,sequence<uintE>& clusterType,size_t k_param) {
  using W = typename Graph::weight_type;
  gbbs_debug(std::cout << "Running fetch edges te" << std::endl;);
  using K = edge;
  using V = edge;
  using KV = std::tuple<K, V>; 
  size_t n = G.n;
  gbbs_debug(std::cout << "num_clusters = " << num_clusters << std::endl;);
  timer count_t;
  count_t.start();
  auto deg_map = sequence<uintE>(n + 1);
  auto pred = [&](const uintE& src, const uintE& ngh, const W& w) {
    uintE c_src = clusters[src];
    uintE c_ngh = clusters[ngh];
    return c_src < c_ngh;
  };
  parallel_for(0, n, 1, [&](size_t i) {
    deg_map[i] = G.get_vertex(i).out_neighbors().count(pred);
  });
  deg_map[n] = 0;
  parlay::scan_inplace(deg_map);
  count_t.stop();
  gbbs_debug(count_t.next("count time"););

  timer ins_t;
  ins_t.start();
  auto empty_edge = std::make_pair(UINT_E_MAX, UINT_E_MAX);
  KV empty = std::make_tuple(empty_edge, empty_edge);
  auto hash_pair = [](const edge& t) {
    size_t l = std::min(t.first, t.second);
    size_t r = std::max(t.first, t.second);
    size_t key = (l << 32) + r;
    return parlay::hash64_2(key);
  };
  auto edge_table = gbbs::make_sparse_table<K, V>(deg_map[n], empty, hash_pair);
  gbbs_debug(std::cout << "sizeof table = " << edge_table.m << std::endl;);
  deg_map.clear();

  auto map_f = [&](const uintE& src, const uintE& ngh, const W& w) {
    uintE c_src = clusters[src];
    uintE c_ngh = clusters[ngh];
    uintE l_src = level[src];
    uintE l_ngh = level[ngh];
    if(c_src!=c_ngh)
    {
    	if((diam[c_src] + diam[c_ngh]) <=(2*(k_param-1)))
   	 {

    		if (c_src < c_ngh) {
    	  		edge_table.insert(std::make_tuple(std::make_pair(c_src, c_ngh),
                                        std::make_pair(src, ngh)));
   		}

  	 } 	
  	 else
  	 {
            if(clusterType[c_src] == 0 && clusterType[c_ngh] == 0)
            {
        	if ((l_src > l_ngh) || ((l_src == l_ngh) && (c_ngh < c_src ) )) {
                	  edge_table.insert(std::make_tuple(std::make_pair(src, c_ngh),
                                        std::make_pair(src, ngh)));
        	 } 
    	    }
            else if(clusterType[c_src] == 1 || clusterType[c_ngh] == 0)
	    {
                	  edge_table.insert(std::make_tuple(std::make_pair(src, c_ngh),
                                        std::make_pair(src, ngh)));
	    }
            else if(clusterType[c_src] == 0 || clusterType[c_ngh] == 1)
	    {                    
		//Edge from Baseline cluster to Predetermined cluster - not added
	    }
            else
	    {
    		if (c_src < c_ngh) {
                	  edge_table.insert(std::make_tuple(std::make_pair(ngh, c_src),
                                        std::make_pair(ngh, src)));
	        }
	    }
  	  }
    }
  };
  parallel_for(0, n, 1,
               [&](size_t i) { G.get_vertex(i).out_neighbors().map(map_f); });
  auto edge_pairs = edge_table.entries();
  ins_t.stop();
  gbbs_debug(ins_t.next("ins time"););
  gbbs_debug(std::cout << "edges.size = " << edge_pairs.size() << std::endl);

  auto edges = sequence<edge>::from_function(
      edge_pairs.size(), [&](size_t i) { return std::get<1>(edge_pairs[i]); });
  return edges;
}

template <class Graph, class C>
sequence<edge> fetch_intercluster(Graph& G, C& clusters, size_t num_clusters) {
  using K = edge;
  using V = edge;
  using KV = std::tuple<K, V>;
  using W = typename Graph::weight_type;

  size_t n = G.n;
  gbbs_debug(std::cout << "num_clusters = " << num_clusters << std::endl;);
  size_t estimated_edges = num_clusters * 5;

  timer ins_t;
  ins_t.start();
  auto empty_edge = std::make_pair(UINT_E_MAX, UINT_E_MAX);
  KV empty = std::make_tuple(empty_edge, empty_edge);
  auto hash_pair = [](const edge& t) {
    size_t l = std::min(t.first, t.second);
    size_t r = std::max(t.first, t.second);
    size_t key = (l << 32) + r;
    return parlay::hash64_2(key);
  };

  auto edge_table =
      gbbs::make_sparse_table<K, V>(estimated_edges, empty, hash_pair);
  gbbs_debug(std::cout << "sizeof table = " << edge_table.m << std::endl;);

  bool abort = false;
  auto map_f = [&](const uintE& src, const uintE& ngh, const W& w) {
    uintE c_src = clusters[src];
    uintE c_ngh = clusters[ngh];
    if (c_src < c_ngh) {
      edge_table.insert_check(std::make_tuple(std::make_pair(c_src, c_ngh),
                                              std::make_pair(src, ngh)),
                              &abort);
    }
  };
  parallel_for(
      0, n, [&](size_t i) { G.get_vertex(i).out_neighbors().map(map_f); }, 1);
  if (abort) {
    gbbs_debug(std::cout << "calling fetch_intercluster_te" << std::endl;);
    return fetch_intercluster_te(G, clusters, num_clusters);
  }
  auto edge_pairs = edge_table.entries();
  ins_t.stop();
  gbbs_debug(ins_t.next("ins time"););
  gbbs_debug(std::cout << "edges.size = " << edge_pairs.size() << std::endl);

  auto edges = sequence<edge>::from_function(
      edge_pairs.size(), [&](size_t i) { return std::get<1>(edge_pairs[i]); });
  return edges;
}

template <class Graph>
sequence<edge> tree_and_intercluster_edges(
    Graph& G, sequence<cluster_and_parent>& cluster_and_parents, sequence<uintE> &diam,sequence<uintE>& level, sequence<uintE>& clusterType, size_t k_param) {
  size_t n = G.n;
  auto edge_list = gbbs::dyn_arr<edge>(2 * n);

  // Compute and add in tree edges.
  auto tree_edges_with_loops = parlay::delayed_seq<edge>(n, [&](size_t i) {
    return std::make_pair(i, cluster_and_parents[i].parent);
  });
  auto tree_edges = parlay::filter(tree_edges_with_loops, [&](const edge& e) {
    return e.first != e.second;
  });
  edge_list.copyIn(tree_edges, tree_edges.size());

  // Compute inter-cluster using hashing.
  auto clusters = parlay::delayed_seq<uintE>(
      n, [&](size_t i) { return cluster_and_parents[i].cluster; });
  sequence<bool> flags(n, false);
  parallel_for(0, n, [&](size_t i) {
    uintE cluster = clusters[i];
    if (!flags[cluster]) {
      flags[cluster] = true;
    }
  });
  auto cluster_size_seq = parlay::delayed_seq<size_t>(
      n, [&](size_t i) { return static_cast<size_t>(flags[i]); });
  size_t num_clusters =
      parlay::reduce(cluster_size_seq, parlay::addm<size_t>());

  auto intercluster = fetch_intercluster_te(G, clusters, num_clusters, diam,level,clusterType, k_param);
  gbbs_debug(std::cout << "num_intercluster edges = " << intercluster.size()
                  << std::endl;);
  edge_list.copyIn(intercluster, intercluster.size());
  size_t edge_list_size = edge_list.size;
  std::cout << "total edges = " << edge_list.size<< std::endl;
  return sequence<edge>::from_function(
      edge_list_size, [&](size_t i) { return edge_list.A[i]; });
}

template <class W>
struct LDD_Parents_F {
  cluster_and_parent* clusters;

  LDD_Parents_F(cluster_and_parent* _clusters) : clusters(_clusters) {}

  inline bool update(const uintE& s, const uintE& d, const W& wgh) {
    clusters[d].cluster = clusters[s].cluster;
    clusters[d].parent = s;
    return true;
  }

  inline bool updateAtomic(const uintE& s, const uintE& d, const W& wgh) {
    if (gbbs::atomic_compare_and_swap(&clusters[d].cluster, UINT_E_MAX,
                                      clusters[s].cluster)) {
      clusters[d].parent = s;
      return true;
    }
    return false;
  }

  inline bool cond(uintE d) { return clusters[d].cluster == UINT_E_MAX; }
};

template <class Graph>
inline sequence<cluster_and_parent> LDD_parents(Graph& G, double beta, sequence<uintE> &diam ,sequence<uintE>& level,sequence<uintE>& clusterType, 
                                                size_t k_param, bool permute = true) {
  using W = typename Graph::weight_type;
  size_t n = G.n;
  sequence<uintE> vertex_perm;
  if (permute) {
    vertex_perm = parlay::random_permutation<uintE>(n);
  }
  auto shifts = generate_shifts_geomcap(n, k_param);
  auto clusters = sequence<cluster_and_parent>(
      n, cluster_and_parent(UINT_E_MAX, UINT_E_MAX));

  size_t round = 0, num_visited = 0;
  vertexSubset frontier(n);  // Initially empty
  size_t num_added = 0;
  uintE round_highcov=0;
  uintE change =0;
  while (num_visited < n) {
    sequence<uintE> centers;
    size_t num_to_add = 0;
   if(change)	
    {
    size_t start = shifts[round];
    size_t end = std::min(static_cast<size_t>(shifts[round + 1]), n);
    num_to_add = end - start;
    if (num_added + num_to_add > n+1)
	num_to_add = n+1 -num_added;
    }
    else
    if(round_highcov ==0)
    {	
  	auto centers1 = HighCov_for_centers(G,k_param);
	num_to_add= (size_t)centers1.size;
        centers = sequence<uintE>::from_function(
        centers1.size, [&](size_t i) { 
                clusterType[centers1.A[i]]=1;
		return centers1.A[i]; });
    }
    else
	num_to_add = 0;
   
    if (num_to_add > 0) {
      	assert((num_added + num_to_add) <= n);
	if(num_added+num_to_add>n)    //Ensure for loop runs till i=n-1 only
		num_to_add = n+1-num_added;

	if (num_to_add <= 0)
		break;
      	auto candidates_f = [&](size_t i) {
        if(change)
	{	
		if (permute)
        		return vertex_perm[num_added + i];
     	  	else
     			return static_cast<uintE>(num_added + i);
	}
        else if (round_highcov == 0)
			return centers[i];
      };

      auto candidates = parlay::delayed_seq<uintE>(num_to_add, candidates_f);
      auto pred = [&](uintE v) { return clusters[v].cluster == UINT_E_MAX; };
      auto new_centers = parlay::filter(candidates, pred);
      add_to_vsubset(frontier, new_centers.begin(), new_centers.size());
      parallel_for(0, new_centers.size(), kDefaultGranularity, [&](size_t i) {
        uintE v = new_centers[i];
        clusters[new_centers[i]] = cluster_and_parent(v, v);

        level[new_centers[i]]=round;   //ADDED
      });
	if(round_highcov ==0)
		num_to_add =0;

      num_added += num_to_add;
    }

    num_visited += frontier.size();
    if (num_visited >= n) 
		break;
   if(change==0 && round_highcov == k_param-1)
		frontier = std::move(0);
    auto ldd_f = LDD_Parents_F<W>(clusters.begin());
    vertexSubset next_frontier =
       		edgeMap(G, frontier, ldd_f, -1, sparse_blocked);

    auto flag = sequence<uintE>(n + 1);
    parallel_for(0, n, 1, [&](size_t i) {
      flag[i] = 0;
    });

    vertexMap(next_frontier, [&](const uintE u) {
		if(change == 0)
			level[u] = round_highcov + 1;
		else
	 		level[u] = round+1;
		
 	 	// Select the parent with least ID
		uintE deg_u = G.get_vertex(u).out_degree();	
		uintE ct_u = 0;
		while (ct_u < deg_u) {
      			auto v_iter = G.get_vertex(u).out_neighbors().get_iter();
    			auto vertex = std::get<0>(v_iter.cur());
       	 		if (v_iter.has_next()) v_iter.next();
				ct_u++;
			if(level[vertex] == round-1 && vertex<clusters[u].parent)
			{
				clusters[u].cluster = clusters[vertex].cluster; 
				clusters[u].parent = vertex;
			}
		}
    });

   auto f = std::move(frontier);

   auto find_cluster = [&](uintE v) { flag[clusters[v].cluster] = 1; };
   vertexMap(f, find_cluster);

   parallel_for(0, n, 1, [&](size_t i) {
    	if(flag[i])
		diam[i] = diam[i] + 1;
   });


   frontier = std::move(next_frontier);

   if(change == 1)	
    	round++;
   else
   if(round_highcov==(k_param))
	change =1;
   else
      	round_highcov ++ ;
  }

   parallel_for(0, n, 1, [&](size_t i) {
		diam[i] *= 2;
   });
  return clusters;
}

template <class Graph>
inline sequence<edge> Spanner_impl(Graph& G, double beta, size_t k_param) {
  bool permute = true;


  timer ldd_t;
  ldd_t.start();
  uintE n=G.n;
  auto diam = sequence<uintE>(n + 1);
  auto level = sequence<uintE>(G.n + 1);
  auto clusterType = sequence<uintE>(G.n + 1); // 0: Baseline cluster, 1: Predetermined cluster
  parallel_for(0, n+1, 1, [&](size_t i) {
    diam[i] = -1;
    clusterType[i]=0;
    level[i] = UINT_E_MAX;
    });
auto clusters_and_parents = LDD_parents(G, beta, diam, level,clusterType, k_param,permute);
  ldd_t.stop();
  gbbs_debug(ldd_t.next("ldd time"););
  timer build_el_t;
  build_el_t.start();
  auto spanner_edges = tree_and_intercluster_edges(G, clusters_and_parents, diam, level, clusterType, k_param);
  build_el_t.stop();
  gbbs_debug(build_el_t.next("build spanner edges time"););

  // return spanner as an edge-list.
  gbbs_debug(std::cout << "Spanner size = " << spanner_edges.size() << std::endl;);
  return spanner_edges;
}

template <class Graph>
inline sequence<edge> Spanner(Graph& G, double beta, size_t k_param) {
  return Spanner_impl(G, beta, k_param);
}

}  // namespace cc
}  // namespace gbbs
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
