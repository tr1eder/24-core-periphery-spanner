from typing import List, Tuple
from collections import defaultdict

import sys
sys.setrecursionlimit(10**6)  # Set to a large enough value

Edge = Tuple[int, int]

def cutEdges(edges: List[Edge]) -> List[Edge]:
    from collections import defaultdict
    
    # Convert edge list to adjacency list
    graph = defaultdict(list)
    for u, v in edges:
        graph[u].append(v)
        graph[v].append(u)

    # Variables for DFS
    timer = 0
    tin = {}  # Discovery time of each node
    low = {}  # Lowest discovery time reachable
    bridges = []
    visited = set()
    
    def dfs(u: int, parent: int):
        nonlocal timer
        visited.add(u)
        tin[u] = low[u] = timer
        timer += 1
        
        for v in graph[u]:
            if v == parent:  # Ignore the edge to the parent
                continue
            if v not in visited:  # Forward edge
                dfs(v, u)
                low[u] = min(low[u], low[v])  # Update low-link value
                
                if low[v] > tin[u]:  # Bridge condition
                    bridges.append((u, v))
            else:  # Back edge
                low[u] = min(low[u], tin[v])
    
    # Run DFS for each component
    for node in graph:
        if node not in visited:
            dfs(node, -1)
    
    return bridges

def loadGraph(filename: str, removeFirst=0) -> List[Edge]:
    edges = []
    with open(filename, 'r') as f:
        for i, line in enumerate(f):
            if line.startswith('#') or i < removeFirst:
                continue
            edge = tuple(map(int, line.strip().split()))
            edges.append(edge)

        return edges

if __name__ == '__main__':
    edges = loadGraph('graphs-sanitized-snap/Slashdot.edges', removeFirst=1)
    # edges = loadGraph('graphs-results/graphs-results-1010/Slashdot_spanner_MPVXbase.txt')
    # edges = loadGraph('graphs-removed-spanner/Slashdot-removed-MPVXbase-1.edges', removeFirst=1)

    bridges = cutEdges(edges)
    print (bridges)
    print (len(bridges))