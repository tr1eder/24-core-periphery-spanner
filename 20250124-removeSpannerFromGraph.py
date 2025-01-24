import re
from typing import List, Tuple, Literal

Edge = Tuple[int, int]

def removeSpannerFromGraph(graph: List[Edge], spanner: List[Edge]) -> List[Edge]:
    spannerSet = set(spanner)
    graphSet = set(graph)

    count = 0
    for edge in spanner:
        if edge not in graphSet and (edge[1], edge[0]) not in graphSet:
            print(f'Edge {edge} not in graph!')
        else: 
            count += 1
    print (f'Edges in spanner and in graph: {count}')


    return [edge for edge in graph if edge not in spannerSet and (edge[1], edge[0]) not in spannerSet]


def loadGraph(filename: str, removeFirst = 0) -> List[Edge]:
    edges = []
    with open(filename, 'r') as f:
        for i, line in enumerate(f):
            if line.startswith('#') or i < removeFirst:
                continue
            edge = tuple(map(int, line.strip().split()))
            edges.append(edge)

        return edges

def storeGraph(graph: List[Edge], nodes: int, filename: str, typ: Literal['snap']):
    with open(filename, 'w') as f:
        f.write(f'{nodes} {len(graph)}\n')
        for edge in graph:
            f.write(f'{edge[0]} {edge[1]}\n')
    
def getNodeEdgeCount(filename: str, typ: Literal['snap', 'spanner']) -> Tuple[int, int]:
    with open(filename, 'r') as f:
        for i, line in enumerate(f):
            if typ == 'spanner':
                match = re.search(r'# Nodes: (\d+) Edges: (\d+)', line)
                if match:
                    return int(match.group(1)), int(match.group(2))
            if typ == 'snap' and i == 0:
                return map(int, line.strip().split())
            
        
        return -1, -1
    
def sanitize(graph: List[Edge]) -> List[Edge]:
    ## remove self loops
    ## remove duplicates (a, b) and (b, a)
    graphSet = set(graph)
    keep = lambda x,y: (y,x) in graphSet and x < y or (y,x) not in graphSet and x != y

    return [edge for edge in graph if keep(edge[0], edge[1])]


if __name__ == '__main__':
    graphname = 'Slashdot'
    spannertype = 'MPVXbase'
    edgesSpanner = loadGraph(f'graphs-results/graphs-results-1010/{graphname}_spanner_{spannertype}.txt')
    edgesGraph = loadGraph(f'graphs-sanitized-snap/{graphname}.edges', removeFirst=1)
    edgesGraphSanitized = sanitize(edgesGraph)
    edgesSpannerSanitized = sanitize(edgesSpanner)
    gNodes, gEdges = getNodeEdgeCount(f'graphs-sanitized-snap/{graphname}.edges', 'snap')
    sNodes, sEdges = getNodeEdgeCount(f'graphs-results/graphs-results-1010/{graphname}_spanner_{spannertype}.txt', 'spanner')

    edgesGraphWithoutSpanner = removeSpannerFromGraph(edgesGraphSanitized, edgesSpanner)

    print(f'Graph: {gNodes} nodes, {gEdges} edges')
    print(f'Graph sanitized: {len(edgesGraphSanitized)} edges')
    print(f'Spanner: {sNodes} nodes, {sEdges} edges')
    print(f'Spanner sanitized: {len(edgesSpannerSanitized)} edges')
    print(f'Graph without spanner: {len(edgesGraphWithoutSpanner)} edges')
    print(f'Should be: {len(edgesGraphSanitized) - sEdges} edges')

    if len(edgesGraphWithoutSpanner) != len(edgesGraphSanitized) - sEdges:
        print('ERROR: Graph without spanner is not the expected size')
    else: 
        print('SUCCESS: Graph without spanner is the expected size')
        print('Writing to file...')
        storeGraph(edgesGraphWithoutSpanner, gNodes, f'graphs-removed-spanner/{graphname}-removed-{spannertype}-1.edges', 'snap')
        print('Done!')


    