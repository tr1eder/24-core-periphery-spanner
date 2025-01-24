from typing import List, Tuple


Edge = Tuple[int, int]

def find(sets: List[set], v: int) -> int:
    for i in range(len(sets)):
        if v in sets[i]:
            return i
    return -1


def isConnected(edges: List[Edge]) -> bool:
    sets = []
    for edge in edges:
        v1 = edge[0]
        v2 = edge[1]
        f1 = find(sets, v1)
        f2 = find(sets, v2)
        if f1 == -1 and f2 == -1:
            sets.append({v1, v2})
        elif f1 == -1 or f2 == -1:
            if f1 == -1:
                sets[f2].add(v1)
            else:
                sets[f1].add(v2)
        else:
            if f1 != f2:
                sets[f1] = sets[f1].union(sets[f2])
                sets.pop(f2)
        
    return len(sets) == 1

def loadGraph(filename: str) -> List[Edge]:
    edges = []
    with open(filename, 'r') as f:
        for line in f:
            if (line.startswith('#')):
                continue
            edge = tuple(map(int, line.strip().split(' ')))
            edges.append(edge)

        return edges
    

if __name__ == '__main__':
    edges = loadGraph('graphs-results/graphs-results-1010/Slashdot_spanner_MPVXbase.txt')
    # edges = loadGraph('graphs-results/graphs-results-1010/soc-hamsterster_spanner_MPVXbase.txt')
    print(isConnected(edges))