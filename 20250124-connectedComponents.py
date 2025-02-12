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
        
    # print (len(sets[0]))
    # print (sets[0])
    print (f"Non-trivial components {len(sets)}")
    print (f"Components with 1 node and no edge {82168 - sum([len(s) for s in sets])}")
    return len(sets) == 1

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
    # edges = loadGraph('graphs-sanitized-snap/Slashdot.edges', removeFirst=1)
    # edges = loadGraph('graphs-results/graphs-results-1010/Slashdot_spanner_MPVXbase.txt')
    edges = loadGraph('graphs-removed-spanner/Slashdot-removed-MPVXbase-1.edges', removeFirst=1)
    # edges = loadGraph('graphs-results/graphs-results-1010/soc-hamsterster_spanner_MPVXbase.txt')
    print(isConnected(edges))

    # bridges = [(1264, 24306), (1264, 24307), (781, 20771), (781, 20772), (781, 20774), (1436, 25681), (432, 17893), (1200, 23909), (411, 17551), (17590, 66943), (411, 17590), (1703, 27599), (1703, 27607), (18959, 67978), (550, 18966), (550, 18967), (675, 19951), (675, 19955), (675, 19959), (2880, 35464), (35468, 75533), (2880, 35468), (35473, 75535), (35473, 75536), (2880, 35473), (2880, 35478), (2880, 35479), (2880, 35484), (2880, 35497), (2880, 35500), (2880, 35501), (2880, 35503), (2880, 35512), (629, 19588), (71602, 81247), (1335, 24882), (1335, 24888), (1467, 25964), (1467, 25967), (1467, 25975), (1467, 25976), (1467, 25983), (309, 15965), (309, 15966), (309, 15969), (309, 15977), (309, 15981), (309, 15987), (309, 15988)]

    # for bridge in bridges:
    #     copy = edges.copy()
    #     copy.remove(bridge)
    #     print(isConnected(copy))                              