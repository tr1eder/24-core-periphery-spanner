import csv
from enum import Enum
from bisect import bisect_left
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, FuncFormatter
import seaborn as sns
import numpy as np
import random as rd

class GraphType(Enum):
    Original = "OG"
    FGVbase = "FB"
    MPVXbase = "MB"
    MPVXcompact = "MC"

previous = lambda x: GraphType.Original if x == GraphType.Original or x == GraphType.FGVbase else GraphType.FGVbase if x == GraphType.MPVXbase else GraphType.MPVXbase
iscloseint = lambda x: abs(x - round(x)) < 1e-5
AVG = lambda lst, default: np.average(lst) if lst else default
MAX = lambda lst, default: np.max(lst) if lst else default
MIN = lambda lst, default: np.min(lst) if lst else default

class Accuracies:
    class Accuracy:
        def __init__(self, name, edges):
            self.name = name
            self.edges = edges
            self.dists = {}


        def avgAdd(self, graphtype):
            return AVG([dist-orig for dist, orig in zip(self.dists[graphtype], self.dists[GraphType.Original]) if dist != -1 and orig != -1], 0)
        def avgMul(self, graphtype):
            return AVG([dist/orig for dist, orig in zip(self.dists[graphtype], self.dists[GraphType.Original]) if dist != -1 and orig != -1], 1)
        def avg(self, graphtype, distOrig):
            return AVG([dist for dist, orig in zip(self.dists[graphtype], self.dists[GraphType.Original]) if orig == distOrig and dist != -1], distOrig)
        def maxAvg(self, graphtype):
            return MAX([self.avg(graphtype, dist) for dist in range(self.minVal(GraphType.Original), self.maxVal(GraphType.Original)+1)], -1)
        def maxVal(self, graphtype):
            return MAX([dist for dist in self.dists[graphtype]], -1)
        def minVal(self, graphtype):
            return MIN([dist for dist in self.dists[graphtype] if dist != -1], -1)
        def failed(self, graphtype):
            return sum([dist == -1 for dist in self.dists[graphtype]])
        def failedRel(self, graphtype):
            return self.failed(graphtype)/len(self.dists[graphtype])
        def invalid(self, graphtype): # number of invalid distances, where the distance is below the original distance
            return sum([dist < orig for dist, orig in zip(self.dists[graphtype], self.dists[GraphType.Original]) if dist != -1 and orig != -1])
        def failedStr(self, graphtype):
            return f'{self.failed(graphtype)}/{len(self.dists[graphtype])}'


        
        def __str__(self): return f'{self.name}'
        def __eq__(self, other): return other.__class__ == self.__class__ and self.name == other.name

    def __init__(self):
        self.graphs = []

    def add(self, name, graphtype, dist, edges):
        if all([name != graph.name for graph in self.graphs]):
            self.graphs.append(self.Accuracy(name, edges))

        graph = next(graph for graph in self.graphs if graph.name == name)
        graph.dists[graphtype] = dist

    def plot(self):

        maxval = max([graph.maxVal(graphtype) for graph in self.graphs for graphtype in graph.dists])
        maxvalorig = max([graph.maxVal(GraphType.Original) for graph in self.graphs])

        fig, ax = plt.subplots()
        # cmap = sns.color_palette("crest", 2)

        for i, graph in enumerate(sorted(self.graphs, key=lambda g: g.edges)):
            failed = graph.failed(GraphType.Original)
            gmin = graph.minVal(GraphType.Original)
            gmax = graph.maxVal(GraphType.Original)

            cmap = sns.color_palette("crest", gmax+1-gmin)

            # if graph.name == "amazon":
            #     print ([f"{graph.avg(GraphType.MPVXcompact, dist):.2f}" for dist in range(gmin, gmax+1)])

            x = lambda i, graphtype: i + (graphtype == GraphType.Original) * -.3 + (graphtype == GraphType.FGVbase) * -.1 + (graphtype == GraphType.MPVXbase) * .1 + (graphtype == GraphType.MPVXcompact) * .3
            y = lambda o, graphtype: graph.avg(graphtype, o) if o >= gmin else 0

            for graphtype in graph.dists.keys():


                for distOrig in range(gmin, gmax+1):
                    pos = x(i, graphtype)
                    top = y(distOrig, graphtype)
                    bottom = y(distOrig-1, graphtype)

                    if graphtype != GraphType.Original:
                        ax.plot([x(i,previous(graphtype)), x(i, graphtype)], [y(distOrig,previous(graphtype)), y(distOrig,graphtype)], color='black', linewidth=.5, zorder=3)


                    if top < bottom: continue
                    color = (*cmap[(distOrig-gmin)][:3], (1-graph.failedRel(graphtype)))
                    ax.bar(pos, top-bottom, bottom=bottom, color=color, width=.15, zorder=2)

                failedText = f"\nF{graph.failed(graphtype)}" if graph.failed(graphtype)!=failed else ''
                invalidText = f"\nI{graph.invalid(graphtype)}" if graph.invalid(graphtype) > 0 else ''
                ax.text(pos, graph.maxAvg(graphtype), f'{graphtype.value}{failedText}{invalidText}', ha='center', va='bottom', color='black', fontsize=7)

                

            ax.text(i, -1, f'{graph.name}', ha='right', va='top', rotation=45)
            ax.text(i+.2, -1, f'm={graph.edges}', ha='right', va='top', rotation=45, fontsize=7)
            if failed > 0: ax.text(i+.38, -1, f'Failed: {graph.failedStr(graphtype)}', ha='right', va='top', rotation=45, fontsize=7)


        
        ax.set_xticks([]) # remove all xticks

        for i in range(0, maxvalorig+1): 
            ax.axhline(i, color='gray', linewidth=.5, linestyle='--', zorder=1, dashes=(2,5))

        ax.set_ylabel('Accuracy')
        ax.set_title("Spanners \u2014 Accuracy")
        plt.legend(handles=[ax.plot([],[])[0] for _ in GraphType], handlelength=0, handleheight=0, handletextpad=0, labels=[f'{stype.value}: {stype.name}' for stype in GraphType], title='Spanner Types')
        plt.tight_layout()
        plt.subplots_adjust(hspace=.5)
        plt.show()




def parseAccuracy(file):
    import regex as re
    import os
    with open(file, 'r') as f:

        filename = os.path.basename(file)
        graphname = re.search(r'(0_10000_)([^_]*)', filename).group(2)
        graphtype = GraphType.FGVbase if 'FGVbase' in filename else GraphType.MPVXbase if 'MPVXbase' in filename else GraphType.MPVXcompact if 'MPVXcompact' in filename else GraphType.Original
        distances = [int(line.strip()) for line in f if line.strip()]

        return graphname, graphtype, distances


# get all files from the directory
def getFiles(folder):
    import os
    files = os.listdir(folder)
    files = [folder + file for file in files if 'distances' in file]
    return files


def main():
    files = (getFiles('graphs-results-accuracy/dists-1027-3/'))
    edgecounts = {row[0]: int(row[1]) for row in csv.reader(open('graphs_edges.csv', 'r'))}

        #    + getFiles('graphs-results-timing/slurmlog1024/')
        #    + getFiles('graphs-results-timing/slurmlog1027/'))
    
    accuracies = Accuracies()
    for file in files:
        graphname, graphtype, distances = parseAccuracy(file)
        accuracies.add(graphname, graphtype, distances, edgecounts[graphname])

        print (f'{graphname:20}\t{graphtype.name:10}')

    accuracies.plot()
    

if __name__ == '__main__':
    main()