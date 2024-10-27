from enum import Enum
from bisect import bisect_left
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, FuncFormatter
import seaborn as sns
import numpy as np
import random as rd

class SpannerType(Enum):
    FGVbase = "FB"
    MPVXbase = "MB"
    MPVXcompact = "MC"

iscloseint = lambda x: abs(x - round(x)) < 1e-5


class Spanners:
    class Spanner:
        def __init__(self, name, type, time, totalEdges, spannerEdges):
            self.name, self.type, self.time, self.totalEdges, self.spannerEdges = name, type, time, totalEdges, spannerEdges
        @property
        def spannerPercentage(self): return self.spannerEdges / self.totalEdges * 100
        def __str__(self): return f'{self.name} {self.type} {self.time}s {self.totalEdges} {self.spannerEdges}'

    def __init__(self):
        self.spanners = {}

    def addSpanner(self, spanner):
        if spanner.name not in self.spanners:
            self.spanners[spanner.name] = []
        self.spanners[spanner.name].append(spanner)
        self.spanners[spanner.name].sort(key=lambda s: s.type.value)

    def plot(self):
        def forward(x: float|np.ndarray) -> float|np.ndarray:
            if not isinstance(x, np.ndarray): x = np.array([x])
            x = np.where(x < fbmap[0,0], fbmap[0,0], x)
            ind = lambda x: bisect_left(fbmap, x, lo=1, hi=len(fbmap)-1, key=lambda x: x[0])
            index = np.array([ind(xi) for xi in x]) if isinstance(x, np.ndarray) else ind(x)
            f = lambda index, vals: np.where('log' == fbmaplinlog[index-1], np.log(vals), vals)
            return (f(index,x) - f(index,fbmap[index-1,0])) * (fbmap[index,1] - fbmap[index-1,1]) / (f(index,fbmap[index,0]) - f(index,fbmap[index-1,0])) + fbmap[index-1,1]

        
        def backward(x):
            if not isinstance(x, np.ndarray): x = np.array([x])
            x = np.where(x < 0, 0, x)
            ind = lambda x: bisect_left(fbmap, x, lo=1, hi=len(fbmap)-1, key=lambda x: x[1])
            index = np.array([ind(xi) for xi in x]) if isinstance(x, np.ndarray) else ind(x)
            fe = lambda index, vals: np.where('log' == fbmaplinlog[index-1], np.exp(vals), vals)
            fl = lambda index, vals: np.where('log' == fbmaplinlog[index-1], np.log(vals), vals)
            return fe(index,((x) - (fbmap[index-1,1])) * (fl(index,fbmap[index,0]) - fl(index,fbmap[index-1,0])) / (fbmap[index,1] - fbmap[index-1,1]) + fl(index,fbmap[index-1,0]))
        
        fig, [ax0,ax1] = plt.subplots(2,1)
        cmap = sns.color_palette("crest", len(self.spanners))
        x = 0
        for i, [name, spanners] in enumerate(sorted(self.spanners.items(), key=lambda s: s[1][0].totalEdges)):
            xTick = x + (len(spanners)-1) / 2
            color = cmap[0]
            cmap.remove(color)

            for j, spanner in enumerate(spanners):
                ax0.bar((x+j)*.8+xTick*.2, spanner.time, color=color, width=.6, zorder=2)
                ax0.text((x+j)*.8+xTick*.2, spanner.time, f'{spanner.type.value}\n{spanner.time:.0f}s', ha='center', va='bottom', color='black', fontsize=7)

                ax1.bar((x+j)*.8+xTick*.2, spanner.spannerPercentage, color=color, width=.6)
                ax1.text((x+j)*.8+xTick*.2, spanner.spannerPercentage, f'{spanner.type.value}\n{spanner.spannerPercentage:.0f}%', ha='center', va='bottom', color='black', fontsize=7)


            ax0.text(xTick, -1, f'{name}', ha='right', va='top', rotation=45)
            ax0.text(xTick+1, -1, f'm={spanners[0].totalEdges}', ha='right', va='top', rotation=45, fontsize=7)
            ax1.text(xTick, -.3, name, ha='right', va='top', rotation=45)
            x += len(spanners)


        fbmap = np.array([[0.0001,0],[10,20],[200,100]]) # map from values to percentages 
        fbmaplinlog = np.array(['log', 'lin', 'lin'])
        
        ax0.set_yscale('function', functions=(lambda x: np.array([forward(xi) for xi in x]), lambda x: np.array([backward(xi) for xi in x])))
        ax0.set_xticks([]) # remove all xticks
        ax1.set_xticks([]) # remove all xticks


        ax0.set_yticks(np.hstack((backward(np.linspace(fbmap[0,1], fbmap[1,1], 6)), np.linspace(fbmap[1,0],fbmap[2,0], 20)[1::2])))
        ax0.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f'{y:.0f}s' if iscloseint(y) else f'{y:.0e}s'))

        for i in range(-5,3): ax0.axhline(10**i, color='gray', linewidth=.5, linestyle='--', zorder=1, dashes=(2,5))

        ax0.set_ylabel('Time (s)')
        ax1.set_ylabel('Spanner Size (%)')
        ax0.set_ylim(0,330)
        ax1.set_ylim(0,53)
        ax1.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f'{y:.0f}%'))

        ax0.set_title("Spanners \u2014 Euler execution time")
        ax1.set_title("Spanners \u2014 Spanner size")
        plt.legend(handles=[ax0.plot([],[])[0] for _ in SpannerType], handlelength=0, handleheight=0, handletextpad=0, labels=[f'{stype.value}: {stype.name}' for stype in SpannerType], title='Spanner Types')
        plt.tight_layout()
        plt.subplots_adjust(hspace=.5)
        plt.show()




def parseLogTiming(logFile):
    import regex as re
    with open(logFile, 'r') as f:
        content = f.read()

        # Extract filename
        filename = re.search(r'(?<=### Graph: (/.*))([^/\n]*)(?=.(txt|edges)\n)', content).group(0)
        # spannertype = re.search(r'(?<=### Application: (.*))([^ \n]*)(?=\)\n)', content).group(0)
        time = re.search(r'(?<=# time per iter: ).*(?=\n)', content).group(0)
        totalEdges = re.search(r'### m: (\d+)', content).group(1)
        spannerEdges = re.search(r'(Total Edges |Spanner size = |total edges = )(\d+)', content).group(2)

        ## dammit spannertype is not correctly set
        spannertype = SpannerType.FGVbase if 'Value of k' in content else SpannerType.MPVXbase if 'Spanner size' in content else SpannerType.MPVXcompact


        return filename, spannertype, float(time), int(totalEdges), int(spannerEdges)


# get all files from the directory graphs-results-timing/slurmlog/
def getFiles(folder='graphs-results-timing/slurmlog/'):
    import os
    files = os.listdir(folder)
    files = [folder + file for file in files]
    return files


def main():
    files = getFiles('graphs-results-timing/slurmlog1010/') + getFiles('graphs-results-timing/slurmlog1024/')
    spanners = Spanners()
    for file in files:
        name, spannertype, time, totalEdges, spannerEdges = parseLogTiming(file)

        print (f'{name:20}\t{spannertype.name:10}\t{time:10}\t{spannerEdges/totalEdges:15}\t{spannerEdges:10}\t{totalEdges:10}')
        spanners.addSpanner(Spanners.Spanner(name, spannertype, time, totalEdges, spannerEdges))

    spanners.plot()
    

if __name__ == '__main__':
    main()