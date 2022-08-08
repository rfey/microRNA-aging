#!/bin/python

# rmf 10.26.2018, original script from Dave, last modified 11.4.2021

import sys, os, re
import numpy as np
import scipy.stats
import matplotlib as mpl
import StringIO
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# USAGE
usage = "python " + sys.argv[0] + " <fileList adjCountsTables> <fileList totalCountsTables>"
if len(sys.argv) != 3:
    print usage
    sys.exit()

# SUBROUTINES
def readFileList_adjCounts(fileListAdjCounts):
    totalExpr = {}
    with open(fileListAdjCounts,'r') as f:
        for line in f:
            acFile = line.strip() # remove newline
            readCountsTable(acFile,totalExpr)
    return totalExpr

def readCountsTable(countsTable,totalExpr):
    pattern = re.compile('adjCountsTable_(.*?).txt')
    match = pattern.search(countsTable)
    sample = match.group(1)
    print sample
    # initializing dict keys
    ZTs = ["ZT0","ZT4","ZT8","ZT12","ZT16","ZT20"]
    totalExpr[sample] = {}  # 2D dict
    for ZT in ZTs:
        totalExpr[sample][ZT] = []
    # initializing counts variables as float
    counts0 = 0.0
    counts4 = 0.0
    counts8 = 0.0
    counts12 = 0.0
    counts16 = 0.0
    counts20 = 0.0
    with open(countsTable,'r') as countsTable:
        for line in countsTable:
            next(countsTable) # skip header line
            mir, ZT0, ZT4, ZT8, ZT12, ZT16, ZT20 = line.strip().split('\t')
            counts0 += float(ZT0)
            counts4 += float(ZT4)
            counts8 += float(ZT8)
            counts12 += float(ZT12)
            counts16 += float(ZT16)
            counts20 += float(ZT20)
    # append summed expression of all mirs to appropriate time point
    totalExpr[sample]["ZT0"] = counts0
    totalExpr[sample]["ZT4"] = counts4
    totalExpr[sample]["ZT8"] = counts8
    totalExpr[sample]["ZT12"] = counts12
    totalExpr[sample]["ZT16"] = counts16
    totalExpr[sample]["ZT20"] = counts20
    return totalExpr

def readFileListTotalCounts(fileListTotalCounts):
    totals = {}
    with open(fileListTotalCounts,'r') as f:
        for line in f:
            tcFile = line.strip()
            totals = readTotalCountsTable(tcFile,totals)
    return totals

def readTotalCountsTable(totalCountsTable,totals):
    pattern = re.compile('totalCountsTable_(.*?\d).*') # totalCountsTable_oldRep1_AS.txt
    match = pattern.search(totalCountsTable)
    sample = match.group(1)
    print sample
    totals[sample] = {}  # 2D dict
    with open(totalCountsTable,'r') as f:
        for line in f:
            timepoint, total = line.strip().split('\t')
            tp = timepoint.strip() # remove extra space
            ZT = "ZT" + tp
            totals[sample][ZT] = float(total)
    return totals

def calculateRatio(totalExpr,totals):
    ratios = {}
    for sample in totals:
        ratios[sample] = {}
        for ZT in totals[sample]:
            total = totals[sample][ZT]
            mirCounts = totalExpr[sample][ZT]
            ratio = mirCounts/total  # reads mapping to miRs / total reads mapped
            ratios[sample][ZT] = ratio
    return ratios

# ARGUMENTS and MAIN
fileListAdjCounts = sys.argv[1]
fileListTotalCounts = sys.argv[2]

expr = readFileList_adjCounts(fileListAdjCounts)
totals = readFileListTotalCounts(fileListTotalCounts)
ratios = calculateRatio(expr,totals)

outfile = 'globalTrendsWithAge_boxplot.pdf'

# wrangle data
young,old = [],[]
for sample in ratios:
    if 'young' in sample:
        for ZT in ratios[sample]:
            young.append(ratios[sample][ZT])
    elif 'old' in sample:
        for ZT in ratios[sample]:
            old.append(ratios[sample][ZT])

# perform significance testing
pval = scipy.stats.ttest_ind(young, old)[1]
print pval


exprData = [np.asarray(young),np.asarray(old)]

# format pvalue for legend label
pvalLabel = StringIO.StringIO()
pvalLabel.write("p-value = %.3g" % pval)

# create a figure instance
plt.figure()
ax = plt.gca()

# create the boxplot
bp = plt.boxplot(exprData, labels = ['young','old'], patch_artist=True, showfliers=False, widths=.6)
plt.ylabel('microRNA-mapped reads /\ntotal genome-mapped reads',size=14)
plt.xticks([],[])
plt.tick_params(labelsize=14) # sets both y axis and x axis label fontsize


# fill with colors
colors = ['red','blue']
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_edgecolor('black')
    patch.set_alpha(0.5)
for median in bp['medians']:
    median.set(color='k')

# make legend
legend_elements = [Patch(facecolor = "red", edgecolor = "red", label = "young", alpha = 0.5),
                   Patch(facecolor = "blue", edgecolor = "blue", label = "old", alpha = 0.5),
                   Patch(facecolor = "white", edgecolor = "white", label = pvalLabel.getvalue())]
ax.legend(handles = legend_elements, loc = "upper right")

y1 = exprData[0]
x1 = np.random.normal(1, 0.02, len(y1))
plt.plot(x1, y1, 'k',linestyle='None',marker='.')
y2 = exprData[1]
x2 = np.random.normal(2, 0.02, len(y2))
plt.plot(x2, y2, 'k',linestyle='None',marker='.')

# save plot
plt.savefig(outfile, bbox_inches='tight')
