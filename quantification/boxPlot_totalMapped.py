#!/bin/python

# rmf 5.17.2019, last updated 11.4.2021

import sys, os, re
import numpy as np
import scipy.stats
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Patch
import matplotlib.ticker as mticker
import StringIO

# USAGE
usage = "python " + sys.argv[0] + " <fileList totalCountsTables>"
if len(sys.argv) != 2:
    print usage
    sys.exit()

# SUBROUTINES
def readFileListTotalCounts(fileListTotalCounts):
    totals = {}
    with open(fileListTotalCounts,'r') as f:
        for line in f:
            tcFile = line.strip()
            totals = readTotalCountsTable(tcFile,totals)
    return totals

def readTotalCountsTable(totalCountsTable,totals):
    pattern = re.compile('totalCountsTable_(.*?).txt') # totalCountsTable_oldRep1_AS.txt
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

# ARGUMENTS and MAIN
fileListTotalCounts = sys.argv[1]
print "Calculating total reads per sample:"

totals = readFileListTotalCounts(fileListTotalCounts)

outfile = 'boxplot_totalReadsMapped_youngVsOld.pdf'

# wrangle data
young,old = [],[]
for sample in totals:
    if 'young' in sample:
        for ZT in totals[sample]:
            young.append(totals[sample][ZT])
    elif 'old' in sample:
        for ZT in totals[sample]:
            old.append(totals[sample][ZT])

# calculate significance
pval = scipy.stats.ttest_ind(young, old)[1]
print pval

exprData = [np.asarray(young),np.asarray(old)]

# format pvalue for legend label
pvalLabel = StringIO.StringIO()
pvalLabel.write("p-value = %.3g" % pval)

# create a figure instance
plt.figure()
ax = plt.gca()
ax.xaxis.set_major_formatter(ticker.ScalarFormatter(useOffset=False))
f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
sciNotation = lambda x,pos : "${}$".format(f._formatSciNotation('%1.10e' % x))
plt.gca().yaxis.set_major_formatter(mticker.FuncFormatter(sciNotation))

# create the boxplot
bp = plt.boxplot(exprData, labels = ['young','old'], patch_artist=True, showfliers=False, widths=.6)
plt.ylabel('genome-mapped reads',size=14)
plt.tick_params(labelsize=14) # sets both y axis and x axis label fontsize
plt.xticks([],[])

# fill with colors
colors = ['red','blue']
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_edgecolor('black')
    patch.set_alpha(0.5)
for median in bp['medians']:
    median.set(color='k')

legend_elements = [Patch(facecolor = "red", edgecolor = "red", label = "young", alpha = 0.5),
                   Patch(facecolor = "blue", edgecolor = "blue", label = "old", alpha = 0.5),
                   Patch(facecolor = "white", edgecolor = "white", label = pvalLabel.getvalue())]

ax.legend(handles = legend_elements, loc = "upper left")

y1 = exprData[0]
x1 = np.random.normal(1, 0.02, len(y1))
plt.plot(x1, y1, 'k',linestyle='None',marker='.')
y2 = exprData[1]
x2 = np.random.normal(2, 0.02, len(y2))
plt.plot(x2, y2, 'k',linestyle='None',marker='.')

# save plot
plt.savefig(outfile, bbox_inches='tight')
