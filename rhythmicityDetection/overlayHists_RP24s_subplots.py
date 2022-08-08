# rmf 4.15.2019, last modified 5.19.2022

import sys, os
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plot
from matplotlib.patches import Patch
import numpy as np
from scipy import stats
import StringIO

# USAGE
usage = 'python ' + sys.argv[0] + ' <young transcript gene group file> <old transcript gene group file> <young mir gene group file> <old mir gene group file>'
if len(sys.argv) != 5 or '-h' in sys.argv or '--help' in sys.argv:
    print usage
    exit()

# SUBROUTINES
def readGroupFile(groupFile):
    RP24s = []
    with open(groupFile,'r') as groupFile:
        next(groupFile)  # skip header
        for line in groupFile:
            info = line.strip().split('\t')
            if info[2] != 'NOTEST':
                RP24 = float(info[2])
                RP24s.append(RP24)
    return RP24s

def plotHist(ax, transcripts, mirs, bin_min, bin_max, pval, title):
    binwidth = 0.5
    bins = np.arange(bin_min, bin_max)-0.5

    ax.hist(transcripts, bins, density = True, histtype='step', fill=True, color='lightpink', alpha =0.5, label = 'transcripts')
    ax.hist(mirs, bins, density = True, histtype='step', fill=True, color='palegreen', alpha =0.5, label = 'microRNAs')

    legend_elements = [Patch(facecolor = 'lightpink', alpha = 0.5),
                       Patch(facecolor = 'palegreen', alpha = 0.5),
                       Patch(facecolor = 'white')]
    legend_labels = ['transcripts', 'microRNAs', pval.getvalue()]
    ax.legend(legend_elements, legend_labels)

    ax.set_title(title)
    ax.set_xlabel('RP24 Values')
    ax.set_ylabel('Density')

# ARGUMENTS and MAIN
yTranscriptFile = sys.argv[1]
oTranscriptFile = sys.argv[2]
yMirFile = sys.argv[3]
oMirFile = sys.argv[4]

yTranscripts = readGroupFile(yTranscriptFile)
oTranscripts = readGroupFile(oTranscriptFile)
yMirs = readGroupFile(yMirFile)
oMirs = readGroupFile(oMirFile)

# perform KS test 
y_pval = stats.ks_2samp(yTranscripts, yMirs)[1]
print('young transcripts vs mirs p-value= '+str(y_pval)+' (KS test)')
o_pval = stats.ks_2samp(oTranscripts, oMirs)[1]
print('old transcripts vs mirs p-value= '+str(o_pval)+' (KS test)')

# format pvalue for legend label 
y_pvalLabel = StringIO.StringIO()
y_pvalLabel.write("KS p-value = %.3g" % y_pval)
o_pvalLabel = StringIO.StringIO()
o_pvalLabel.write("KS p-value = %.3g" % o_pval)

fig, axs = plot.subplots(2,1, sharex = True)

this_min = min([min(yTranscripts), min(oTranscripts), min(yMirs), min(oMirs)])
this_max = max([max(yTranscripts), max(oTranscripts), max(yMirs), max(oMirs)])

plotHist(axs[0], yTranscripts, yMirs, this_min, this_max, y_pvalLabel, 'young')
plotHist(axs[1], oTranscripts, oMirs, this_min, this_max, o_pvalLabel, 'old')

plot.tight_layout()
plot.savefig('histRP24s_mirsVsTranscripts.pdf')
