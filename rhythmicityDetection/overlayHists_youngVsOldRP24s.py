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
usage = 'python ' + sys.argv[0] + ' <young geneGroups file> <old geneGroups file> <outfile base> <binwidth (try 0.1 for genes, 0.3 for miRs)>'
if len(sys.argv) != 5 or '-h' in sys.argv or '--help' in sys.argv:
    print usage
    sys.exit()

# SUBROUTINES
def readGroupFile(groupFile):
    RP24s = []
    with open(groupFile,'r') as groupFile:
        next(groupFile)  # skip header
        for line in groupFile:
            info = line.strip().split('\t')
            RP24 = info[2]
            if RP24 != 'NOTEST':
                RP24s.append(float(RP24))
    print len(RP24s)
    return RP24s


# ARGUMENTS and MAIN
fileYoung = sys.argv[1]
fileOld = sys.argv[2]
outbase = sys.argv[3]
binwidth = float(sys.argv[4])

outfile = 'histRP24s_' + outbase + '.pdf'

youngRP24s = readGroupFile(fileYoung)
oldRP24s = readGroupFile(fileOld)

# perform KS test
pval = stats.ks_2samp(youngRP24s,oldRP24s)[1]
print('old vs young p-value= '+str(pval)+' (KS test)')

# format pvalue for legend label
pvalLabel = StringIO.StringIO()
pvalLabel.write("KS p-value = %.3g" % pval)

this_min = min(min(youngRP24s),min(oldRP24s))
this_max = max(max(youngRP24s),max(oldRP24s))
bins = np.arange(this_min, this_max + binwidth, binwidth)

plot.hist(youngRP24s, bins = bins, color = 'red', alpha = 0.5)
plot.hist(oldRP24s, bins = bins, color = 'blue', alpha = 0.5)

plot.xlabel('RP24 Values')
plot.ylabel('Count')

legend_elements = [Patch(facecolor = 'red', alpha = 0.5),
                   Patch(facecolor = 'blue', alpha = 0.5),
                   Patch(facecolor = 'white')]
legend_labels = ['young', 'old', pvalLabel.getvalue()]
plot.legend(legend_elements, legend_labels)

plot.savefig(outfile)


