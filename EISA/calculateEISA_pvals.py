# rmf 2.28.2019, last modified 5.20.2022

import sys

usage = 'python ' + sys.argv[0] + ' <EISA counts table> <BH qval threshold>'

if len(sys.argv) != 3 or '--help' in sys.argv or '-h' in sys.argv:
    print usage
    exit()

print 'Importing modules...'

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plot
import numpy as np
from scipy.stats import norm
from scipy.optimize import curve_fit
import statsmodels.stats.multitest as mt
from math import sqrt
from scipy import asarray as ar,exp

# SUBROUTINES
def readTableOfCounts(tableOfCounts):
    valsScatter = {}
    allInfo = {}
    symbols, FBIDs, scores = ([] for i in range(3))
    with open(tableOfCounts) as table:
        print 'Reading table...'
        next(table) # skip header row
        for line in table:
            FBID, symbol, exonY, exonO, exonFC, intronY, intronO, intronFC, EISA = line.strip().split('\t')
            valsScatter[symbol] = (float(exonFC),float(intronFC))
            symbols.append(symbol)
            FBIDs.append(FBID)
            scores.append(float(EISA))
            allInfo[FBID] = (symbol, exonY, exonO, exonFC, intronY, intronO, intronFC, EISA)
        valsHist = zip(symbols,FBIDs,scores)
    return valsScatter, valsHist, allInfo

def gaus(x,a,mu,sigma):
    return a*exp(-(x-mu)**2/(2*sigma**2))

def plotHist(valsHist,outbase,gaus,alpha):
    # plot hist
    binWidth = 0.1
    bins = np.arange(-6.0,9.0,binWidth) # hard coded for now so all plots are comparable
    symbols, FBIDs, scores = zip(*valsHist)
    print 'Plotting histogram...'
    counts,bins,bars = plot.hist(scores,bins,color='orange',alpha=0.5,edgecolor='orange',linewidth=1,density=True)
    # fit curve
    xvals = np.asarray(scores) # EISA scores
    mean = np.mean(xvals)
    sigma = sqrt(sum((xvals-mean)**2)/len(xvals)) # standard deviation
    p0 = [max(counts),mean,sigma] # for fitting 1 gaussian
    popt,pcov = curve_fit(gaus,np.delete(bins,len(bins)-1),counts,p0=p0)
    # calculate pvals
    pVals_raw = calculatePvals(valsHist,mean,sigma) # feed in zipped list
    pVals_adj = correctPvals(pVals_raw,alpha)
    # plot curve
    pg1 = popt[0:3]
    g1 = gaus(bins, *pg1)
    plot.plot(bins,g1,'r',label='Gaussian')
    # configure plot
    plot.title('Distribution of EISA Scores ' + outbase,fontsize=18)
    plot.xlabel("EISA Score\n$log_2 (exon_{old}/exon_{young}) - log_2 (intron_{old}/intron_{young})$",fontsize=14)
    plot.ylabel("Density",fontsize=14)
    plot.tight_layout()
    # save plot
    plot.savefig('histEISA_' + outbase + '.pdf')
    # clear plot memory in preparation for second plot
    plot.clf()
    return pVals_adj

def calculatePvals(vals,mean,sigma):
    pvals = []
    symbols, FBIDs, scores = zip(*vals)
    for symbol in symbols:
        i = symbols.index(symbol)
        EISA = scores[i]
        pval = 2*min(norm.sf(EISA,loc=mean,scale=sigma),norm.cdf(EISA,loc=mean,scale=sigma))
        pvals.append(pval)
    rawPvals = zip(symbols,FBIDs,scores,pvals)
    return rawPvals

def correctPvals(pvals_raw,alpha):
    symbols, FBIDs, scores, pvals = zip(*pvals_raw)
    pvals_adj = mt.multipletests(pvals,alpha=alpha,method='fdr_bh')[1] # return second array in tuple
    adjPvals = zip(symbols,FBIDs,scores,pvals,pvals_adj)
    return adjPvals

def getSigPvals(adjPvals,alpha,outbase,allInfo):
    print 'Retrieving significant p-values...'
    sig = {}
    count = 0
    symbols,FBIDs,scores,pvals_raw,pvals_adj = zip(*adjPvals)
    pvals_raw = list(pvals_raw) # convert from numpy array
    pvals_adj = list(pvals_adj)
    scores = list(scores)
    writePvalsFile(adjPvals,outbase, allInfo) # pass in zipped lists
    for symbol in symbols:  # loop through these to get index bc they are unique, pvals are not necessarily
        i = symbols.index(symbol)
        if pvals_adj[i] <= alpha:
            count += 1
            sig[symbol] = (pvals_adj[i],scores[i])
    return sig

def writePvalsFile(adjPvals,outbase,allInfo):
    outfile = 'pvalsEISA_' + outbase + '.txt'
    with open(outfile,'w') as outfile:
        outfile.write('symbol\tFBID\texonY\texonO\texonFC\tintronY\tintronO\tintronFC\tEISA_Score\tpval_raw\tpval_BH\n')
        symbols,FBIDs,scores,pvals,pvals_adj = zip(*adjPvals)
        for symbol in symbols:
            i = symbols.index(symbol)
            FBID = FBIDs[i]
            pval = pvals[i]
            pval_adj = pvals_adj[i]
            score = scores[i]
            exonY = allInfo[FBID][1]
            exonO = allInfo[FBID][2]
            exonFC = allInfo[FBID][3]
            intronY = allInfo[FBID][4]
            intronO = allInfo[FBID][5]
            intronFC = allInfo[FBID][6]
            outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%f\t%f\t%f\n' % (symbol,FBID,exonY,exonO,exonFC,intronY,intronO,intronFC,score,pval,pval_adj))
    outfile.close()

def plotScatter(valsScatter,sig,outbase,alpha):
    x = []
    y = []
    labels = []
    colors = []
    pval = []
    fitx = []
    fity = []
    # arrange data for plotting
    for symbol in valsScatter:
        exonFC = valsScatter[symbol][0]
        intronFC = valsScatter[symbol][1]
        x.append(intronFC) ## intron on x axis (control)
        y.append(exonFC)
        labels.append(symbol)
        if symbol not in sig:
            fitx.append(intronFC)
            fity.append(exonFC)
            colors.append('k')
            pval.append(1.0) # used for annotation purposes only
        elif symbol in sig:
            colors.append('r')
            pval.append(sig[symbol][0])
    # calculate linear fit
    print 'Calculating linear fit...'
    xLinFit = np.array(fitx) # fit only non-significant values
    m, b = np.polyfit(xLinFit, np.array(fity),1) # calculate linear fit model
    yLinFit = xLinFit*m+b  # y = mx+b
    # plot
    print 'Plotting scatter plot...'
    plot.scatter(x,y, s=1**1, color=colors)
    # annotate
    info = zip(x,y,labels,pval)
    info.sort(key = lambda x: x[3]) # sort on pval
    toPlot = info[0:10] # get top 10 significant symbols and info
    for item in toPlot:
        xval,yval,symbol,pval=item
        print symbol
        plot.annotate(symbol,(xval,yval))
    # plot linear fit line
    plot.plot(xLinFit,yLinFit,c='r',label='linear fit: y='+str(round(m,6))+'x + '+str(round(b,6)))
    # configure plot
    plot.title(outbase,fontsize=24)
    plot.xlabel("$Log_2 (Intron_{old}/Intron_{young})$",fontsize=20)
    plot.ylabel("$Log_2 (Exon_{old}/Exon_{young})$",fontsize=20)
    xLeft, xRight = plot.xlim()
    plot.xlim(xLeft-1,xRight+1)
    yBottom, yTop = plot.ylim()
    plot.ylim(yBottom-1,yTop+1)
    plot.legend()
    plot.tight_layout()
    # save plot as png bc all the points take up too much space as pdf
    plot.savefig("scatterEISA_" + outbase + ".png")
    # clear plot memory in preparation for second plot
    plot.clf()

# ARGUMENTS and MAIN
tableOfCounts = sys.argv[1]
alpha = float(sys.argv[2])

#infile = os.path.basename(tableOfCounts) # remove path
#filebase = infile.split('.')[0] # remove file extension
#outbase = '_'.join(filebase.split('_')[1:]) # remove prefix

outbase = 'youngVsOld'

valsScatter, valsHist, allInfo = readTableOfCounts(tableOfCounts)
pvals_adj = plotHist(valsHist,outbase,gaus,alpha)
sig = getSigPvals(pvals_adj,alpha,outbase,allInfo)
plotScatter(valsScatter,sig,outbase,str(alpha))
