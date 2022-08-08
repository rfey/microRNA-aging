#!/usr/bin/python

# DH added fixed random seed 3/4/21

# PM Last updated 3/2/21

'''
Flexible, vectorized calculation of all RP Scores from tabulated time series expression data

# Current major differences from original scripts/classifyExpression_filter script:

# Advantages/what's new

-Significantly faster (even without multiprocessing)
-Can calculate RP score of any valid period (not just 24)
-Works for any (even) number of time points (may not work for odd, not tested)
-Option to compute RP scores with or without permutation tests

-Specifiable minimum median and max/minimum expression filtering criteria to 
conduct permutation test
-Specifiable option to exclude expression profiles that contain 0 from 
permutation test

# To do:

-Does not support global shuffling or z-scores for significance tests
-Some values in calculateSP/PVEC/phases are hard-coded for RP24 in context of 48 hrs in 4 hr intervals
'''

###############
### IMPORTS ###
###############

import sys, argparse
import multiprocessing as mp
from multiprocessing import Process, Manager
import numpy as np
np.random.seed(222)

###################
### SUBROUTINES ###
###################

EPS = 1e-5


def calculateSpectralParity(PSD):
    # row-wise sums
    sumEvenTerms = np.sum(PSD[:, 2::2], axis=1)+EPS
    sumOddTerms = np.sum(PSD[:, 1::2], axis=1)+EPS
    return np.log2(sumEvenTerms/sumOddTerms)

def calculatePVEC(expArray):
    # may not work if odd number of values per row
    nRows, nCols = expArray.shape
    half = nCols/2
    # 12 in 12/pi hardcoded for 24 hours 
    phases = np.angle(np.fft.fft(expArray, axis=1))*(12.0/np.pi)
    return np.std(phases[:, 2:half+1:2], axis=1)**2

def calculatePhase(expArray):
    components = np.fft.fft(expArray, axis=1)
    # index of 2 & 12/pi hard-coded for RP24
    # could not vectorize np.angle/complex functions but list comp. here is fast
    phiList = [np.angle(complex(row[2].real, row[2].imag)) for row in components]
    rawPhaseList = [-12/np.pi*phi for phi in phiList]
    phaseList = [24+rawPhase if rawPhase < 0 else rawPhase for rawPhase in rawPhaseList]
    return np.array(phaseList)
        

def computePSDVectorized(expArray):
    nRows, nCols = expArray.shape
    half = nCols/2
    DFTArray = np.abs(np.fft.fft(expArray, axis=1)/nCols)**2
    PSDArray = np.zeros(DFTArray.shape)
    # if even # of data points/columns, every DFT val from 1 to half point except half point is doubled (e.g. if 6 columns, PSD = [0, 2*DFT, 2*DFT, DFT 0, 0]).
    if nCols % 2 == 0: 
        PSDArray[:, 1:half] = 2*DFTArray[:, 1:half]
        PSDArray[:, half] = DFTArray[:, half]
    # if odd, every DFT val from 1 to half way point rounded up is doubled (e.g. if 5 columns, PSD = [0, 2*DFT, 2*DFT, 0, 0])
    else: 
        PSDArray[:, 1:half+1] = 2*DFTArray[:, 1:half+1]
    return PSDArray

def computeRPVectorized(PSDArray, k): 
    # k = frequency. If 48-hour time series data, a period of 24 = frequency of 48/24 = 2

    nRows, nCols = PSDArray.shape
    # all columns/frequencies in PSD that do not correspond to frequency k
    otherColumns = [freq for freq in range(nCols) if freq != k] 
    RPArray = np.log2((PSDArray[:, k]+EPS)/(np.sum(PSDArray[:, otherColumns], axis=1)+EPS))
    # round to 6 decimal places otherwise the output is very long
    RPArray = np.around(RPArray, 6) 
    return RPArray

def shuffleByAxis(array, axis, rep):
    '''
    extremely fast 2D array shuffling
    Reference: 'https://stackoverflow.com/questions/5040797/shuffling-numpy-array-along-a-given-axis/55317373#55317373' with modification for with rep.
    '''
    # shuffle without replacement
    if rep is False: 
        indices = np.random.rand(*array.shape).argsort(axis=axis)
    # shuffle with replacement, modification of line above
    else: 
        indices = np.random.randint(array.shape[axis], size=(array.shape))
    return np.take_along_axis(array, indices, axis=axis)


def computePvalueVectorized(expArray, RPArray, k, nTrials, useRepFlag, divideFlag=True, cut=50):    
    '''
    expArray must be 2D
    RPArray must be 1D 
    k = freq, must be an integer
    if useRepFlag is True, shuffle with replacement
    if divideFlag is False, returns array of raw counts of RP scores above RP scores in RPArray, undivided by nTrials
    ^ useful for multiprocessing where only dividing total counts at the end yields p-values
    slowest step is calculating the PSDs, but I don't think this can be easily optimized further
    '''
    # copy array
    temp = list(expArray)

    # initialize empty array of counts size of expArray
    countsArray = np.zeros(len(expArray))

    # cut the number of trials to work with at a time to save memory
    # otherwise takes a huge chunk of memory for a long time if nTrials is very large
    # cutting has very little impact on how fast the function runs
    nTrialsChunk = nTrials/cut
    if nTrials % cut == 0:
        sizes = [nTrialsChunk for x in range(cut)] 
    else:
        lastChunk = nTrials % cut # remainder
        sizes = [nTrialsChunk for x in range(cut)] 
        sizes += [lastChunk] 
        cut += 1 # add additional idx for remainder when iterating through (below)
                
    # sanity check
    if sum(sizes) == nTrials:
        for idx in range(cut):
            tempAll = np.tile(temp, (sizes[idx], 1)) # np.tiles repeats an input array
            tempAll = shuffleByAxis(tempAll, axis=1, rep=useRepFlag) # shuffle the columns of each row -> axis=1
            tempRPArray = computeRPVectorized(computePSDVectorized(tempAll), k) 
            for i in range(len(expArray)):
                # all indices that correspond to expression data at expArray[i], based on the number of trials, in tempRPArray
                allIndices = np.arange(i, len(tempRPArray), len(expArray))
                countsArray[i] += len(np.where(tempRPArray[allIndices]>RPArray[i])[0])    
        if divideFlag is True:
            pValArray = countsArray/float(nTrials)
            return pValArray
        else:
            return countsArray

def chunks(lst, n):
    # https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def computeCountsVectorizedMP(expArr, rpArr, IDs, k, nPerm, useRepFlag, outDict, divideFlag=False):
    # multiprocessing-friendly w/manager-shared outDict
    # yields counts, not p-values, since multiprocessing
    countsArray = computePvalueVectorized(expArr, rpArr, k, nPerm, useRepFlag, divideFlag=divideFlag) 
    countDict = {ID:count for ID,count in zip(IDs, countsArray)}
    outDict.update(countDict)

# copied directly from classifyExpression_filter.py
def computeQvalues(pValues):
    pValuesUnsorted = []
    for geneID in pValues:
        if pValues[geneID] != "NOTEST":
            pValuesUnsorted.append((geneID,pValues[geneID]))

    # Sort list by p value
    pValuesSorted = sorted(pValuesUnsorted,key=lambda p:p[1])

    # calculate q values from sorted list, going backwards so q values can be adjusted if needed
    qPrev = 1.0
    qValues = {}
    for i in reversed(range(len(pValuesSorted))):
        geneID,pValue = pValuesSorted[i]
        r=float(i+1.0)
        qTemp=pValue*len(pValuesSorted)/r
        qValue=min(qTemp,qPrev)
        qPrev=qValue
        qValues[geneID] = qValue

    for geneID in pValues:
        if pValues[geneID] == "NOTEST":
            qValues[geneID] = "NOTEST"

    return qValues

def calculateMaxMinRatio(expList):
    EPS = 1e-5
    maxx = max(expList)+EPS
    minn = min(expList)+EPS
    mmRatio = maxx/float(minn)
    return mmRatio

def isNonZero(expList):
    status = True
    if 0 in expList:
        status = False
    return status
     
def passesFilter(dataList, medT, maxMinT, nonZeroFlag):
    EPS = 1e-5
    filters = {'med': medT, 'mm': maxMinT, 'nz': nonZeroFlag}

    filterStatus = {}
    for criteria, val in filters.items():
        if val is not None and val is not False:
            filterStatus[criteria] = False

    # check max/min
    if filters['mm'] is not None:
        mmRatio = calculateMaxMinRatio(dataList)
        if mmRatio >= filters['mm']:
            filterStatus['mm'] = True
    # check median
    if filters['med']is not None:
        med = np.median(dataList)
        if med >= filters['med']:
            filterStatus['med'] = True
    # check nonzero condition
    if filters['nz'] is True:
        filterStatus['nz'] = isNonZero(dataList)

    filterStatuses = filterStatus.values()

    # checks if every entry in filterStatuses is True. 
    # If even one is False, returns False. If all are True, returns True
    return all(filterStatuses)

def calculateParams(expList):
    # Max/Min Ratio, Median, Nonzero
    mmRatio = calculateMaxMinRatio(expList)
    median = np.median(expList)
    nonZero = isNonZero(expList)
    if nonZero is True:
        nonZero = 1
    else:
        nonZero = 0
    return mmRatio, median, nonZero

def fillDataDicts(expArray, psdArray, IDs, freqs, periods, minRP, filters, computeSigFlag, nProcessors):
    overallRPDict = {}
    overall_pValDict = {}
    overall_qValDict = {}

    ## fill RP, pVal, and qVal dictionaries
    for freq, period in zip(freqs, periods):
        rpArray = computeRPVectorized(psdArray, freq)
        rpDict = {ID:rpScore for ID,rpScore in zip(IDs, rpArray)}
        if computeSigFlag is True:
            pValDict = {}
            pValDictToTest = {}
            pValDictToNotTest = {}

            # indices where RP scores are above min threshold for computing significance
            indicesAbove = np.where(rpArray>=minRP)[0]
            # indices that pass filters AND are above RP thresh
            indicesToTest = [idx for idx in indicesAbove if passesFilter(expArray[idx],*filters)]
            
            # use filtered indices to filter RP array, Exp array, and IDs
            # could be empty; none could be above thresh
            if len(indicesToTest) > 0: 
                rpArrayToTest = rpArray[indicesToTest]
                expArrayToTest = expArray[indicesToTest]
                IDsToTest = IDs[indicesToTest]

                if nProcessors > 1:
                    size = int(len(rpArrayToTest)/float(nProcessors))+1
                    # divide into chunks approx. equally sized for each processor
                    rpChunks, expChunks, IDsChunks = [list(chunks(arr, size)) for arr in (rpArrayToTest, expArrayToTest, IDsToTest)]

                    manager = Manager()
                    # proxy dict object so that each processor can work with same dict. Does not behave identically like Py dict but close
                    countDictToTest = manager.dict()
                    processList = []

                    for rpChunk, expChunk, IDsChunk in zip(rpChunks, expChunks, IDsChunks):
                        process = Process(target=computeCountsVectorizedMP, args=(expChunk, rpChunk, IDsChunk, freq, nPerm, useRepFlag, countDictToTest, False))
                        process.start()
                        processList.append(process)

                    for process in processList:
                        process.join()

                    # processes have finished by this point
                    # proxy dict object cannot be iterated through, so convert to Python dict
                    countDictToTest = dict(countDictToTest) 

                    # divide counts by # permutations to yield p-values
                    for ID in countDictToTest: 
                        pValDictToTest[ID] = countDictToTest[ID]/float(nPerm) 

                else: # if not multiprocessing
                    pValArrayToTest = computePvalueVectorized(expArrayToTest, rpArrayToTest, freq, nPerm, useRepFlag, divideFlag=True)
                    pValDictToTest = {IDA:pValA for IDA,pValA in zip(IDsToTest, pValArrayToTest)}

            # indices below thresholds that are not tested for significance
            # indices where RP scores are below min threshold
            indicesBelow = list(np.where(rpArray<minRP)[0])

            # indices that fail filter
            indicesFailFilter = [idx for idx,row in enumerate(expArray) if not passesFilter(row, *filters)]

            # indices that are below min RP thresh OR fail filter
            indicesToNotTest = np.array(tuple(set(indicesBelow+indicesFailFilter)))

            # could be empty; all could pass criteria
            if len(indicesToNotTest) > 0:
                IDsToNotTest = IDs[indicesToNotTest]
                pValDictToNotTest = {IDB: 'NOTEST' for IDB in IDsToNotTest}

            pValDict.update(pValDictToTest)
            pValDict.update(pValDictToNotTest)

            overall_pValDict[period] = pValDict

            qValDict = computeQvalues(pValDict)        
            overall_qValDict[period] = qValDict

        overallRPDict[period] = rpDict

    return overallRPDict, overall_pValDict, overall_qValDict



###############################
### USAGE/INPUT SUBROUTINES ###
###############################

def strToBool(s): # returns booleans for strings beginning with corresponding substrings, for use with commandline arguments
    status = None
    s = s.lower()

    trueStrings = ["t", "true", "1"]
    falseStrings = ["f", "false", "0"]
    for string in trueStrings:
        if s.startswith(string):
            status = True
    for string in falseStrings:
        if s.startswith(string):
            status = False
    if status is not None:
        return status
    else: # exit condition
        print "Error - Invalid boolean string. Expected t/true/1 for True or f/false/0 for False (case insensitive)"
        sys.exit()

def isFloatable(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def readExpressionFile(inFile, symbolsFlag):
    expDict = {}
    with open(inFile) as F:
        for i, line in enumerate(F):
            if i > 1: # skip header
                columns = line.strip().split('\t')
                if symbolsFlag is True:
                    ID, symbol = columns[:2]
                    expVals = map(float, filter(isFloatable, columns[2:]))
                else:
                    ID = columns[0]
                    symbol = ID
                    expVals = map(float, filter(isFloatable, columns[1:]))
                
                key  = '_'.join((ID,symbol))
                expDict[key] = expVals
                      
    nCols = len(expVals)
    return expDict, nCols

def readSymbolFile(inFile):
    idToSymbolDict = {}
    with open(inFile) as F:
        for line in F:
            ID, symbol = line.strip().split('\t')
            idToSymbolDict[ID] = symbol
    return idToSymbolDict

def configureIDsSymbols(idsAndSymbols, symbolsFlag, convertIDsFlag):
    IDs, symbols = zip(*[idS.split('_') for idS in idsAndSymbols])
    if symbolsFlag is False and convertIDsFlag is True:
        # IDs in symbolDict
        inSym = [IDin for IDin in IDs if IDin in symbolDict]
        # IDs not in symbolDict
        outSym = [IDout for IDout in IDs if IDout not in symbolDict]

        symbols = [symbolDict[ID] for ID in inSym]
        # check
        if len(IDs) != len(symbols):
            print "Error - not all IDs could be converted to symbols"
            print "List of IDs that could not be converted:"
            print outSym
            sys.exit()

    if symbols == IDs:
        symbols = None

    IDs = np.array(IDs)
    if symbols is not None:
        symbols = np.array(symbols)    

    return IDs, symbols

def configureFreqsPeriods(freqs, nCols, maxTime):
    # get frequencies and periods based on number of columns in exp profiles
    if freqs != 'all':
        # split by comma, convert to integers
        if type(freqs) == str: 
            freqs = map(int, freqs.split(','))
        else:
            freqs = [freqs]
            # filter out freqs equal to or below 0 and above n/2 (Nyquist freq)
            # not sure if this works properly for odd number of columns
        freqs = filter(lambda n:0 < n <=nCols/2, freqs)

    else: # all valid
        freqs = range(1,nCols/2+1)

    if freqs:
        freqs.sort() # sort frequencies from lowest to highest 
    else:
        print "Error - invalid frequency specified. Must be a positive integer less than or equal to half the number of time points"
        sys.exit()

    # since freqs are sorted from low to high, periods are sorted from high to low 
    periods = [maxTime/float(freq) for freq in freqs]
    # convert x.0 floats to ints e.g. 24.0 to 24 for cleaner labels
    periods = map(cleanFloat, periods) 
    return freqs, periods

def cleanFloat(val): # ex: if float value is 5.0, convert to int 5
    if float(val) - int(val) == 0:
        return int(val)
    else:
        return float(val)


def getArgs():
    parser = argparse.ArgumentParser(description="Compute RP scores for specified periods of time series expression data. Time interval between exp. points assumed to be equal.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # required
    parser.add_argument("--exp", "-e", type=str, required=True, help="Expression file",default=False)
    parser.add_argument("--outBase", "-oB", type=str, required=True, help="outBase",default=False)
    
    # optional; multiprocessing
    parser.add_argument("--nProcessors", "-mp", '-nPro', type=int, required=False, help="Number of processors if more than 1. Solely used for permutation tests.",default=1)
    
    # optional; specify frequences
    parser.add_argument('--freqs', '-f', '-k', type=str, required=False, help="Frequencies (k) to calculate RP score of. Must be below Nyquist Freq (n/2). Comma separated if more than one (e.g. 1,2,3,4) Can pass 'all' to do all up to k=n/2.", default=2)
    parser.add_argument('--maxTime', '-mt', type=float, required=False, help="Maximum time; used to calculate periods (maxTime/freq = period)", default=48)

    # optional; related to q-value computation
    parser.add_argument("--computeSig", "-sig", type=strToBool, required=False, help="T/True/1 to compute significance of RP score(s) with permutation tests.", default=False)
    parser.add_argument("--minRP", "-mRP", type=float, required=False, help="Minimum RP score to compute permutation test.", default=0)
    parser.add_argument("--nPerm", "-n", type=str, required=False, help="Number of permutations to conduct.", default=10000)
    parser.add_argument("--useRep", "-r", type=strToBool, required=False, help="Shuffle expression with replacement. T/True/1 to shuffle with replacement", default=False)

    # filtering exp before q-value computation
    parser.add_argument("--nonZero", "-nz", type=strToBool, required=False, help="T/True/1 to exclude expression profiles that contain 0 from significance tests.", default=False)
    parser.add_argument("--med", "-med", type=float, required=False, help="Minimum median expression to compute permutation test for significance. Rows above threshold excluded.", default=None)
    parser.add_argument("--minMax", "-mm", type=float, required=False, help="Minimum max/min ratio of expression to compute permutation test for significance. Rows above threshold.", default=None)

    # optional; related to conversion of ID to symbol
    parser.add_argument("--symbol", "-S", type=strToBool, required=False, help='F/False/0 if symbols NOT present in second column', default=True)    
    parser.add_argument("--symbolFile", "-sf", type=str, required=False, help='ID to Symbol map file', default='/nfs0/BB/Hendrix_Lab/Drosophila/Annotations/symbols_r6.21.txt')    
    parser.add_argument("--convertIDs", "-cIDs", type=strToBool, required=False, help="Convert IDs in first column to symbols and add second column of symbols to outfile. F/False/0 to not do so", default=True)
    args = parser.parse_args()
    return args

##############
### INPUTS ###
##############

args = getArgs()

# required, no defaults
expFile = args.exp
outBase = args.outBase
outFile = outBase+'.txt'

# time, freqs
maxTime = args.maxTime # float
freqs = args.freqs

# multiprocessing
nProcessors = args.nProcessors
minPro = 1
maxPro = mp.cpu_count()

if nProcessors > maxPro:
    print 'Error - # of processors specified:', nProcessors, 'exceeds available:', maxPro
    sys.exit()
if nProcessors < minPro:
    print 'Error - # of processors specified:', nProcessors, 'must be at least 1'
    sys.exit()

# permutation tests
computeSigFlag = args.computeSig
minRP = args.minRP

nPerm = int(float(args.nPerm)) # supports scientific notation e.g. 1e4 for 10,000
useRepFlag = args.useRep

# pre-test expression filtering
medThresh = args.med
minMaxThresh = args.minMax
nonZeroFlag = args.nonZero

filters = (medThresh, minMaxThresh, nonZeroFlag)

# IDs, symbols
symbolsFlag = args.symbol
symbolFile = args.symbolFile
convertIDsFlag = args.convertIDs

############
### MAIN ###
############

symbolDict = readSymbolFile(symbolFile)
expDict, nCols = readExpressionFile(expFile, symbolsFlag)
freqs, periods = configureFreqsPeriods(freqs, nCols, maxTime)
  
# unpack expDict 
idsAndSymbols, expArray = zip(*expDict.items())
expArray = np.array(expArray)

IDs, symbols = configureIDsSymbols(idsAndSymbols, symbolsFlag, convertIDsFlag)

psdArray = computePSDVectorized(expArray)

PVECArray = calculatePVEC(expArray)
spectralParityArray = calculateSpectralParity(psdArray)
phaseArray = calculatePhase(expArray)


PVECDict = {IDs[i]:PVECArray[i] for i in range(len(IDs))}
spectralParityDict = {IDs[i]:spectralParityArray[i] for i in range(len(IDs))}
phaseDict = {IDs[i]:phaseArray[i] for i in range(len(IDs))}
# Max/Min Ratio, Median, NonZero status of expression
expValsParams = {IDs[i]:calculateParams(expArray[i]) for i in range(len(IDs))}

# calculate RP Scores, and p/q values with permutations (if applicable)
overallRPDict, overall_pValDict, overall_qValDict = fillDataDicts(expArray, psdArray, IDs, freqs, periods, minRP, filters, computeSigFlag, nProcessors)


## write to outfile    
with open(outFile, 'w') as O:    
    # prepare & write header first
    if symbols is None:
        headerList = ['ID']
    else:
        headerList = ['ID', 'Symbol']
    rpHeaderList = ['RP' + str(period) for period in periods]

    if computeSigFlag is True:
        qHeaderList = ['qVal' for period in periods]
        pHeaderList = ['pVal' for period in periods]
        for rpH, qH, pH in zip(rpHeaderList, qHeaderList, pHeaderList):
            headerList += [rpH, qH, pH]
    else:
        for rpH in rpHeaderList:
            headerList.append(rpH)
        
    headerList.append('SP')
    headerList.append('PVEC')
    headerList.append('Phase')
    headerList.append("Max/Min")
    headerList.append("Median")
    headerList.append("Non-Zero")

    O.write('\t'.join(headerList)+'\n')
    
    # write data    
    for i, ID in enumerate(IDs):
        # non-RP expression params 

        # number of digits to round to
        n = 6
        SP = str(round(spectralParityDict[ID],n))
        PVEC = str(round(PVECDict[ID],n))
        Phase = str(round(phaseDict[ID],n))
        mmRatio, median, nonZero = expValsParams[ID]
        mmRatio = round(mmRatio,n)
        median = round(median,n)
        mmRatio, median, nonZero = map(str, (mmRatio,median,nonZero))
        
        data = []
        for period in periods:
            rpScore = overallRPDict[period][ID]
            if computeSigFlag is True:
                pVal = overall_pValDict[period][ID]
                qVal = overall_qValDict[period][ID]
                if pVal != 'NOTEST':
                    pVal = round(pVal, 6)
                if qVal != 'NOTEST':
                    qVal = round(qVal, 6)

                data.append('\t'.join(map(str, (rpScore, qVal, pVal))))
            else:

                data.append(str(rpScore))

        data.append(SP)
        data.append(PVEC)
        data.append(Phase)
        data.append(mmRatio)
        data.append(median)
        data.append(nonZero)

        if symbols is None:
            O.write('\t'.join((ID, '\t'.join(data))))            
        else: 
            symbol = symbols[i]
            O.write('\t'.join((ID, symbol, '\t'.join(data))))

        if i < len(IDs):
            O.write('\n')

                            

            


