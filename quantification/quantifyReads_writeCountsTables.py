# rmf 9.21.2018 Eileen's birthday!, last modified 7.5.2019
#!/local/cluster/bin/python
import sys, os, re
import pysam
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plot
import numpy as np

#########
# USAGE #
#########

# USAGE
if len(sys.argv) != 5 or "-h" in sys.argv or "--help" in sys.argv:
    print "USAGE: python " + os.path.basename(__file__) + " <BAM file list> <GFF file> <outbase> <window around start site (in bp)>"
    sys.exit()

###############
# SUBROUTINES #
###############

# Reads in GFF file and splits on tab to assign values to dictionary
def readGFFfile(gffFile):
    namePattern = re.compile('Name=(.*?);')
    IDpattern = re.compile('ID=(.*?);') # for dme-miR-10404
    miRs = {}
    mirList = []
    for line in open(gffFile,'r'):
        if '#' not in line: # skip headers
            seqID, source, gtfType, start, stop, score, strand, phase, attributes = line.strip().split('\t')
            if gtfType == 'miRNA':
                match1 = namePattern.search(attributes) # search attributes for miR name according to pattern, above
                name = match1.group(1) # assign miR name to a variable
                match2 = IDpattern.search(attributes) # some mirs have forms from multiple chromosomes
                ID = match2.group(1) # this ID ensures each mirID will be unique
                mirName = name + '_' + ID
                if seqID not in miRs:
                    miRs[seqID] = []
                miRs[seqID].append((int(start),int(stop),strand,mirName))
                mirList.append(mirName) # get simple list for writing all names to outfile
    return miRs, mirList

def readBamListFile(bamListFile,miRs,window):
    rawCounts = {}
    adjCounts = {}
    adjRPM = {}
    totalCounts = {}
    heterogeneity = {}
    with open(bamListFile,'r') as bamListFile:
        for line in bamListFile:
            ZT, bamFile = line.strip().split('\t')
            rawCounts[ZT],adjCounts[ZT],adjRPM[ZT],totalCounts[ZT],heterogeneity[ZT] = getBamReadCounts(bamFile,miRs,window)
            print 'BAM file for timepoint ZT', ZT, 'successfully parsed!'
    bamListFile.close()
    return rawCounts, adjCounts, adjRPM, totalCounts, heterogeneity

def getBamReadCounts(bamFile,miRs,window):
    rawCounts = {} # store raw read counts
    adjCounts = {} # store read counts adjusted using NH tags
    adjRPM = {} # store RPM calculated using NH tag adjusted counts
    heterogeneity = {}
    # Create object to read input file as a BAM file; 'rb' is read as bam                                                                                          
    bam = pysam.AlignmentFile(bamFile, 'rb')
    totalCounts = bam.mapped  # get total number of reads in the library (reads that aligned to the genome)
    print 'Total reads in BAM file successfully mapped!'
    # getting reads that map to the mirbase GFF annotation file
    for seqID in miRs:
        for miR in miRs[seqID]:
            positionList = []
            rawCount = 0.0 # start as a decimal to make sure it's float
            adjCount = 0.0
            # unpacking the tuple
            refStart,refStop,refStrand,mirID = miR
            for read in bam.fetch(seqID,refStart,refStop): # get reads matching miRbase annotation
                readStrand = getReadStrand(read)
                if readStrand == refStrand: # if the read strand matches the annotation strand
                    if refStrand == '+':
                        if refStart - window <= read.reference_start + 1 and read.reference_start +1 <= refStart + window:
                            positionList.append(read.reference_start + 1) # 5' start position
                            rawCount += 1 # count 1 for each strand-specific mapped read
                            numHits = read.get_tag("NH")
                            adjCount += 1.0/float(numHits)
                    elif refStrand == '-':
                        if refStop - window <= read.reference_end and read.reference_end <= refStop + window:
                            positionList.append(read.reference_end) # 3' end position
                            rawCount += 1
                            numHits = read.get_tag("NH")
                            adjCount += 1.0/float(numHits)
            # don't calculate if annotated miR has no reads
            if len(positionList) > 0:
                heterogeneity[mirID] = computeHeterogeneity(positionList)
            else:
                heterogeneity[mirID] = 0.0
            RPM = adjCount/(totalCounts/float(1000000))  # calculate reads per million reads in library using adjusted counts
            rawCounts[mirID] = rawCount
            adjCounts[mirID] = adjCount
            adjRPM[mirID] = RPM
    bam.close()
    return rawCounts,adjCounts,adjRPM, totalCounts, heterogeneity

# helper function to get strand of each read
def getReadStrand(read):
    if read.is_reverse:
        readStrand = '-'
    else:
        readStrand = '+'
    return readStrand

def computeHeterogeneity(positionList):
    # get most common start position
    consensusPosition = max(set(positionList), key=positionList.count)
    # count number of times that position occurs
    consensusCount = positionList.count(consensusPosition)
    # calculate fraction
    heterogeneity = 1 - (consensusCount/float(len(positionList)))
    return heterogeneity

def writeTotalCounts(readCounts,outfile):
    orderedZT = []
    orderedZT = readCounts.keys()
    orderedZT.sort(key=int)
    with open(outfile,'w') as outfile:
        for ZT in orderedZT:
            if ZT in readCounts:
                print>>outfile, ZT, '\t', readCounts[ZT]

def writeTableOfCounts(readCounts,outfile,mirList):
    orderedZT = []
    orderedZT = readCounts.keys()
    orderedZT.sort(key=int)
    header = '\t'.join(orderedZT) # gets all ZT numbers and joins with tab
    with open(outfile,'w') as outfile:
        outfile.write('#mirID\t' + header + '\n')
        for miR in mirList:
            toWrite = []
            toWrite.append(miR)
            for ZT in orderedZT:
                if ZT in readCounts:
                    if miR in readCounts[ZT]:
                        toWrite.append(str(readCounts[ZT][miR]))
                    else:
                        toWrite.append(str(0))
            outfile.write('\t'.join(toWrite) + '\n')
    outfile.close()

def writeHeterogeneityTable(heterogeneity,readCounts,outfile,mirList):
    orderedZT = []
    orderedZT = readCounts.keys()
    orderedZT.sort(key=int)
    # write header
    header = []
    for ZT in orderedZT:
        countZT = 'ZT' + ZT + '_count'
        hetZT = 'ZT' + ZT + '_heterogeneity'
        header.append(countZT)
        header.append(hetZT)
    header = '\t'.join(header) # gets all ZT numbers and joins with tab                       
    with open(outfile,'w') as outfile:
        outfile.write('#mirID\t' + header + '\n')
        for miR in mirList:
            toWrite = []
            toWrite.append(miR)
            for ZT in orderedZT:
                # write count first
                if ZT in readCounts:
                    if miR in readCounts[ZT]:
                        toWrite.append(str(readCounts[ZT][miR]))
                    else:
                        toWrite.append(str(0))
                # next write heterogeneity score
                if ZT in heterogeneity:
                    if miR in heterogeneity[ZT]:
                        toWrite.append(str(heterogeneity[ZT][miR]))
                    else:
                        toWrite.append(str(heterogeneity[ZT][miR]))
            outfile.write('\t'.join(toWrite) + '\n')
    outfile.close()

########
# MAIN #
########

# Assigns variable names to user input files
bamListFile = sys.argv[1] # tab delim file with TP and bam file
gffFile = sys.argv[2]
outbase = sys.argv[3]
window = int(sys.argv[4])

rawCountsTable = 'rawCountsTable_' + outbase + '.txt'
adjCountsTable = 'adjCountsTable_' + outbase + '.txt'
RPM_table = 'adjRPMtable_' + outbase + '.txt'
totalCountsTable = 'totalCountsTable_' + outbase + '.txt'
het_table = 'heterogeneity_' + outbase + '.txt'

# Calls subroutines, above
miRs,mirList = readGFFfile(gffFile)
print 'GFF file successfully read!'
rawCounts, adjCounts, adjRPM, totalCounts, heterogeneity = readBamListFile(bamListFile,miRs,window)
writeTotalCounts(totalCounts,totalCountsTable) # write table with total mapped counts for each ZT (not necessarily mirs)
writeTableOfCounts(rawCounts,rawCountsTable,mirList) # write raw counts table
writeTableOfCounts(adjCounts,adjCountsTable,mirList) # write adjusted counts table
writeTableOfCounts(adjRPM,RPM_table,mirList) # write table with RPM 
# write heterogeneity table with raw counts, bc these are what is used to calculate heterogeneity
writeHeterogeneityTable(heterogeneity,rawCounts,het_table,mirList)
