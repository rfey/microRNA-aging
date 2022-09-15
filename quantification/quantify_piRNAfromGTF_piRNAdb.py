# rmf 9.21.2018 Eileen's birthday!, last modified 11.19.2021
#!/local/cluster/bin/python
import sys, re
import pysam
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plot
import numpy as np

#########
# USAGE #
#########

# USAGE
if len(sys.argv) != 4 or "-h" in sys.argv or "--help" in sys.argv:
    print "USAGE: python " + sys.argv[0] + " <BAM file list> <GFF file> <outbase>"
    exit()

###############
# SUBROUTINES #
###############

# Reads in GFF file and splits on tab to assign values to dictionary
def readGFFfile(gffFile):
    piRNAs = {}
    piRNA_clusters = []
    for line in open(gffFile,'r'):
        if 'Cluster_Acess' not in line: # skip headers
            clusterID, chrom, start, stop, strand, buildCode = line.strip().split('\t')
            if chrom not in piRNAs:
                piRNAs[chrom] = []
            piRNAs[chrom].append((int(start),int(stop),strand,clusterID))
            piRNA_clusters.append(clusterID) # get simple list for writing all names to outfile
    return piRNAs, piRNA_clusters

def readBamListFile(bamListFile, piRNAs):
    rawCounts = {}
    adjCounts = {}
    adjRPM = {}
    totalCounts = {}
    senseLengths, antisenseLengths = {},{}
    allLengths = []
    with open(bamListFile,'r') as bamListFile:
        for line in bamListFile:
            ZT, bamFile = line.strip().split('\t')
            print ZT, bamFile
            rawCounts[ZT],adjCounts[ZT],adjRPM[ZT],totalCounts[ZT],senseLengths[ZT],antisenseLengths[ZT],allLengths = getBamReadCounts(bamFile,piRNAs,allLengths)
            print 'BAM file for timepoint ZT', ZT, 'successfully parsed!'
    bamListFile.close()
    return rawCounts, adjCounts, adjRPM, totalCounts, senseLengths, antisenseLengths, allLengths

def getBamReadCounts(bamFile, piRNAs, allLengths):
    rawCounts = {} # store raw read counts
    adjCounts = {} # store read counts adjusted using NH tags
    adjRPM = {} # store RPM calculated using NH tag adjusted counts
    totalSense, totalAntisense = 0.0, 0.0
    senseLengths, antisenseLengths = {},{}
    # Create object to read input file as a BAM file; 'rb' is read as bam
    bam = pysam.AlignmentFile(bamFile, 'rb')
    totalCounts = bam.mapped  # get total number of reads in the library (reads that aligned to the genome)
    print 'Total reads in BAM file successfully mapped!'
    # getting reads that map to the piRNA GFF annotation file
    for chrom in piRNAs:
        for cluster in piRNAs[chrom]:
            rawSense, rawAntisense = 0.0, 0.0  # start as a decimal to make sure it's float
            adjSense, adjAntisense = 0.0, 0.0
            # unpacking the tuple
            refStart,refStop,refStrand,clusterID = cluster
            for read in bam.fetch(chrom,refStart,refStop): # get reads matching GFF annotation
                length = read.query_length
                if length not in allLengths:
                    allLengths.append(length)
                readStrand = getReadStrand(read)
                if readStrand == refStrand: # if the read strand matches the annotation strand
                    rawSense += 1 # count 1 for each strand-specific mapped read
                    numHits = read.get_tag("NH")
                    adjSense += 1.0/float(numHits)
                    totalSense += 1
                    if length not in senseLengths:
                        senseLengths[length] = 0
                    senseLengths[length] += 1
                else:
                    rawAntisense += 1
                    numHits = read.get_tag("NH")
                    adjAntisense += 1.0/float(numHits)
                    totalAntisense += 1
                    if length not in antisenseLengths:
                        antisenseLengths[length] = 0
                    antisenseLengths[length] += 1
            senseRPM = adjSense/(totalCounts/float(1000000))  # calculate reads per million reads in library using adjusted counts
            antisenseRPM = adjAntisense/(totalCounts/float(1000000))
            rawCounts[clusterID] = (rawSense, rawAntisense)
            adjCounts[clusterID] = (adjSense, adjAntisense)
            adjRPM[clusterID] = (senseRPM, antisenseRPM)
    bam.close()
    return rawCounts,adjCounts,adjRPM, (totalCounts, totalSense, totalAntisense), senseLengths, antisenseLengths, allLengths

# helper function to get strand of each read
def getReadStrand(read):
    if read.is_reverse:
        readStrand = '-'
    else:
        readStrand = '+'
    return readStrand

def writeLengthTable(lengths, outfile, allLengths):
    orderedZT = []
    orderedZT = lengths.keys()
    orderedZT.sort(key=int)
    allLengths.sort(key=int)
    with open(outfile, 'w') as outfile:
        outfile.write('#length\t' + '\t'.join(orderedZT) + '\n')
        for length in allLengths:
            toWrite = [str(length)]
            for ZT in orderedZT:
                if length in lengths[ZT]:
                    n = lengths[ZT][length]
                else:
                    n = 0
                toWrite.append(str(n))
            outfile.write('\t'.join(toWrite) + '\n')
    outfile.close()

def writeTotalCounts(readCounts,outfile):
    orderedZT = []
    orderedZT = readCounts.keys()
    orderedZT.sort(key=int)
    with open(outfile,'w') as outfile:
        outfile.write('ZT\ttotalCounts\ttotalSenseCounts\ttotalAntisenseCounts\n')  # header
        for ZT in orderedZT:
            if ZT in readCounts:
                total, sense, antisense = readCounts[ZT]
                outfile.write(ZT + '\t' + str(total) + '\t' + str(sense) + '\t' + str(antisense) + '\n')
    outfile.close()

def writeTableOfCounts(readCounts,outfile,piRNA_clusters):
    orderedZT = []
    orderedZT = readCounts.keys()
    orderedZT.sort(key=int)
    h = orderedZT * 2
    header = '\t'.join(h) # gets all ZT numbers and joins with tab
    with open(outfile,'w') as outfile:
        outfile.write('#clusterID\t' + header + '\n')
        for clusterID in piRNA_clusters:
            toWrite = []
            toWrite.append(clusterID)
            # sense counts
            for ZT in orderedZT:
                if ZT in readCounts:
                    if clusterID in readCounts[ZT]:
                        toWrite.append(str(readCounts[ZT][clusterID][0]))
                    else:
                        toWrite.append(str(0))
            # antisense counts
            for ZT in orderedZT:
                if ZT in readCounts:
                    if clusterID in readCounts[ZT]:
                        toWrite.append(str(readCounts[ZT][clusterID][1]))
                    else:
                        toWrite.append(str(0))
            outfile.write('\t'.join(toWrite) + '\n')
    outfile.close()

########
# MAIN #
########

# Assigns variable names to user input files
bamListFile = sys.argv[1] # tab delim file with TP and bam file
gffFile = sys.argv[2]
outbase = sys.argv[3]

rawCountsTable = 'rawCountsTable_' + outbase + '.txt'
adjCountsTable = 'adjCountsTable_' + outbase + '.txt'
RPM_table = 'adjRPMtable_' + outbase + '.txt'
totalCountsTable = 'totalCountsTable_' + outbase + '.txt'
senseLengthsTable = 'senseLengthsTable_' + outbase + '.txt'
antisenseLengthsTable = 'antisenseLengthsTable_' + outbase + '.txt'

# Calls subroutines, above
piRNAs, piRNA_clusters = readGFFfile(gffFile)
print 'Read in ' + str(len(piRNA_clusters)) + ' piRNA clusters.'

rawCounts, adjCounts, adjRPM, totalCounts, senseLengths, antisenseLengths, allLengths = readBamListFile(bamListFile, piRNAs)

writeLengthTable(senseLengths, senseLengthsTable, allLengths)
writeLengthTable(antisenseLengths, antisenseLengthsTable, allLengths)
writeTotalCounts(totalCounts,totalCountsTable) # write table with total mapped counts for each ZT (not necessarily mirs)
writeTableOfCounts(rawCounts,rawCountsTable,piRNA_clusters) # write raw counts table
writeTableOfCounts(adjCounts,adjCountsTable,piRNA_clusters) # write adjusted counts table
writeTableOfCounts(adjRPM,RPM_table,piRNA_clusters) # write table with RPM 
