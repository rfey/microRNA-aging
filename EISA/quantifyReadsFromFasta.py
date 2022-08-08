# rmf 7.17.2018, last modified 5.1.2019

import os, sys, re
import pysam
from Bio import SeqIO
from Bio.Seq import Seq

# USAGE: takes two BAM files and quantifies intronic and exonic reads
usage = 'python ' + os.path.basename(__file__) + ' <transcript (exon) BAM file> <genomic (intron) BAM file> <flybase transcript fasta file> <flybase intron fasta file> <transcriptIDtoSymbol.txt> <out prefix>\nUSAGE: takes two BAM files and quantifies exonic and intronic reads.'
if len(sys.argv) != 7 or '-h' in sys.argv or '--help' in sys.argv:
    print usage
    sys.exit()

# SUBROUTINES
def readMapFile(mapFile):
    maps = {}
    with open(mapFile,'r') as f:
        for line in f:
            trID, symbol = line.strip().split('\t')
            maps[trID] = symbol
    f.close()
    return maps

# get list of transcript IDs from transcript fasta
def readExonFasta(exonFasta):
    exonLengths = {}
    transcriptIDs = []
    lengthPattern = re.compile('length=(\d*?);')
    for record in SeqIO.parse(exonFasta,'fasta'):
        match = lengthPattern.search(record.description)
        length = match.group(1)
        transcriptIDs.append(record.id)
        exonLengths[record.id] = int(length)
    return transcriptIDs,exonLengths

def readIntronFasta(intronFasta,intronBam):
    introns = {}
    intronLengths = {}
    for record in SeqIO.parse(intronFasta,'fasta'):
        chrom,start,stop,strand,transcripts = parseRecord(record) # transcripts is a list
        length = abs(int(start)-int(stop)) # get length of intron
        readCount = getIntronCounts(intronBam,chrom,int(start),int(stop),strand) # feed in one intron at a time, get all reads mapping to that intron
        for transcript in transcripts:
            if transcript not in introns:
                introns[transcript] = 0  # initialize with zero counts 
            if transcript not in intronLengths:
                intronLengths[transcript] = 0
            introns[transcript] += readCount # key by transcript ID; each ID gets total reads mapped to this intron
            intronLengths[transcript] += int(length)
    return introns,intronLengths

def parseRecord(record):
    transcripts = []
    coordsPattern = re.compile('loc=(.*):.*?(\d+)..(\d+)') # chrom, start, stop
    chrom, loc, parent, MD5, release, species, length = record.description.strip().split('; ') # split on semicolon and space
    if 'complement' not in loc:  # get strand info
        strand = '+'
    else:
        strand = '-'
    match = coordsPattern.search(loc)
    chrom = match.group(1)
    start = match.group(2)
    stop = match.group(3)
    transcripts = parent.strip().split(',')[1:] # split on comma, return all but first element (first element is parent gene ID)
    return chrom, start, stop, strand, transcripts

# exons have been aligned to a transcript fasta file
def getExonCounts(exonSam,transcriptIDs):
    exons = {}
    samFile = pysam.AlignmentFile(exonSam, 'rb')
    # get exonic reads
    for transcriptID in transcriptIDs: # loop through transcriptIDs and pass in as chrom
        count = 0
        # the chrom is FBID in the bam file
        # transcript bam file also have NH tags added
        for read in samFile.fetch(transcriptID):  # give just FBID, which is chrom in bamfile (gets all reads aligning to that transcript)
            # check strand-- require is_reverse == false (all transcripts in alignment fasta were forward)
            if not read.is_reverse:
                count += 1
            exons[transcriptID] = count # store FBID and number of reads mapping to it
    samFile.close()
    return exons

# introns have been aligned to the genome
def getIntronCounts(intronBam,chrom,start,stop,strand):
    bamFile = pysam.AlignmentFile(intronBam, 'rb')
    # get intronic reads
    readCount = 0
    for read in bamFile.fetch(chrom,start,stop):
        readStrand = getReadStrand(read)
        # check strand and make sure it's the same as the annotation; don't count the read otherwise
        if readStrand == strand:
            readCount += 1
    bamFile.close()
    return readCount

def getReadStrand(read):
    if read.is_reverse:
        readStrand = '-'
    else:
        readStrand = '+'
    return readStrand

def writeCountsTable(exons,exonLengths,introns,intronLengths,transcriptIDs,outfile,maps):
    with open(outfile,'w') as outfile:
        outfile.write("%s\t%s\t%s\t%s\t%s\n" % ("Flybase ID","Exon Count","Exon Length","Intron Count","Intron Length"))
        for ID in transcriptIDs:
            if ID in maps:
                symbol = maps[ID]
            else:
                symbol = 'symbol_not_found'
            exonCount = 0
            intronCount = 0
            if ID in exons:
                exonCount = exons[ID]
            else: 
                exonCount = 0
            if ID in introns:
                intronCount = introns[ID]
            else:
                intronCount = 0
            if ID in exonLengths:
                exonLength = exonLengths[ID]
            else:
                exonLength = 0
            if ID in intronLengths:
                intronLength = intronLengths[ID]
            else:
                intronLength = 0
            outfile.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (ID,symbol,exonCount,exonLength,intronCount,intronLength))
    outfile.close()


# ARGUMENTS and MAIN
exonBam = sys.argv[1]
intronBam = sys.argv[2]
exonFasta = sys.argv[3]
intronFasta = sys.argv[4]
mapFile = sys.argv[5]
outPrefix = sys.argv[6]

outfile = 'tableOfCounts_' + outPrefix + '.txt'

maps = readMapFile(mapFile)
transcriptIDs,exonLengths = readExonFasta(exonFasta) # get exon coordinates from flybase exon fasta file
introns,intronLengths = readIntronFasta(intronFasta,intronBam) # get intron read counts
exons = getExonCounts(exonBam,transcriptIDs)
writeCountsTable(exons,exonLengths,introns,intronLengths,transcriptIDs,outfile,maps)
