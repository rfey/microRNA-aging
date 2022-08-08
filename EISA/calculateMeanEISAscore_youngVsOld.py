# rmf 2.28.2019, last modified 5.20.2022

import sys
import numpy as np

# USAGE
usage = 'python ' + sys.argv[0] + ' file lists: <tables of mean counts young> <tables of mean counts old>'
if len(sys.argv) != 3 or '-h' in sys.argv or '--help' in sys.argv:
    print usage
    exit()

# SUBROUTINES
def readFileList(f):
    counts = {}
    with open(f, 'r') as f:
        for line in f:
            table = line.strip()
            counts = readTable(table, counts)
    f.close()
    return counts

def readTable(table,counts):
    print 'Reading', table, '...'
    with open(table,'r') as table:
        next(table) # skip header row
        for line in table:
            FBID, symbol, exon, exonLen, intron, intronLen = line.strip().split('\t')
            key = (FBID, symbol)
            if key not in counts:
                counts[key] = []
            counts[key].append((float(exon), float(intron)))
    return counts

def sumCounts(counts):
    sums = {}
    for key in counts:
        exons, introns = zip(*counts[key])
        exon_sum = sum(exons)
        intron_sum = sum(introns)
        sums[key] = (exon_sum, intron_sum)
    return sums

def calculateEISAscore(ySums, oSums):
    all_counts = {}
    keys = list(set(ySums.keys() + oSums.keys()))
    for key in keys:
        # add pseudocount so score can be computed for every transcript
        yExon = 1 + float(ySums[key][0])
        oExon = 1 + float(oSums[key][0])
        yIntron = 1 + float(ySums[key][1])
        oIntron = 1 + float(oSums[key][1])

        # calculate log fold changes
        exonFC = np.log2(oExon/yExon)
        intronFC = np.log2(oIntron/yIntron)
        EISA = exonFC - intronFC

        all_counts[key] = (yExon,oExon,yIntron,oIntron,exonFC,intronFC,EISA)

    return all_counts

def writeOutfile(outfile,counts):
    with open(outfile,'w') as outfile:
        outfile.write('#FBID\tsymbol\tyoungExon\toldExon\texonFC(o/y)\tyoungIntron\toldIntron\tintronFC(o/y)\tEISA_score\n')
        for key in counts:
            FBID, symbol = key
            to_write = [str(x) for x in counts[key]]  # convert everything to string for writing
            yExon,oExon,yIntron,oIntron,exonFC,intronFC,EISA = to_write
            outfile.write(FBID+'\t'+symbol+'\t'+yExon+'\t'+oExon+'\t'+exonFC+'\t'+yIntron+'\t'+oIntron+'\t'+intronFC+'\t'+EISA+'\n')
    outfile.close()

# ARGUMENTS and MAIN
yFileList = sys.argv[1]
oFileList = sys.argv[2]

outfile = 'youngVsOldEISAtable.txt'

yCounts = readFileList(yFileList)
oCounts = readFileList(oFileList)
print 'Counts saved for', len(yCounts), 'transcripts in young.'
print 'Counts saved for', len(yCounts), 'transcripts in old.'

print 'Summing counts...'
ySums = sumCounts(yCounts)
oSums = sumCounts(oCounts)

print 'Calculating EISA scores...'
EISAcounts = calculateEISAscore(ySums, oSums)

print 'Writing output file...'
writeOutfile(outfile,EISAcounts)
