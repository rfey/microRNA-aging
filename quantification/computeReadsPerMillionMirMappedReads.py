# rmf 6.20.2019, last modified 1.10.2022

import sys, os

usage = 'python ' +sys.argv[0] + ' <adjCountsTable>'
if len(sys.argv) != 2 or '-h' in sys.argv or '--help' in sys.argv:
    print usage
    sys.exit()
if 'combined' not in sys.argv[1]:
    print 'Must supply adjCountsTable with both known and novel mirs.'
    exit()

# SUBROUTINES
def readInfile(table):
    mirs, sums = {},{}
    vals0,vals4,vals8,vals12,vals16,vals20 = [],[],[],[],[],[]
    with open(table,'r') as table:
        next(table)  # skip header
        for line in table:
            mir, ZT0, ZT4, ZT8, ZT12, ZT16, ZT20 = line.strip().split('\t')
            vals0.append(float(ZT0))
            vals4.append(float(ZT4))
            vals8.append(float(ZT8))
            vals12.append(float(ZT12))
            vals16.append(float(ZT16))
            vals20.append(float(ZT20))
            mirs[mir] = (float(ZT0),float(ZT4),float(ZT8),float(ZT12),float(ZT16),float(ZT20))
    sums['ZT0'] = sum(vals0)
    sums['ZT4'] = sum(vals4)
    sums['ZT8'] = sum(vals8)
    sums['ZT12'] = sum(vals12)
    sums['ZT16'] = sum(vals16)
    sums['ZT20'] = sum(vals20)
    table.close()
    return mirs, sums

def writeOutfile(mirs, sums, outfile):
    with open(outfile,'w') as outfile:
        outfile.write('#mirID\t0\t4\t8\t12\t16\t20\n') # header
        for mir in mirs:
            ZT0, ZT4, ZT8, ZT12, ZT16, ZT20 = mirs[mir]
            val0 = ZT0/(sums['ZT0']/float(1000000))
            val4 = ZT4/(sums['ZT4']/float(1000000))
            val8 = ZT8/(sums['ZT8']/float(1000000))
            val12 = ZT12/(sums['ZT12']/float(1000000))
            val16 = ZT16/(sums['ZT16']/float(1000000))
            val20 = ZT20/(sums['ZT20']/float(1000000))
            outfile.write('%s\t%f\t%f\t%f\t%f\t%f\t%f\n' % (mir,val0,val4,val8,val12,val16,val20))
    outfile.close()

# ARGUMENTS and MAIN
infile = sys.argv[1]  # ex: adjCountsTable_oldRep2_combined.txt

outbase = '_'.join(infile.split('_')[1:])
outfile = 'adjRPMMirReads_' + outbase

mirs, sums = readInfile(infile)
writeOutfile(mirs, sums, outfile)
