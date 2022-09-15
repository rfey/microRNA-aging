# rmf 7.31.2022

# quantification scripts:
quantifyReads_writeCountsTables.py
quantify_piRNAfromGTF_piRNAdb.py
usage: python quantifyReads_writeCountsTables.py <BAM file list> <GFF file> <outbase> <window around start site (in bp)>
<BAM file list> example: bamListFile_d5youngRep1
<GFF file>: mirbase GFF with chromosome names changed to the GenBank format
<window>: allows for a flexible start site for mapping reads to coordinates

4 table types produced by this script:
raw counts table: counts not adjusted for number of hits to the genome
adj counts table: counts adjusted for number of hits to the genome using NH tags
adj RPM table: Reads Per Million calculated using adj counts
total counts table: total number of reads that map to genome for each time point


# script to compute RPMMM (reads per million mapped mirs)
computeReadsPerMillionMirMappedReads.py
usage: python computeReadsPerMillionMirMappedReads.py <adjCountsTable>
<adjCountsTable>: should include both known and novel mirs for accurate normalization


# plot scripts
# used to produce box plots in Figure 3A and 3B
boxPlot_globalExpr.py
usage: python boxPlot_globalExpr.py <fileList adjCountsTables> <fileList totalCountsTables>

boxPlot_globalExpr.py
usage: python boxPlot_globalExpr.py <fileList adjCountsTables> <fileList totalCountsTables>
