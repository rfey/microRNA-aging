# rmf 7.31.2022

# input file example:
inputFileExample_young.txt
 
header must include '#geneID' as shown
numbers correspond to ordered ZTs for rep1 followed by rep2


# rhythmicity detection scripts
classifyExpression_update.py
this was run with the following parameters:
-sig T --symbol T --minRP 0 --nonZero F --med 1 --minMax 1.4 --freqs 2 --nProcessors 24

classifyGenes_filter_update.py
this was run with the following parameters:
--a True --rt 0.05 --at 0.05


# plot scripts
used to plot histograms in Figure 2A:
overlayHists_RP24s_subplots.py

used to plot histograms in Supplementary Figure 1:
overlayHists_youngVsOldRP24s.py
