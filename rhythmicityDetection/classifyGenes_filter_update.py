#!/local/cluster/bin/python
import numpy,csv,argparse

# for classifyExpression9_filter.py
# adapted from classifyGenes4.py

###########
## INPUT ##
###########

def readGeneGroupFile(geneGroupFile,printAll):
    expressionInfo = {}
    for line in open(geneGroupFile, 'r'):
        if not line.startswith('ID'):
            geneID,symbol,RP24,qVal,pVal,SP,PVEC,phase,foldChange,median,nonZero  = line.strip().split('\t')
            detectable = 0 
            if float(median) >= 1.0 and float(foldChange) >= 1.5:
                detectable = 1
            if qVal != "NOTEST":
                expressionInfo[geneID] = (symbol,float(qVal),float(pVal),float(RP24),float(SP),float(PVEC),float(phase),float(foldChange),float(median),detectable,int(nonZero))
            else:
                expressionInfo[geneID] = (symbol,qVal,pVal,float(RP24),float(SP),float(PVEC),float(phase),float(foldChange),float(median),detectable,int(nonZero))
    return expressionInfo

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def getArgs():
    parser = argparse.ArgumentParser(description="classify groups", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--ggf1",type=str,required=True,help="classified group file 1",default=False)
    parser.add_argument("--ggf2",type=str,required=True,help="classified group file 2",default=False)
    parser.add_argument("--a",type=str2bool,required=False,help="print both detectable and non-detectable genes", nargs='?',const=True,default=False)
    parser.add_argument("--rt",type=float,required=False,help="rhythmicity threshold",default=0.05)
    parser.add_argument("--at",type=float,required=False,help="arrhythmic threshold",default=0.3)
    parser.add_argument("--spt",type=float,required=False,help="spectral parity threshold",default=2.0)
    parser.add_argument("--pvect",type=float,required=False,help="phase variance of even coefficients",default=1.0)

    args = parser.parse_args()
    return args

#################
## SUBROUTINES ##
#################

def classify(expressionInfo,rhythmicThresh,arrhythmicThresh,SPthresh,PVECthresh):
    groups = {}
    for gene in expressionInfo:
        symbol, Qval, Pval, RP24, SP, PVEC, phase, foldChange, median, detect, nonZero = expressionInfo[gene]
        if Qval != "NOTEST":
            if Qval <= rhythmicThresh:
                groups[gene] = "rhythmic"
            elif Qval > arrhythmicThresh:
                groups[gene] = "arrhythmic"
            else:
                groups[gene] = "indeterminant"
        else:
            # either high q-value or min/max foldChange is < 1.4
            groups[gene] = "arrhythmic"
        # every gene should be labeld at this point. 
        if SP >= SPthresh and PVEC <= PVECthresh and groups[gene] != "rhythmic":
            groups[gene] = "staccato"
    return groups

def identifyGeneGroups(yGroups,oGroups,yExpressionInfo,oExpressionInfo,printAll):
    R2R,R2A,R2S = {},{},{}
    A2R,A2A,A2S = {},{},{}
    S2R,S2A,S2S = {},{},{}
    for gene in yExpressionInfo:
        if gene in oExpressionInfo:
            symbol, yQval, yPval, yRP24, ySP, yPVEC, yPhase, yFold, yMedian, yDetect, yNonZero = yExpressionInfo[gene]
            symbol, oQval, oPval, oRP24, oSP, oPVEC, oPhase, oFold, oMedian, oDetect, oNonZero = oExpressionInfo[gene]
            expInfo = (symbol,yQval,oQval,yPval,oPval,yRP24,oRP24,ySP,oSP,yPVEC,oPVEC,yPhase,oPhase,yMedian,oMedian,yDetect,oDetect,yNonZero,oNonZero)
            # if printAll, then print regardless of whether detectable.
            if printAll:                
                #Rhythmic to Rhythmic case
                if yGroups[gene] == "rhythmic" and oGroups[gene] == "rhythmic":
                    R2R[gene] = expInfo                    
                #Rhythmic to Arrhythmic case
                if yGroups[gene] == "rhythmic" and oGroups[gene] == "arrhythmic":
                    R2A[gene] = expInfo
                #Rhythmic to staccato case
                if yGroups[gene] == "rhythmic" and oGroups[gene] == "staccato":
                    R2S[gene] = expInfo                        
                #Arrhythimc to Rhythmic case
                if yGroups[gene] == "arrhythmic" and oGroups[gene] == "rhythmic":
                    A2R[gene] = expInfo
                #Arrhythmic to Arrhythmic case
                if yGroups[gene] == "arrhythmic" and oGroups[gene] == "arrhythmic":
                    A2A[gene] = expInfo
                #Arrhythmic to staccato case
                if yGroups[gene] == "arrhythmic" and oGroups[gene] == "staccato":
                    A2S[gene] = expInfo
                #Staccato to Rhythmic case
                if yGroups[gene] == "staccato" and oGroups[gene] == "rhythmic":
                    S2R[gene] = expInfo
                #Staccato to Arrhythmic case
                if yGroups[gene] == "staccato" and oGroups[gene] == "arrhythmic":
                    S2A[gene] = expInfo
                #Staccato to Staccato case
                if yGroups[gene] == "staccato" and oGroups[gene] == "staccato":
                    S2S[gene] = expInfo
            else:
                #Rhythmic to Rhythmic case
                if yGroups[gene] == "rhythmic" and oGroups[gene] == "rhythmic":
                    if yDetect == 1 or oDetect == 1:
                        R2R[gene] = expInfo
                    else:
                        print(gene)
                #Rhythmic to Arrhythmic case
                elif yGroups[gene] == "rhythmic" and oGroups[gene] == "arrhythmic":
                    if yDetect == 1: 
                        R2A[gene] = expInfo
                    else:
                        print(gene)
                #Rhythmic to staccato case
                elif yGroups[gene] == "rhythmic" and oGroups[gene] == "staccato":
                    if yDetect == 1: 
                        R2S[gene] = expInfo                        
                    else:
                        print(gene)
                #Arrhythimc to Rhythmic case
                elif yGroups[gene] == "arrhythmic" and oGroups[gene] == "rhythmic":
                    if oDetect == 1:
                        A2R[gene] = expInfo
                    else: 
                        print(gene)
                #Arrhythmic to Arrhythmic case
                elif yGroups[gene] == "arrhythmic" and oGroups[gene] == "arrhythmic":
                    A2A[gene] = expInfo
                #Arrhythmic to staccato case
                elif yGroups[gene] == "arrhythmic" and oGroups[gene] == "staccato":
                    A2S[gene] = expInfo
                #Staccato to Rhythmic case
                elif yGroups[gene] == "staccato" and oGroups[gene] == "rhythmic":
                    S2R[gene] = expInfo
                #Staccato to Arrhythmic case
                elif yGroups[gene] == "staccato" and oGroups[gene] == "arrhythmic":
                    S2A[gene] = expInfo
                #Staccato to Staccato case
                elif yGroups[gene] == "staccato" and oGroups[gene] == "staccato":
                    S2S[gene] = expInfo
                else:
                    print(gene, yGroups[gene], oGroups[gene])


    return R2R,R2A,R2S,A2R,A2A,A2S,S2R,S2A,S2S

def printList(genes, groupType):
    outFile = groupType+'.txt'
    OUT = open(outFile, 'w')
    OUT.write('geneID\tsymbol\tyQvalue\toQvalue\tyPvalue\toPvalue\tyRP24\toRP24\tySP\toSP\tyPVEC\toPVEC\tyPhase\toPhase\tyMedian\toMedian\tyDetectable\toDetectable\tyNonZero\toNonzero\n')
    for gene in genes:
        symbol,yQval,oQval,yPval,oPval,yRP24,oRP24,ySP,oSP,yPVEC,oPVEC,yPhase,oPhase,yMedian,oMedian,yDetect,oDetect,yNonZero,oNonZero = genes[gene]
        OUT.write('%s\t%s\t%s\t%s\t%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\n' % (gene,symbol,str(yQval),str(oQval),str(yPval),str(oPval),yRP24,oRP24,ySP,oSP,yPVEC,oPVEC,yPhase,oPhase,yMedian,oMedian,yDetect,oDetect,yNonZero,oNonZero))


##########
## MAIN ##
##########

args = getArgs()

youngGeneGroupFile = args.ggf1
oldGeneGroupFile = args.ggf2
printAll = args.a
print(printAll)
rhythmicThresh = args.rt
arrhythmicThresh = args.at
SPthresh = args.spt
PVECthresh = args.pvect

yExpressionInfo = readGeneGroupFile(youngGeneGroupFile,printAll)
oExpressionInfo = readGeneGroupFile(oldGeneGroupFile,printAll)

yGroups = classify(yExpressionInfo,rhythmicThresh,arrhythmicThresh,SPthresh,PVECthresh)
oGroups = classify(oExpressionInfo,rhythmicThresh,arrhythmicThresh,SPthresh,PVECthresh)

R2R,R2A,R2S,A2R,A2A,A2S,S2R,S2A,S2S = identifyGeneGroups(yGroups,oGroups,yExpressionInfo,oExpressionInfo,printAll)

printList(R2R,"R2R")
printList(R2A,"R2A")
printList(R2S,"R2S")
printList(A2R,"A2R")
printList(A2A,"A2A")
printList(A2S,"A2S")
printList(S2R,"S2R")
printList(S2A,"S2A")
printList(S2S,"S2S")
