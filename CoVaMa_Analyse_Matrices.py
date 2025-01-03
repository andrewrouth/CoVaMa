##      Copyright (c) 2014, 2015 Andrew Laurence Routh
##      
##      Permission is hereby granted, free of charge, to any person obtaining a copy
##      of this software and associated documentation files (the "Software"), to deal
##      in the Software without restriction, including without limitation the rights
##      to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
##      copies of the Software, and to permit persons to whom the Software is
##      furnished to do so, subject to the following conditions:
##      
##      The above copyright notice and this permission notice shall be included in
##      all copies or substantial portions of the Software.
##      
##      THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
##      IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
##      FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
##      AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
##      LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
##      OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
##      THE SOFTWARE.
##
##      ----------------------------------------------------------------------------------------
print '\n-------------------------------------------------------------------------------------------'
print 'Co-Variation Mapper 0.5 - Analysis Module - written by Andrew Routh'
print 'Last modified 17/07/2017'
print '-------------------------------------------------------------------------------------------'
##      ----------------------------------------------------------------------------------------
import time
start = time.time()
import numpy as np
np.seterr(all='ignore')
import cPickle as pickle
import config as cfg
import argparse
from collections import defaultdict

##      -------------------------------------------------------------------------------------------
##      Take variable names from Command Line.  Variables are stored in the config.py file in the
##      same folder as the parent script.
##      -------------------------------------------------------------------------------------------
if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("PickleFile", help="Enter name of PickleFile")
	parser.add_argument("Output_File", help="Enter name of desired output file")
	parser.add_argument("Mode", help="Enter Matrix Mode, either nucleotides: NT or Amino Acids: AA")
	parser.add_argument("--Min_Coverage", help="Enter minimum coverage over pairs of associated nucleotides. Default value is 1000")
	parser.add_argument("--Min_Fusion_Coverage", help="Enter minimum coverage over recombination event. Default value is 1000")
	parser.add_argument("--Min_Pop_Fraction", help="Enter a float between 0 and 1 for the required association fraction. Default value is 0.001.")
	parser.add_argument("-OutArray",  action='store_true', help="Write array values into output file if LD is found. Default value is false.")
	parser.add_argument("-OutAllLDs",  action='store_true', help="Write out all LD values found. Default value is output highest value only.")
	parser.add_argument("-Merge", action='store_true', help="If input Matrix is from merged matrices, select -Merge. Default is False")
	parser.add_argument("-Weighted", action='store_true', help="Weight each LD value according to fraction of reads from total contingency table. Default is False")
	parser.add_argument("--Min_Merged_Matrices", help="Enter a minimum number of matrices required that were used to generated merged matrix. Default value is 1.")
	parser.add_argument("--LD_Precision", help="Enter number of floating point digits for LD or R2 output. Default value is 6.")
	args = parser.parse_args()

        ## Required Files
	cfg.Output_File = str(args.Output_File)
	cfg.PickleFile = str(args.PickleFile)
	cfg.Mode = str(args.Mode)
	if cfg.Mode != 'AA' and cfg.Mode != 'NT':
                print "ERROR! Incorrect mode selected. Must be either 'AA' or 'NT'. Defaulting to NT."
                cfg.Mode = 'NT'
        else:
                pass
	
	##Parameters
	if args.Min_Coverage:
		cfg.Min_Coverage = int(args.Min_Coverage)
	else:
		cfg.Min_Coverage = 1000
	if args.Min_Fusion_Coverage:
		cfg.Min_Fusion_Coverage = int(args.Min_Fusion_Coverage)
	else:
		cfg.Min_Fusion_Coverage = cfg.Min_Coverage
	if args.Min_Pop_Fraction:
		cfg.Min_Pop_Fraction = float(args.Min_Pop_Fraction)
	else:
		cfg.Min_Pop_Fraction = 0.001
	if args.LD_Precision:
		cfg.LD_Precision = int(args.LD_Precision)
	else:
		cfg.LD_Precision = 6
	if args.OutArray:
                cfg.OutArray = True
        else:
                cfg.OutArray = False
        if args.OutAllLDs:
                cfg.OutAllLDs = True
        else:
                cfg.OutAllLDs = False
        if args.Merge:
                cfg.Merge = True
        else:
                cfg.Merge = False
        if args.Weighted:
                cfg.Weighted = True
        else:
                cfg.Weighted = False
	if args.Min_Merged_Matrices:
		cfg.Min_Merged_Matrices = int(args.Min_Merged_Matrices)
	else:
		cfg.Min_Merged_Matrices = 1
##      -------------------------------------------------------------------------------------------
##      The Basics 
##      -------------------------------------------------------------------------------------------

Codons = {'ACC': 'T', 'GTC': 'V', 'ACA': 'T', 'ACG': 'T',
          'GTT': 'V', 'AAC': 'N', 'CCT': 'P', 'TGG': 'W',
          'AGC': 'S', 'ATC': 'I', 'CAT': 'H', 'AAT': 'N',
          'AGT': 'S', 'ACT': 'T', 'GTG': 'V', 'CAC': 'H',
          'AAA': 'K', 'CCG': 'P', 'CCA': 'P', 'CAA': 'Q',
          'CCC': 'P', 'GGT': 'G', 'TCT': 'S', 'GCG': 'A',
          'TGC': 'C', 'CAG': 'Q', 'TGA': 'STOP', 'TAT': 'Y',
          'CGG': 'R', 'TCG': 'S', 'AGG': 'R', 'GGG': 'G',
          'TCC': 'S', 'TCA': 'S', 'GAA': 'E', 'TAA': 'STOP',
          'GGA': 'G', 'TAC': 'Y', 'CGT': 'R', 'TAG': 'STOP',
          'ATA': 'I', 'GCA': 'A', 'CTT': 'L', 'GGC': 'G',
          'ATG': 'M', 'CTG': 'L', 'GAG': 'E', 'CTC': 'L',
          'AGA': 'R', 'CTA': 'L', 'GCC': 'A', 'AAG': 'K',
          'GAT': 'D', 'TTT': 'F', 'GAC': 'D', 'GTA': 'V',
          'CGA': 'R', 'GCT': 'A', 'TGT': 'C', 'ATT': 'I',
          'TTG': 'L', 'TTA': 'L', 'CGC': 'R', 'TTC': 'F'}

Bases = ['A','T','G','C']
AAIndex = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

AACoords = [(i,j) for i in range(20) for j in range(20) if i<j]
AAAllCoords = [(i+j) for i in AACoords for j in AACoords]
NucCoords = [(i,j) for i in range(4) for j in range(4) if i<j]
AllCoords = [(i+j) for i in NucCoords for j in NucCoords]

def dd400():
        return [0]*400
def dd40():
        return [0]*40
def dd16():	
	return [0]*16
def dd8():	
	return [0]*8
def dd4():	
	return [0]*4
def dd2():	
	return [0]*2

##      -------------------------------------------------------------------------------------------
##      
##      -------------------------------------------------------------------------------------------
##@profile

def Linkage(Matrix, ArrayTot, Total):
        ArrayTot = float(ArrayTot)
        pA, pa = np.sum(Matrix, axis = 1)/ArrayTot
        pB, pb = np.sum(Matrix, axis = 0)/ArrayTot
	pAB = Matrix[0,0]/ArrayTot
	pAb = Matrix[0,1]/ArrayTot
	paB = Matrix[1,0]/ArrayTot
	pab = Matrix[1,1]/ArrayTot
	
	## Linkage Disequilibrium Calc
	LD = (pAB*pab)-(pAb*paB)

        ## R2 Calc
        RSquare = (LD*LD)/(pA * pa * pB * pb)

        ## LDmax Calc
        if LD <= 0:
                LDmax = min([(pA*pB), (pa*pb)])
        else:
                LDmax = min([(pA*pb), (pA*pb)])
        
        ## Optional Weighting
        if cfg.Weighted:
                LD = LD*(ArrayTot/Total)
                RSquare = RSquare*(ArrayTot/Total)
                LDmax = LDmax*(ArrayTot/Total)
        else:
                pass
        
        LD = round(LD, cfg.LD_Precision)
        RSquare = round(RSquare, cfg.LD_Precision)
        LDmax = round(LDmax, cfg.LD_Precision)
	return LD, RSquare, LDmax


##      -------------------------------------------------------------------------------------------
##      
##      -------------------------------------------------------------------------------------------

def AnyLessThan(List, Value):
        for i in List:
                if i < Value:
                        return True
        return False

##@profile
def AllLDs(Matrix, Type, Total):
        LDs = []
        RSquares = []
        LDMaxes = []
        Total = float(Total)
        if cfg.Mode == 'NT':
                TypeCoords = NucCoords
                Elements = Bases
        else:
                TypeCoords = AACoords
                Elements = AAIndex
        if Type == 'Nucs':
                OutCoords = []
                for i in AllCoords:
                        ##Extract Flat SubArray from Matrix
                        ## Coord in format (e.g.): (0,1,1,2)
                        ## Coord[0] and Coord[1] give y values (5' Nuc) e.g. A and T
                        ## Coord[2] and Coord[3] give x values (3' Nuc) e.g. T and G
                        Array = Matrix[[i[0],i[0],i[1],i[1]],[i[2],i[3],i[2],i[3]]]
                        ArrayTot = Array[0]+Array[1]+Array[2]+Array[3]
                        if ArrayTot > cfg.Min_Coverage:
                                Sums = [Array[0]+Array[1], Array[2]+Array[3], Array[0]+Array[2], Array[1]+Array[3]]
                                Fractions = [x/Total for x in Sums]
                                ##Test that each allele is sufficiently populated in table
                                if not AnyLessThan(Sums, cfg.Min_Coverage):
                                        ##Test that each allele is sufficiency frequent in Total
                                        ##Prevents acceptance of highly populated table, but with rare mutation
                                        if not AnyLessThan(Fractions, cfg.Min_Pop_Fraction):
                                                ##Reshape flat array into numpy array after tests complete
                                                Array = Array.reshape(2,2)
                                                LD, RSquare, LDMax = Linkage(Array, ArrayTot, Total)
                                                LDs.append(abs(LD))
                                                RSquares.append(RSquare)
                                                LDMaxes.append(abs(LDMax))
                                                if LD > 0:
                                                        k = (Elements[i[0]], Elements[i[2]], Elements[i[1]], Elements[i[3]])
                                                        OutCoords.append(k)
                                                else:
                                                        k = (Elements[i[0]], Elements[i[3]], Elements[i[1]], Elements[i[2]])
                                                        OutCoords.append(k)
                                        else:
                                                pass
                                else:
                                        pass
                        else:
                                pass
        elif Type == 'Mix':
                OutCoords = []
                for i in TypeCoords:
                        Array = Matrix[[i[0],i[0],i[1],i[1]],[0,1,0,1]]
                        ## Coord in format (e.g.): (0,2)
                        ## Coord[0] and Coord[1] give y values (5' Nuc) e.g. A and G
                        ## Coord[0] and Coord[1] give y values (5' AA) e.g. A and D
                        ArrayTot = Array[0]+Array[1]+Array[2]+Array[3]
                        if ArrayTot > cfg.Min_Fusion_Coverage:
                                Sums = [Array[0]+Array[1], Array[2]+Array[3], Array[0]+Array[2], Array[1]+Array[3]]
                                Fractions = [x/Total for x in Sums]
                                ##Test that each allele is sufficiently populated in table
                                if not AnyLessThan(Sums, cfg.Min_Fusion_Coverage):
                                        ##Test that each allele is sufficiency frequent in Total
                                        ##Prevents acceptance of highly populated table, but with rare mutation
                                        if not AnyLessThan(Fractions, cfg.Min_Pop_Fraction):
                                                Array = Array.reshape(2,2)
                                                LD, RSquare, LDMax = Linkage(Array, ArrayTot, Total)
                                                LDs.append(abs(LD))
                                                RSquares.append(RSquare)
                                                LDMaxes.append(abs(LDMax))
                                                if LD > 0:
                                                        k = (Elements[i[0]], 'Neg', Elements[i[1]], 'Pos')
                                                        OutCoords.append(k)
                                                else:
                                                        k = (Elements[i[0]], 'Pos', Elements[i[1]], 'Neg')
                                                        OutCoords.append(k)
                                        else:
                                                pass
                                else:
                                        pass
                        else:
                                pass
        else:
                ##Rec-to-Rec array
                x1 = np.sum(Matrix,axis=0)
                x2 = np.sum(Matrix,axis=1)
                xf1 = [x/Total for x in x1]
                xf2 = [x/Total for x in x2]
                OutCoords = ['Rec-Rec']
                if np.all(x1 > cfg.Min_Fusion_Coverage) and np.all(x2 > cfg.Min_Fusion_Coverage):
                        if np.all(xf1 > cfg.Min_Pop_Fraction) and np.all(xf2 > cfg.Min_Pop_Fraction):
                                LD, RSquare, LDMax = Linkage(Matrix, Total, Total)
                                LDs.append(abs(LD))
                                RSquares.append(RSquare)
                                LDMaxes.append(abs(LDMax))
                                if LD > 0:
                                        k = ('Neg', 'Neg', 'Pos', 'Pos')
                                        OutCoords.append(k)
                                else:
                                        k = ('Neg', 'Pos', 'Pos', 'Neg')
                                        OutCoords.append(k)
                        else:
                                pass
                else:
                        pass
        return LDs, OutCoords, RSquares, LDMaxes
                
##      ------------------------------------------------------------------------------------------
##      Analyse_Popultation_Matrix() will assess every covariation matrix for evidence of non-independent 
##      variation.
##      Matrices exhibiting non-independence are written into the output file specified in the cmdline.
##      -------------------------------------------------------------------------------------------
##@profile
def Analyse_Population_Matrix(Libs, i, j):
        global Output
        global ContingencyTableCount
        global CurrentTable
        global StatsLD
        global StatsR2
        global StatsLDMax
        CurrentTable += 1
##        if CurrentTable % ContingencyTable5 == 0:
##                                print " .... %s Percent of Contingency Tables Processed .... " % (CurrentTable/ContingencyTable5*5)
##        else:
##                pass
        if '_to_Rec' in Libs:
                if Libs == 'Rec_to_Rec':
                        Type ="Fusions"
                        k = Total_Matrices[Libs][i][j]
                        if cfg.Merge:
                                ## k[1] Gives number of contributors to merged matrix
                                NoMerged = k[1]
                                if NoMerged >= cfg.Min_Merged_Matrices:
                                        k = np.array(k[0], dtype = int).reshape(2,2)
                                        MergeOK = True
                                else:
                                        MergeOK = False
                        else:
                                k = np.array(k, dtype = int).reshape(2,2)
                                MergeOK = True
                else:
                        Type = "Mix"
                        k = Total_Matrices[Libs][i][j]
                        if cfg.Merge:
                                ## k[1] Gives number of contributors to merged matrix
                                NoMerged = k[1]
                                if NoMerged >= cfg.Min_Merged_Matrices:
                                        if cfg.Mode == 'NT':
                                                k = np.array(k[0], dtype = int).reshape(4,2)
                                        else:
                                                k = np.array(k[0], dtype = int).reshape(20,2)
                                        MergeOK = True
                                else:
                                        MergeOK = False
                        else:
                                if cfg.Mode == 'NT':
                                        k = np.array(k, dtype = int).reshape(4,2)
                                else:
                                        k = np.array(k, dtype = int).reshape(20,2)
                                MergeOK = True
        else:
                Type = "Nucs"
                k = Total_Matrices[Libs][i][j]
                if cfg.Merge:
                        ## k[1] Gives number of contributors to merged matrix
                        if k[1] >= cfg.Min_Merged_Matrices:
                                if cfg.Mode == 'NT':
                                        k = np.array(k[0], dtype = int).reshape(4,4)
                                else:
                                        k = np.array(k[0], dtype = int).reshape(20,20)
                                MergeOK = True
                        else:
                                MergeOK = False
                else:
                        if cfg.Mode == 'NT':
                                k = np.array(k, dtype = int).reshape(4,4)
                        else:
                                k = np.array(k, dtype = int).reshape(20,20)
                        MergeOK = True
        if MergeOK:
                kTotal = np.sum(k)
                ## Min_Fusion_Coverage cannot be greater than Min_Coverage and therefore is more permissive
                if kTotal >= cfg.Min_Fusion_Coverage:
                        ##FindLD
                        LDs, OutCoords, RSquares, LDMax = (AllLDs(k, Type, kTotal))
                        if Type =='Nucs':
                                [StatsLD.append(x) for x in LDs]
                                [StatsR2.append(x) for x in RSquares]
                                [StatsLDMax.append(x) for x in LDMax]
                        else:
                                pass
                        if LDs:
                                if cfg.OutArray:
                                        ##Flatten Array for output
                                        Array = str([m for m in k.ravel()])
                                else:
                                        Array = ''
                                if cfg.OutAllLDs:
                                        Output = (str(Libs) + "\t" + str(i) + "\t" + str(j) + "\t" + str(LDs) + "\t" + str(RSquares) + "\t" + str(OutCoords) + "\t" + Array + "\n")
                                        Out.write(''.join(Output))
                                else:   
                                        MaxLD = max(LDs)
                                        Index = LDs.index(MaxLD)
                                        MaxRSquares = RSquares[Index]
                                        MaxLDMax = LDMax[Index]
                                        MaxCoord = OutCoords[Index]
                                        Output = (str(Libs) + "\t" + str(i) + "\t" + str(j) + "\t" + str(MaxLD) + "\t" + str(MaxRSquares) + "\t" + str(MaxLDMax) + "\t" + str(MaxCoord) + "\t" + Array + "\n")
                                        Out.write(''.join(Output))
                        else:
                                pass
                else:
                        pass
        else:
                pass
##      -------------------------------------------------------------------------------------------
##      Run Main Script
##      -------------------------------------------------------------------------------------------

if __name__ == "__main__":
        ##Analyse Association Matrix. 
	print 'Reading in Pickle file...'
	Total_Matrices = pickle.load(open(cfg.PickleFile,'rb'))
	print 'Analysing Populated Association Matrix'
	Out = open(cfg.Output_File, 'w')
	ContingencyTableCount = 0
	StatsLD = []
	StatsR2 = []
	StatsLDMax = []
	for Libs in Total_Matrices:
                for i in Total_Matrices[Libs]:
                        for j in Total_Matrices[Libs][i]:
                                ContingencyTableCount += 1
        ContingencyTable5 = ContingencyTableCount/20
        CurrentTable = 0
        print '%s Contingency Tables to be processed...' % (ContingencyTableCount)
	[Analyse_Population_Matrix(Libs, i, j) for Libs in Total_Matrices for i in Total_Matrices[Libs] for j in Total_Matrices[Libs][i]]
	print 'Saving results to %s' % (cfg.Output_File)

	##Calculate Final LDs Stats:
	NumValues = len(StatsLD)
        SDLD = np.std(StatsLD)
        MeanLD = np.mean(StatsLD)
        print "For LDs: Std and mean are: ", SDLD, MeanLD
        if NumValues >= 370:
                ## Calculate 3 sigma Stats
                Sigma3LD = 3*SDLD + MeanLD
                ExpectedSigma3 = NumValues/370
                RealSigma3 = sum(1 for i in StatsLD if i > Sigma3LD)
                FDR3 = (ExpectedSigma3/float(RealSigma3)) * 100
                print "%s events found with greater than 3sigma, with a false discovery rate of: %s percent" % (RealSigma3, FDR3)
        else:
                print "Too few events found to determine FDR with 3 sigma."
        if NumValues >= 15787:
                ## Calculate 4 sigma Stats
                Sigma4LD = 4*SDLD + MeanLD
                ExpectedSigma4 = NumValues/15787
                RealSigma4 = sum(1 for i in StatsLD if i > Sigma4LD)
                FDR4 = (ExpectedSigma4/float(RealSigma4)) * 100
                print "%s events found with greater than 4sigma, with a false discovery rate of: %s percent" % (RealSigma4, FDR4)
        else:
                pass
        if NumValues >= 1744278:
                ## Calculate 5 sigma Stats
                Sigma5LD = 5*SDLD + MeanLD
                ExpectedSigma5 = NumValues/1744278
                RealSigma5 = sum(1 for i in StatsLD if i > Sigma5LD)
                FDR5 = (ExpectedSigma5/float(RealSigma5)) * 100
                print "%s events found with greater than 5sigma, with a false discovery rate of: %s percent" % (RealSigma5, FDR5)
        else:
                pass
##        SDR2 = np.std(StatsR2)
##        MeanR2 = np.mean(StatsR2)
##        Sigma3R2 = 3*SDR2 + MeanR2
##        print "For RSquares: Std, mean, and three sigma are: ", SDR2, MeanR2, Sigma3R2
	Out.close()

finish = time.time()
print "Time to complete in seconds: ", float(finish - start)

