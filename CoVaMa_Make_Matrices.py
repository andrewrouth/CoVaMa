##      Copyright (c) 2014-2025 Andrew Laurence Routh
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
print('\n-------------------------------------------------------------------------------------------')
print('Co-Variation Mapper 0.8 - written by Andrew Routh')
print('Last modified 3/01/2025')
print('-------------------------------------------------------------------------------------------')
##      ----------------------------------------------------------------------------------------
import time
start = time.time()
import cPickle as pickle
import gzip
import config as cfg
import argparse
try:
        import pysam
except:
        print("Pysam not installed; no support for .BAM files.")
from re import findall
from collections import defaultdict

##Testi = 'NC_004144_FHV_RNA2.seq_to_NC_004144_FHV_RNA2.seq_@_1349_AA_1349_@_Ins'
##Testj = 'NC_004144_FHV_RNA2.seq_to_NC_004144_FHV_RNA2.seq_@_1384_to_1387_@_Rec'

##      -------------------------------------------------------------------------------------------
##      Take variable names from Command Line.  Variables are stored in the config.py file in the
##      same folder as the parent script.
##      -------------------------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("Output_Tag", help="Enter name of desired output tag: e.g. 'FHV' gives FHV.Total_Matrices.py.pi")
    parser.add_argument("FASTAFile", help="Enter name of Fasta reference genome used in original sequence alignment")
    parser.add_argument("Mode", help="Enter Matrix Mode, either nucleotides: NT or Amino Acids: AA")
    parser.add_argument("--Mode2", help="Enter variant types: 'Nucs', 'Recs' or 'Both'. Default == Both")
    parser.add_argument("--SAM1", help="Enter name of SamFile")
    parser.add_argument("--SAM2", help="Enter name of second/paired SamFile")
    parser.add_argument("--ViReMa_Output1", help="Enter location of ViReMa 'CoVaMa_Output.txt' file")
    parser.add_argument("--ViReMa_Output2", help="Enter location of paired ViReMa 'CoVaMa_Output.txt' file")
    parser.add_argument("--PileUp_Fraction", help="Only create Matrices for mutants present above this frequency. Default == 0.0")
    parser.add_argument("--Min_Coverage_Output", help="Number of mapped nucleotides per contigency table required for output. Default == 100")
    parser.add_argument("--MicroInDel_Length", help="Required number of nucleotides to constitute MicroInDel event. Default == 5")
    parser.add_argument("--InDel_Exclusion", help="Required number of nucleotides to exclude MicroInDel event. Default == 5")
    parser.add_argument("--Rec_Exclusion", help="Required number of nucleotides to exclude recombination event. Default == 10")
    parser.add_argument("--Min_Fusion_Coverage", help="Enter minimum coverage over pairs of associated nucleotides. Default value is 100")
    parser.add_argument("--NtStart", help="Enter start range for pairs of associated nucleotides.")
    parser.add_argument("--NtFinish", help="Enter start range for pairs of associated nucleotides.")
    parser.add_argument("--ORFStartNuc", help="Enter start range for pairs of associated AAs.")
    parser.add_argument("--ORFFinishNuc", help="Enter end range for pairs of associated AAs.")
    parser.add_argument("--Ends", help="Enter number of nucleotides to ignore from 5' and 3' extremities (Soft-padding). Default value is 0")
    parser.add_argument("-PileUpOnly", action='store_true', help="Only calculate Pileup tables then quit")
    parser.add_argument("--PrintReadCount", help="Print read count during read-in (helps debug or approximate run-time). Enter chunk amount")

    args = parser.parse_args()

    ##File Handling
    cfg.Output_Tag = str(args.Output_Tag)
    cfg.FASTAFile = str(args.FASTAFile)
    cfg.Mode = str(args.Mode)
    if cfg.Mode != 'AA' and cfg.Mode != 'NT':
        print("ERROR! Incorrect mode selected. Must be either 'AA' or 'NT'. Defaulting to NT.")
        cfg.Mode = 'NT'
    else:
        pass

    if args.Mode2:
                cfg.Mode2 = str(args.Mode2)
    else:
                cfg.Mode2 = 'Both'
    if cfg.Mode2 != 'Nucs' and cfg.Mode2 != 'Recs' and cfg.Mode2 != 'Both':
        print("ERROR! Incorrect variant mode selected. Must be either 'Nucs', 'Recs' or 'Both'. Defaulting to Both.")
        cfg.Mode2 = 'Both'
    else:
        pass

    if args.SAM1:
        cfg.SAMIN1 = str(args.SAM1)
    else:
        cfg.SAMIN1 = None
    if args.SAM2:
        cfg.SAMIN2 = str(args.SAM2)
    else:
        cfg.SAMIN2 = None
    if args.ViReMa_Output1:
        cfg.ViReMaIn1 = str(args.ViReMa_Output1)
    else:
        cfg.ViReMaIn1 = None
    if args.ViReMa_Output2:
        cfg.ViReMaIn2 = str(args.ViReMa_Output2)
    else:
        cfg.ViReMaIn2 = None
    if args.PileUpOnly:
        cfg.PileUpOnly = True
    else:
        cfg.PileUpOnly = False
    if args.PrintReadCount:
        cfg.PrintReadCount = int(args.PrintReadCount)
    else:
        cfg.PrintReadCount = 0

    ##Parameters
    if args.PileUp_Fraction:
        cfg.PileUp_Fraction = float(args.PileUp_Fraction)
    else:
        cfg.PileUp_Fraction = 0.0
    if args.Min_Coverage_Output:
        cfg.Min_Coverage_Output = int(args.Min_Coverage_Output)
    else:
        cfg.Min_Coverage_Output = 100
    if args.Ends:
        cfg.Ends = int(args.Ends)
    else:
        cfg.Ends = 0

    ##Mapping ranges
    if args.ORFStartNuc:
        cfg.ORFStartNuc = int(args.ORFStartNuc)
    else:
        cfg.ORFStartNuc = False
    if args.ORFFinishNuc:
        cfg.ORFFinishNuc = int(args.ORFFinishNuc)
        cfg.ProtLength = (cfg.ORFFinishNuc - cfg.ORFStartNuc + 1)/3
    else:
        cfg.ORFFinishNuc = False

    if args.NtStart:
        cfg.NtStart = int(args.NtStart)
    else:
        cfg.NtStart = 0
    if args.NtFinish:
        cfg.NtFinish = int(args.NtFinish)
    else:
        cfg.NtFinish = 1000000000

    ###Fusion Variables
    if args.MicroInDel_Length:
        cfg.uInDel = int(args.MicroInDel_Length)
    else:
        cfg.uInDel = 5
    if args.InDel_Exclusion:
        cfg.XInDel = int(args.InDel_Exclusion)
    else:
        cfg.XInDel = 5
    if args.Rec_Exclusion:
        cfg.XRec = int(args.Rec_Exclusion)
    else:
        cfg.XRec = 10
    if args.Min_Fusion_Coverage:
        cfg.Min_Fusion_Coverage = int(args.Min_Fusion_Coverage)
    else:
        cfg.Min_Fusion_Coverage = 100

##      -------------------------------------------------------------------------------------------
##      ReadFasta() takes a fasta file and reads in the gene entries.
##      -------------------------------------------------------------------------------------------

def ReadFasta():
    with open(cfg.FASTAFile,'r') as Input_Genome:
        Name = Input_Genome.readline().rstrip()[1:]
        Name = Name.split()[0]
        cfg.Genes_Lib[Name] = ''
        print(Name)
        Line = Input_Genome.readline().rstrip()
        while Line:
            if '>' not in Line:
                cfg.Genes_Lib[Name] += Line.upper()
            else:
                Name = Line[1:]
                Name = Name.split()[0]
                cfg.Genes_Lib[Name] = ''
                print(Name)
            Line = Input_Genome.readline().rstrip()

##      -------------------------------------------------------------------------------------------
##        ##The Basics
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

AAIndex = ['A','C','D','E','F',
       'G','H','I','K','L',
       'M','N','P','Q','R',
       'S','T','V','W','Y']

##      -------------------------------------------------------------------------------------------
##      ddn() functions allow pickling of defaultDicts.
##      -------------------------------------------------------------------------------------------

def dd400():
    return [0]*400
def dd40():
    return [0]*40
def dd20():
    return [0]*20
def dd16():
    return [0]*16
def dd8():
    return [0]*8
def dd4():
    return [0]*4
def dd2():
    return [0]*2

##      -------------------------------------------------------------------------------------------
##      Make_Gene_Dicts Generates an Association Matrix for each pair of Genes (in fact a list of lists)
##      -------------------------------------------------------------------------------------------

def Make_Gene_Dicts(Ref):
        if len(cfg.Genes_Lib[Ref]) > cfg.NtFinish:
                x = cfg.NtFinish + 1
        else:
                x = len(cfg.Genes_Lib[Ref]) + 1  #must have +1 as .SAM alignments are in leftmost position = 1.
        if cfg.Mode == 'NT':
            #Generate gene matrix List of Lists for NT
            cfg.Total_Dicts[Ref + '_to_' + Ref] = {i:defaultdict(dd16) for i in range(cfg.NtStart, x)}
            cfg.PileUp[Ref] = defaultdict(dd4)
            cfg.Filtered_Dicts[Ref + '_to_' + Ref] = {i:{} for i in range(cfg.NtStart, x)}
            cfg.Coverage_Dicts[Ref + '_to_' + Ref] = {i:defaultdict(int) for i in range(cfg.NtStart, x)}
        else:
            #Generate gene matrix List of Lists for AA
            x = cfg.ProtLength + 1 #must have +1 as .SAM alignments are in leftmost position = 1.
            cfg.Total_Dicts[Ref + '_to_' + Ref] = {i:defaultdict(dd400) for i in range(cfg.NtStart, x)}
            cfg.PileUp[Ref] = defaultdict(dd20)
            cfg.Filtered_Dicts[Ref + '_to_' + Ref] = {i:{} for i in range(cfg.NtStart, x)}
            cfg.Coverage_Dicts[Ref + '_to_' + Ref] = {i:defaultdict(int) for i in range(cfg.NtStart, x)}
        cfg.AllowedCoords[Ref] = set()

##      -------------------------------------------------------------------------------------------
##      Make_Fusion_Dicts Generates an Association Matrix for each pair of Genes (in fact a list of lists)
##      -------------------------------------------------------------------------------------------

def Make_Fusion_Dicts(Ref):
    #Generate gene matrix List of Lists
        if len(cfg.Genes_Lib[Ref]) > cfg.NtFinish:
                x = cfg.NtFinish + 1
        else:
                x = len(cfg.Genes_Lib[Ref]) + 1  #must have +1 as .SAM alignments are in leftmost position = 1.
        if cfg.Mode == 'NT':
            cfg.Total_Dicts[Ref + '_to_Rec'] = {i:defaultdict(dd8) for i in range(cfg.NtStart, x)}
            cfg.Filtered_Dicts[Ref + '_to_Rec'] = {i:defaultdict(dd8) for i in range(cfg.NtStart, x)}
            cfg.Coverage_Dicts[Ref + '_to_Rec'] = {i:defaultdict(dd8) for i in range(cfg.NtStart, x)}
        else:
            x = cfg.ProtLength + 1 #must have +1 as .SAM alignments are in leftmost position = 1.
            cfg.Total_Dicts[Ref + '_to_Rec'] = {i:defaultdict(dd40) for i in range(cfg.NtStart, x)}
            cfg.Filtered_Dicts[Ref + '_to_Rec'] = {i:defaultdict(dd40) for i in range(cfg.NtStart, x)}
            cfg.Coverage_Dicts[Ref + '_to_Rec'] = {i:defaultdict(dd40) for i in range(cfg.NtStart, x)}

##      -------------------------------------------------------------------------------------------
##      Yield successive n-sized chunks from l
##      -------------------------------------------------------------------------------------------

def Chunks(l, n):
    #l is list
    #n is chunksize.
    for i in range(0, len(l), n):
        yield l[i:i + n]

##      -------------------------------------------------------------------------------------------
##      Extract_Data_fromSAM() builds a dictionary of read names where each entry describes the mapping or mismatch events.
##
##      For example:
##      Name_Events['D3NJ6HQ1:256:C0V23ACXX:3:1109:12501:23122'] =
##      ['gi|1906382|gb|K03455.1|HIVHXB2CG_1478_to_1543', 'gi|1906382|gb|K03455.1|HIVHXB2CG_1544_M_A', 'gi|1906382|gb|K03455.1|HIVHXB2CG_1545_to_1560']
##
##      This tells us that read 'D3NJ6HQ1:256:C0V23ACXX:3:1109:12501:23122' mapped
##      to gene gi|1906382|gb|K03455.1|HIVHXB2CG from 1478 to 1543,
##      then there was a mismatched A, then from 1545 to 1560
##      -------------------------------------------------------------------------------------------

def Extract_Data_fromSAM(SAMFILE):
    global Name_Events
    global FusionLibCount
    SamHeaders = ['@HD', '@SQ', '@RG', '@PG', '@CO']
    InDels = ['D','N','I']
    Pads = ['S', 'H']
    if SAMFILE[-6:] == 'sam.gz':
                OPEN = gzip.open(SAMFILE,'rb')
    elif SAMFILE[-4:] == '.bam':
                OPEN = pysam.AlignmentFile(SAMFILE,'rb')
    else:
                OPEN = open(SAMFILE,'r')
    with OPEN as InputSAM:
        print("Reading Reads from SAMFile")
        Count = 0
        for line in InputSAM:
            if line[:3] not in SamHeaders:
                line = line.split('\t')
                Count +=1
                if cfg.PrintReadCount:
                    if Count%cfg.PrintReadCount == 0:
                        print(Count)
                    else:
                        pass
                else:
                    pass
                if line[5] != '*':
                    ReadName = line[0]
                    Ref = line[2]
                    Position = line[3]
                    if int(Position) < cfg.NtFinish:
                                                Read = line[9]
                                                Cigar = list(Chunks(findall(r"[^\W\d_]+|\d+", line[5]), 2))
                                                #Check padded and correct
                                                if Cigar[0][1] in Pads:
                                                        #Remove 5' pad:
                                                        Read = Read[int(Cigar[0][0]):]
                                                        ##Don't correct position, as position is reported from first mapped Nt, not including Pad
        ##                        ##Position = str(int(Position) + int(Cigar[0][0]))
                                                        Cigar = Cigar[1:]
                                                else:
                                                        #No 5' Pad
                                                        pass
                                                if Cigar[-1][1] in Pads:
                                                        #Remove 3' pad
                                                        Read = Read[:-int(Cigar[-1][0])]
                                                        Cigar = Cigar[:-1]
                                                else:
                                                        #No 3' Pad
                                                        pass

                                                #Check InDels
                                                if any([i in line[5] for i in InDels]):
                                                        #iterate by segment and turn InDels into fusion events.
                                                        CurrPosition = int(Position)
##                                                        while CurrPosition < cfg.NtFinish:
                                                        for i in Cigar:
                                                                        if i[1] == 'D' or i[1] == 'N':
                                                                                #"D" or "N" is a deletion: e.g. 'NC_004146.1_FHV_RNA1.seq_to_NC_004146.1_FHV_RNA1.seq_@_310_to_959_@_Rec'
                                                                                Positions = str(CurrPosition - 1) + '_to_' + str(CurrPosition + int(i[0]))
                                                                                if int(i[0]) <= cfg.uInDel:
                                                                                        Event = Ref + '_to_' + Ref  + "_@_" + Positions + "_@_uDel"
                                                                                else:
                                                                                        Event = Ref + '_to_' + Ref  + "_@_" + Positions + "_@_Rec"
                                                                                CurrPosition += int(i[0])
                                                                                if cfg.Mode2 != 'Nucs':
                                                                                        if Event not in FusionLibCount:
                                                                                                FusionLibCount[Event][1] = 1
                                                                                        else:
                                                                                                FusionLibCount[Event][1] += 1
                                                                                else:
                                                                                        pass
                                                                                #Don't change Read
                                                                        elif i[1] == 'I':
                                                                                #"I" is an insertion e.g. 'NC_004146.1_FHV_RNA1.seq_to_NC_004146.1_FHV_RNA1.seq_@_310_ACAGCATCACG_310_@_Rec'
                                                                                Insertion = Read[:int(i[0])]
                                                                                Positions = str(CurrPosition) + '_' + Insertion + '_' + str(CurrPosition)
                                                                                if len(Insertion) <= cfg.uInDel:
                                                                                        Event = Ref + '_to_' + Ref  + "_@_" + Positions + "_@_uIns"
                                                                                else:
                                                                                        Event = Ref + '_to_' + Ref  + "_@_" + Positions + "_@_Ins"
                                                                                Read = Read[int(i[0]):]
                                                                                if cfg.Mode2 != 'Nucs':
                                                                                        if Event not in FusionLibCount:
                                                                                                FusionLibCount[Event][1] = 1
                                                                                        else:
                                                                                                FusionLibCount[Event][1] += 1
                                                                                else:
                                                                                        pass
                                                                                #Don't change Ref Position
                                                                        elif i[1] == 'M':
                                                                                #"M"
                                                                                Segment = Read[:int(i[0])]
                                                                                if cfg.Ends > 0:
                                                                                        Segment = Segment[cfg.Ends:-cfg.Ends]
                                                                                        CurrPosition = CurrPosition + cfg.Ends
                                                                                else:
                                                                                        pass
                                                                                Event = Ref + "_@_" + str(CurrPosition) + "_@_" + Segment
                                                                                # e.g. 'NC_004146.1_FHV_RNA1.seq_@_284_@_CCGTGACGCCTAAGCGTGTCGCAGACAA'
                                                                                Read = Read[int(i[0]):]
                                                                                CurrPosition += int(i[0])
                                                                        else:
                                                                                print("Unknown CIGAR entry: ", i)
                                                                                break
                                                                        if CurrPosition > cfg.NtStart:
                                                                                if ReadName in Name_Events:
                                                                                        Name_Events[ReadName].append(Event)
                                                                                else:
                                                                                        Name_Events[ReadName] = [Event]
                                                                        else:
                                                                                pass
                                                else:
                                                        if cfg.Ends > 0:
                                                                Read = Read[cfg.Ends:-cfg.Ends]
                                                                Position = str(int(Position) + cfg.Ends)
                                                        else:
                                                                pass
                                                        if int(Position) + len(Read) > cfg.NtStart:
                                                                Event = Ref + "_@_" + Position + "_@_" + Read
                                                                if ReadName in Name_Events:
                                                                        Name_Events[ReadName].append(Event)
                                                                else:
                                                                        Name_Events[ReadName] = [Event]
        print("A total of %s reads in SAMFile were read in." %(Count))

##      -------------------------------------------------------------------------------------------
##      Take arguments from command line, and send them to the config file for cross-module access
##      -------------------------------------------------------------------------------------------

def Extract_Data_fromViReMa_Output(ViReMaFile):
        ##NEED TO ADJUST FOR NtStart / NtFinish
    global Name_Events
    global FusionLibCount
    if ViReMaFile[-3:] == '.gz':
                OPEN = gzip.open(ViReMaFile,'rb')
    else:
                OPEN = open(ViReMaFile,'r')
    with OPEN as Input:
        print("Making FusionLib")
        Data = Input.readline().split()
        while Data:
                for i in range(len(Data[1:])/3):
                    n=1
                    Ref= Data[i*3 + n]
                    n+=1
                    Position = Data[i*3 + n]
                    n+=1
                    Read = Data[i*3 + n]
                    n+=1
                    if Read not in Fusions:
                        if cfg.Ends > 0:
                            Read = Read[cfg.Ends:-cfg.Ends]
                            Position = str(int(Position) + cfg.Ends)
                        else:
                            pass
                    else:
                        pass
                    Event = Ref + "_@_" + Position + "_@_" + Read
                    if Data[0] in Name_Events:
                        Name_Events[Data[0]].append(Event)
                    else:
                        Name_Events[Data[0]] = [Event]
                    if Read in Fusions and cfg.Mode2 != 'Nucs':
                        if Event not in FusionLibCount:
                            FusionLibCount[Event][1] = 1
                        else:
                            FusionLibCount[Event][1] += 1
                    else:
                        pass
                Data = Input.readline().split()

##      -------------------------------------------------------------------------------------------
##      Take arguments from command line, and send them to the config file for cross-module access
##      -------------------------------------------------------------------------------------------

#       -----------------------------------------------------
#       -----------------------------------------------------

def AATest(Nuc):
    if cfg.ORFStartNuc > Nuc:
        return False
    elif Nuc > cfg.ORFFinishNuc:
        return False
    else:
        if (Nuc - cfg.ORFStartNuc+3)/3.0 == (Nuc - cfg.ORFStartNuc+3)/3:
            return (Nuc - cfg.ORFStartNuc+3)/3
        else:
            return False

def Translate(Seq):
    ProteinSeq = ''
    for i in range(len(Seq)/3):
        Codon = Seq[(i*3):(i*3)+3]
        if 'N' in Codon:
            AA = 'X'
            ##This may result in Translate() reading through a stop codon.
        else:
            AA = Codons[Codon]
            if AA == 'STOP':
                break
            else:
                pass
        ProteinSeq += AA
    return ProteinSeq

def TranslateSegment(Gene, Segment, FromNuc):
    x= False
    Trim = 0
    if FromNuc < cfg.ORFFinishNuc - 3:
        ## reduce cfg.ORFFinishNuc by three to prevent translation that begins at STOP codon
        while not x and FromNuc < cfg.ORFFinishNuc:
            x = AATest(FromNuc)
            if not x:
                Trim += 1
                FromNuc += 1
        NewSegment=Segment[Trim:]
        NewSegment=NewSegment[:((len(NewSegment)/3)*3)]
        ## Returns Gene, x (which is either AA Start Coord or false/0), and Peptide Seq
        return [Gene, x, Translate(NewSegment)]
    else:
        return [Gene, 0, '']

#       -----------------------------------------------------
#       -----------------------------------------------------

def FindNegativeFusionsFromRead(Gene, Pos, ReadLen):
    Negs = set()
    ReadTo = Pos + ReadLen - 1
    for i in FusionLibCount:
        ## i in format Ref + "_@_" + Positions + "_@_" + FusionType.
        ## Example is 'NC_004146.1_FHV_RNA1.seq_to_NC_004146.1_FHV_RNA1.seq_@_345_to_915_@_Rec'
        Event = i
        i = i.split("_@_")
        Ref = i[0]
        try:
            RefFrom, RefTo = Ref.split('_to_')
        except:
            RefFrom, RefTo = Ref, Ref
            #In this case is uDel or Ins
        RecFrom = int(i[1].split("_")[0])
        RecTo = int(i[1].split("_")[2])
        if RefFrom == RefTo and Gene == RefFrom:
            if RecFrom < RecTo:
                ##Deletion - Read must map inside deletion to negate
                if cfg.uInDel > (RecTo - RecFrom):
                    # Small deletion
                    if ReadTo - cfg.XInDel  >= RecTo and Pos + cfg.XInDel <= RecFrom:
                        ##Read maps over small deletion by more than X nucs on either side.
                        if FusionLibCount[Event][1] >= cfg.Min_Fusion_Coverage:
                            Negs.add(Event)
                        else:
                            #fewer detected events than will be reported
                            pass
                    else:
                        #No Data Overlap
                        pass
                else:
                    ## Large Deletion
                    Count = min(ReadTo, RecTo) - max(Pos, RecFrom)
                    if Count >= cfg.XRec:
                        ##Read maps inside deletion by more than cfg.X nucs.
                        if FusionLibCount[Event][1] >= cfg.Min_Fusion_Coverage:
                            Negs.add(Event)
                        else:
                            #fewer detected events than will be reported
                            pass
                    else:
                        #No Data overlap
                        pass
            else:
                ##Insertion - Read must continuously map over insertion sites to negate
                if len(i[1].split("_")[1]) <= cfg.uInDel:
                    #Small Insertion
                    if Pos <= RecTo - cfg.XInDel and ReadTo >= RecFrom + cfg.XInDel:
                        if FusionLibCount[Event][1] >= cfg.Min_Fusion_Coverage:
                            Negs.add(Event)
                        else:
                            #fewer detected events than will be reported
                            pass
                    else:
                        #No Data overlap
                        pass
                else:
                    ##Insertion - Read must continuously map over insertion sites to negate
                    if Pos <= RecTo - cfg.XRec and ReadTo >= RecFrom + cfg.XRec:
                        if FusionLibCount[Event][1] >= cfg.Min_Fusion_Coverage:
                            Negs.add(Event)
                        else:
                            #fewer detected events than will be reported
                            pass
                    else:
                        #No Data overlap
                        pass
        else:
            ## Inter-genic recombination event
            ## If multiple inter-genic recombination event are in a single read this may throw up erroneous results eg:
            ## Mapping 1, gene1-gene2, mapping 2, gene2-gene1, mapping 1.
            ## Suggest limiting number of recombination events per read, or trimming reads.
            if RefFrom == Gene and RefTo !=Gene:
                ## Read mapping is over first gene upstream of recombination event
                ## So, if mapping is downstream it negates recomb (ReadTo gives end of mapping).
                if ReadTo > RecFrom:
                    if FusionLibCount[Event][1] >= cfg.Min_Fusion_Coverage:
                        Negs.add(Event)
                    else:
                        #fewer detected events than will be reported
                        pass
                else:
                    ##Compatible
                    pass
            elif RefTo == Gene and RefFrom != Gene:
                ## Read mapping is over second gene downstream of recombination event
                ## So, if mapping is upstream it negates Recomb (Pos gives start of mapping).
                if Pos < RecTo:
                    if FusionLibCount[Event][1] >= cfg.Min_Fusion_Coverage:
                        Negs.add(Event)
                    else:
                        #fewer detected events than will be reported
                        pass
                else:
                    ##Compatible
                    pass

            else:
                pass
    return Negs

#       -----------------------------------------------------
#       -----------------------------------------------------

def AddToNucDict(Dict,i,jCoords, Count):
##    Dict = Gene1 + '_to_' + Gene2
    #Dict is same name as generated in Make_Gene_Dicts()
        if cfg.Mode =='NT':
                iNum = i[1]*4
        else:
                iNum = i[1]*20
        for j in jCoords:
            if i[0] < j[0]:
                cfg.Total_Dicts[Dict][i[0]][j[0]][iNum + j[1]] += Count
            else:
                break

def AddToMixDict(Dict,iCoords,jRec, Count, FusionCoord):
        ##Dict is Dict name, i is from iCoords, j is jRec.
        for i in iCoords:
            cfg.Total_Dicts[Dict][i[0]][jRec][i[1]*2 + FusionCoord] += Count

def AddToTotal_Matrices(i, j, Count):
    ##i and j are events in format: [Gene, Position, Read] for end-to-end mapping
    ## e.g. 'NC_004146.1_FHV_RNA1.seq_@_284_@_CCGTGACGCCTAAGCGTGTCGCAGACAA'
    iRec = i
    jRec = j
    i = i.split('_@_')
    j = j.split('_@_')
    if i[2] not in Fusions and j[2] not in Fusions:
                #Mapping vs mapping
        if int(i[1]) > int(j[1]):
            i,j = j,i
        else:
            pass
        iCoords = []
        jCoords = []
        Gene1 = i[0]
        Gene2 = j[0]
        if Gene1 == Gene2:
            ##Not handling non-homologous recombination events
            #Coords in format: [(nuc1, baseindex), (nuc2, baseindex), (nuc3, baseindex), ...]
            n = 0
            if cfg.Mode == 'AA':
                i = TranslateSegment(Gene1, i[2], int(i[1]))
                j = TranslateSegment(Gene2, j[2], int(j[1]))
                CoordLib = AAIndex
            else:
                CoordLib = Bases
            for k in i[-1]:
                if k in CoordLib:
                    Coord = int(i[-2])+n
                    if Coord in cfg.AllowedCoords[Gene1]:
                        iCoords.append((Coord, CoordLib.index(k)))
                    else:
                        ## Coord not sufficiently variable in pileup data
                        pass
                else:
                    pass
                n+=1
            n=0
            for k in j[-1]:
                if k in CoordLib:
                    Coord = int(j[-2])+n
                    if Coord in cfg.AllowedCoords[Gene2]:
                        jCoords.insert(0,(Coord, CoordLib.index(k)))
                    else:
                        ## Coord not sufficiently variable in pileup data
                        pass
                else:
                    pass
                n+=1
            #Dict is same name as generated in Make_Gene_Dicts()
            Dict = Gene1 + '_to_' + Gene2
            [AddToNucDict(Dict,i,jCoords, Count) for i in iCoords]
        else:
            pass
    elif i[2] in Fusions and j[2] in Fusions and cfg.Mode2 != 'Nucs':
        #Two Recombination Events
        if iRec != jRec:
            #A recombination event cannot associate with itself
            ##To prevent duplication of matrices in Recombination Dicts:
            if int(i[1].split("_")[2]) > int(j[1].split("_")[2]):
                iRec, jRec = jRec, iRec
            else:
                pass
            cfg.Total_Dicts['Rec_to_Rec'][iRec][jRec][3] += Count
        else:
            pass
    elif cfg.Mode2 != 'Nucs':
        ##Mix
        if i[2] in Fusions:
            #Order is mapping vs Rec
            i,j = j,i
            jRec, iRec = iRec, jRec
        else:
            pass
        iCoords = []
        Gene1 = i[0]
        n = 0
        if cfg.Mode == 'AA':
            i = TranslateSegment(Gene1, i[2], int(i[1]))
            CoordLib = AAIndex
        else:
            CoordLib = Bases
        for k in i[-1]:
            if k in CoordLib:
                Coord = int(i[-2])+n
                if Coord in cfg.AllowedCoords[Gene1]:
                    iCoords.append((Coord, CoordLib.index(k)))
                else:
                    pass
            else:
                pass
            n+=1
        Dict = Gene1 + '_to_Rec'
        for m in iCoords:
            cfg.Total_Dicts[Dict][m[0]][jRec][m[1]*2 + 1] += Count
    else:
                pass


##      -------------------------------------------------------------------------------------------
##      g(Event)
##      -------------------------------------------------------------------------------------------

def AddToPileup(Event, Count):
    ## Example Event :      ['NC_004146.1_FHV_RNA1.seq_@_284_@_CCGTGACGCCTAAGCGTGTCGCAGACAA',
    ##                       'NC_004146.1_FHV_RNA1.seq_to_NC_004146.1_FHV_RNA1.seq_@_310_to_959_@_Rec',
    ##                       'NC_004146.1_FHV_RNA1.seq_@_960_@_ACCTGGTTGGAACAGATTGGAGTATGTGACCGACAAGAATGAACTGCTGGTTTCCATCGGTCGAGAA']
    ## Do not add fusions, this is done elsewhere in FusionsLibs
    ## Care not to duplicate paired reads
    ## Pileup in format: PileupDict[Ref][Nuc] = [A,T,G,C]
    for i in Event:
        Record = []
        #prevents recording nt twice in paired reads
        Data = i.split("_@_")
        if Data[-1] not in Fusions:
            Ref = Data[0]
            n = int(Data[1])
            if cfg.Mode == 'NT':
                for j in Data[-1]:
                    if n not in Record and j != 'N':
                        cfg.PileUp[Ref][n][Bases.index(j)] += Count
                        Record.append(n)
                    else:
                        pass
                    n+=1
            else:
                AASeq = TranslateSegment(Data[0], Data[-1], n)
                n = AASeq[1]
                for j in AASeq[-1]:
                    if j == 'STOP':
                        break
                    else:
                        pass
                    if n not in Record and j != 'X' and n <= cfg.ProtLength:
                        cfg.PileUp[Ref][n][AAIndex.index(j)] += Count
                        Record.append(n)
                    else:
                        pass
                    n+=1
        else:
            pass

##      -------------------------------------------------------------------------------------------
##      EvaluatePileUp()
##      -------------------------------------------------------------------------------------------

def EvaluatePileUp():
    global FusionLibCount
    for Event in Events_List:
        Count = Event_Counts[Event]
        AddToPileup(Events_List[Event], Count)
    print('Writing Out PileUp list')
    pickle.dump(cfg.PileUp, open(cfg.Output_Tag + '.PileUp.py.pi','wb'))
    pickle.dump(FusionLibCount, open(cfg.Output_Tag + '.Fusion_Pileup.py.pi','wb'))
    ## Generate ErrorRates and allowed Coordinates
    for i in cfg.PileUp:
        for j in cfg.PileUp[i]:
            if j >= cfg.NtStart and j <= cfg.NtFinish:
                Data = cfg.PileUp[i][j]
                Coverage = sum(Data)
                if Coverage > cfg.Min_Coverage_Output and cfg.Mode2 != 'Recs':
                    ErrorRate = 1 - (max(Data)/float(Coverage))
                    if ErrorRate >= cfg.PileUp_Fraction:
                                    cfg.AllowedCoords[i].add(j)
                    else:
                            pass
                else:
                        pass
            else:
                    pass
    ##Filter FusionLib
    cfg.temp = defaultdict(dd2)
    for i in FusionLibCount:
        if FusionLibCount[i][1] >= cfg.Min_Fusion_Coverage:
            cfg.temp[i] = FusionLibCount[i]
        else:
            ##Should reduce low frequency InDel burden
            pass
    FusionLibCount = cfg.temp

    for i in cfg.AllowedCoords:
        print('Number of Allowed Coords in', i, 'is', len(cfg.AllowedCoords[i]))

##      -------------------------------------------------------------------------------------------
##      Populate()
##      -------------------------------------------------------------------------------------------
#@profile
def Populate():
    n = 0
    for Event in Events_List:
        ## Example Event :      ['NC_004146.1_FHV_RNA1.seq_@_284_@_CCGTGACGCCTAAGCGTGTCGCAGACAA',
        ##                       'NC_004146.1_FHV_RNA1.seq_to_NC_004146.1_FHV_RNA1.seq_@_310_to_959_@_Rec',
        ##                       'NC_004146.1_FHV_RNA1.seq_@_960_@_ACCTGGTTGGAACAGATTGGAGTATGTGACCGACAAGAATGAACTGCTGGTTTCCATCGGTCGAGAA']
        TotalUniqueEvents = len(Events_List)
        n5 = TotalUniqueEvents / 20
        n += 1
        if cfg.PrintReadCount:
            print(n)
        else:
            pass
        try:
            if n % n5 == 0:
                print(" .... %s Percent of Alignments Processed .... " % (n/n5*5))
            else:
                pass
        except:
            print(" Fewer than 20 reads being processed: Read %s " % (n))
        x = list(Events_List[Event])
        for i in x:
            if i.split('_@_')[2] in Fusions and i not in FusionLibCount:
                x.remove(i)
        Events_List[Event] = tuple(x)
        ##Populate Straight-forward Matrices
        Count = Event_Counts[Event]
        NumEvents = len(Events_List[Event])
        for i in range(NumEvents):
            for j in Events_List[Event][i:]:
                AddToTotal_Matrices(Events_List[Event][i], j, Count)
        ##If Recs Mode2 selected
        if cfg.Mode2 != 'Nucs':
                        ##Find Negative Recombination Events
                        kNegs = set()
                        for i in range(NumEvents):
                                k = Events_List[Event][i].split("_@_")
                                if k[2] not in Fusions:
                                        ##i.e. normal mapping eg.CCGTGACGCCTAAGCGTGTCGCAGACAA
                                        kNegs.update(FindNegativeFusionsFromRead(k[0], int(k[1]), len(k[2])))
                                else:
                                        pass

                        ##Populate Negative rec data
                        Dict = 'Rec_to_Rec'
                        for kNeg1 in kNegs:
                                for kNeg2 in kNegs:
                                        Pos1 = int(kNeg1.split("_@_")[1].split("_")[2])
                                        Pos2 = int(kNeg2.split("_@_")[1].split("_")[2])
                                        if Pos1 < Pos2 and kNeg1 != kNeg2:
                                                cfg.Total_Dicts[Dict][kNeg1][kNeg2][0] += Count
                                        else:
                                                pass
                        for i in range(NumEvents):
                                k = Events_List[Event][i].split("_@_")
                                ##Mix, Neg Rec to Normal Coord
                                if k[2] not in Fusions:
                                        kCoords = []
                                        m = 0
                                        if cfg.Mode == 'NT':
                                                for l in k[2]:
                                                        if l != 'N':
                                                                Coord = int(k[-2])+m
                                                                if Coord in cfg.AllowedCoords[k[0]]:
                                                                        kCoords.append((Coord, Bases.index(l)))
                                                                else:
                                                                        pass
                                                        else:
                                                                pass
                                                        m+=1
                                        else:
                                                kAAseq = TranslateSegment(k[0], k[2], int(k[1]))
                                                for l in kAAseq[-1]:
                                                        if l != 'X':
                                                                Coord = int(kAAseq[1])+m
                                                                if Coord in cfg.AllowedCoords[kAAseq[0]]:
                                                                        kCoords.append((Coord, AAIndex.index(l)))
                                                                else:
                                                                        pass
                                                        else:
                                                                pass
                                                        m+=1
                                        Dict = k[0] + '_to_Rec'
                                        [AddToMixDict(Dict,kCoords,kNeg, Count, 0) for kNeg in kNegs]
                                else:
                                        ##Recombination to Negative Recombination
                                        kRec = Events_List[Event][i]
                                        Dict = 'Rec_to_Rec'
                                        for kNeg in kNegs:
                                                if int(kRec.split("_@_")[1].split("_")[2]) > int(kNeg.split("_@_")[1].split("_")[2]):
                                                        cfg.Total_Dicts[Dict][kNeg][kRec][1] += Count
                                                else:
                                                        cfg.Total_Dicts[Dict][kRec][kNeg][2] += Count
                        for i in kNegs:
                                FusionLibCount[i][0] += Count
        else:
                        pass
##            except:
##                    print('Error in Event:/n', Event

##      -------------------------------------------------------------------------------------------
##      Run Main Script
##      -------------------------------------------------------------------------------------------

if __name__ == "__main__":

    ##Global variable made
    Bases = ['A','T','G','C']    # Nucleotide Index Numbers for order in Association Matrix given by Bases.index() command.
    Fusions = ['uDel','uIns','Ins','Rec']
    cfg.Genes_Lib = {}
    cfg.Total_Dicts = {}
    cfg.Filtered_Dicts = {}
    cfg.Coverage_Dicts = {}
    cfg.PileUp = {}
    cfg.AllowedCoords = {}
    Name_Events = {}
    Event_Counts = {}
    Events_List = {}
    Differences = []
    FusionLibCount = defaultdict(dd2)

    ##Extract Sequence Data
    print('Reading in FASTA File...')
    ReadFasta()

    ##Generate dicts of dicts for each type of matrix
    print('Generating unpopulated gene matrix/matrices')
    [Make_Gene_Dicts(i) for i in cfg.Genes_Lib]

    ##Generate dicts of dicts for fusion matrices
    print('Generating unpopulated fusion matrix/matrices')
    cfg.Total_Dicts['Rec_to_Rec'] = {}
    cfg.Filtered_Dicts['Rec_to_Rec'] = {}
    cfg.Coverage_Dicts['Rec_to_Rec'] = {}
    [Make_Fusion_Dicts(i) for i in cfg.Genes_Lib]

    ## Extract data from first SAMFile
    if cfg.SAMIN1:
        print('Extracting data from SAM file')
        Extract_Data_fromSAM(cfg.SAMIN1)
    else:
        pass

    ## Extract data from second SAMFile
    if cfg.SAMIN2:
        print('Extracting data from second/paired SAM file')
        Extract_Data_fromSAM(cfg.SAMIN2)
    else:
        pass

    ## Extract data from ViReMa Output
    if cfg.ViReMaIn1:
        print('Extracting data from ViReMa output')
        Extract_Data_fromViReMa_Output(cfg.ViReMaIn1)
    else:
        pass
    if cfg.ViReMaIn2:
        print('Extracting data from paired ViReMa output')
        Extract_Data_fromViReMa_Output(cfg.ViReMaIn2)
    else:
        pass

    ##Uniquify Dicts and count entries. Not essential, but otherwise mulitple identical events will repeatedly processed
    print('Compressing and counting non-unique alignments')
    for i in Name_Events:
        ##i is read name. Example:  'HISEQ:305:D2D4KACXX:4:1208:4266:49184'
        j = str(Name_Events[i])
        ## Example j :  ['NC_004146.1_FHV_RNA1.seq_@_284_@_CCGTGACGCCTAAGCGTGTCGCAGACAA',
        ##               'NC_004146.1_FHV_RNA1.seq_to_NC_004146.1_FHV_RNA1.seq_@_310_to_959_@_Rec',
        ##               'NC_004146.1_FHV_RNA1.seq_@_960_@_ACCTGGTTGGAACAGATTGGAGTATGTGACCGACAAGAATGAACTGCTGGTTTCCATCGGTCGAGAA']
        if j in Events_List:
            ##Events_Counts is a dictionary
            Event_Counts[j] += 1
        else:
            Events_List[j] = Name_Events[i]
            Event_Counts[j] = 1
    print("%s aligned reads were detected, %s of these alignments were unique" %(len(Name_Events), len(Events_List)))

    ## Generate Pileup List, ErrorRates and allowed Coordinates
    print('Generating PileUp list')
    EvaluatePileUp()
    
    cfg.Total_Dicts['Rec_to_Rec'] = {i:defaultdict(dd4) for i in FusionLibCount}
    cfg.Filtered_Dicts['Rec_to_Rec'] = {i:defaultdict(dd4) for i in FusionLibCount}
    cfg.Coverage_Dicts['Rec_to_Rec'] = {i:defaultdict(dd4) for i in FusionLibCount}

    if not cfg.PileUpOnly:
        ## Enter counts for each event into association matrix.  Slow on Python, >30x faster with pypy.
        print("Populating Association Matrix with alignment data...")
        Populate()

        ##Remove Under-populated Contigency tables
        Count = 0
        print("Finished Populating Association Matrix.")
        print("Filtering final Matrix.")
        for Lib in cfg.Total_Dicts:
            for i in cfg.Total_Dicts[Lib]:
                for j in cfg.Total_Dicts[Lib][i]:
                    if sum(cfg.Total_Dicts[Lib][i][j]) >= cfg.Min_Coverage_Output:
                        cfg.Filtered_Dicts[Lib][i][j] = cfg.Total_Dicts[Lib][i][j]
                        cfg.Coverage_Dicts[Lib][i][j] = 2
                    else:
                        Count += 1
                        cfg.Coverage_Dicts[Lib][i][j] = 1
        print(Count, "contingencies removed for having fewer that %s events mapped" %(cfg.Min_Coverage_Output))

        ##Save matrix for future (re)analyses.
        print('Pickle-dumping Populated Matrix/matrices')
        pickle.dump(cfg.Filtered_Dicts, open(cfg.Output_Tag + '.Total_Matrices.py.pi','wb'))
        pickle.dump(cfg.Coverage_Dicts, open(cfg.Output_Tag + '.Coverage_Matrices.py.pi','wb'))
    else:
        print("'PileUpOnly' option selected. Exiting without populating contingency tables.")


finish = time.time()
print("Time to complete in seconds: ", float(finish - start))

