#!/usr/bin/python

import csv
import os
import sys
from matplotlib import pyplot as plt
import numpy as np
from termcolor import colored


class StatisticPair:

    RawReadCount = 0
    TrimReadCount1 = 0
    TrimReadCount2 = 0
    ShortReadCount1 = 0
    ShortReadCount2 = 0
    TotGCSum1 = 0
    TotGCSum2 = 0
    TotLenSum1 = 0
    TotLenSum2 = 0
    TotGCSuma1 = 0
    TotGCSuma2 = 0
    TotLenSuma1 = 0
    TotLenSuma2 = 0
    CleanReadCount = 0
    NreadCount1 = 0         ##Reads containing more than Nbase % (given)
    NreadCount2 = 0
    AdtrimCount1 = 0
    AdtrimCount2 = 0
    FailQualCount1 = 0
    FailQualCount2 = 0
    TotACount1 = 0; TotACount2 = 0
    TotTCount1 = 0; TotTCount2 = 0
    TotGCount1 = 0; TotGCount2 = 0
    TotCCount1 = 0; TotCCount2 = 0
    TotNCount1 = 0; TotNCount2 = 0
    TotAaCount1 = 0; TotAaCount2 = 0
    TotTaCount1 = 0; TotTaCount2 = 0
    TotGaCount1 = 0; TotGaCount2 = 0
    TotCaCount1 = 0; TotCaCount2 = 0
    TotNaCount1 = 0; TotNaCount2 = 0
    TotQual1_1 = 0; TotQual2_1 = 0
    TotQual1_2 = 0; TotQual2_2 = 0
    TotQual1_3 = 0; TotQual2_3 = 0
    TotQual1_4 = 0; TotQual2_4 = 0
    TotTrimLen1 = 0; TotTrimLen2 = 0
    hash_gcList1_1 = 0; hash_gcList2_1 = 0
    hash_gcList1_2 = 0; hash_gcList2_2 = 0
    hash_gcList1_3 = 0; hash_gcList2_3 = 0
    hash_gcList1_4 = 0; hash_gcList2_4 = 0
    hash_gcList1_5 = 0; hash_gcList2_5 = 0
    hash_gcList1_6 = 0; hash_gcList2_6 = 0
    hash_gcList1_7 = 0; hash_gcList2_7 = 0
    hash_gcList1_8 = 0; hash_gcList2_8 = 0
    hash_gcList1_9 = 0; hash_gcList2_9 = 0
    hash_gcList1_10 = 0; hash_gcList2_10 = 0
    hash_gcList1a_1 = 0; hash_gcList2a_1 = 0
    hash_gcList1a_2 = 0; hash_gcList2a_2 = 0
    hash_gcList1a_3 = 0; hash_gcList2a_3 = 0
    hash_gcList1a_4 = 0; hash_gcList2a_4 = 0
    hash_gcList1a_5 = 0; hash_gcList2a_5 = 0
    hash_gcList1a_6 = 0; hash_gcList2a_6 = 0
    hash_gcList1a_7 = 0; hash_gcList2a_7 = 0
    hash_gcList1a_8 = 0; hash_gcList2a_8 = 0
    hash_gcList1a_9 = 0; hash_gcList2a_9 = 0
    hash_gcList1a_10 = 0; hash_gcList2a_10 = 0
    hash_qual1_1 = 0; hash_qual2_1 = 0
    hash_qual1_2 = 0; hash_qual2_2 = 0
    hash_qual1_3 = 0; hash_qual2_3 = 0
    hash_qual1_4 = 0; hash_qual2_4 = 0
    hash_qual1_5 = 0; hash_qual2_5 = 0
    hash_qual1_6 = 0; hash_qual2_6 = 0
    hash_qual1_7 = 0; hash_qual2_7 = 0
    hash_qual1_8 = 0; hash_qual2_8 = 0
    hash_qual1_9 = 0; hash_qual2_9 = 0
    hash_qual1_10 = 0; hash_qual2_10 = 0
    hash_qual1_11 = 0; hash_qual2_11 = 0
    hash_qual1_12 = 0; hash_qual2_12 = 0
    hash_qual1_13 = 0; hash_qual2_13 = 0
    hash_qual1_14 = 0; hash_qual2_14 = 0
    hash_qual1_15 = 0; hash_qual2_15 = 0
    hash_qual1_16 = 0; hash_qual2_16 = 0
    hash_qual1_17 = 0; hash_qual2_17 = 0
    hash_qual1_18 = 0; hash_qual2_18 = 0
    hash_qual1_19 = 0; hash_qual2_19 = 0
    hash_qual1_20 = 0; hash_qual2_20 = 0
    hash_qual1_21 = 0; hash_qual2_21 = 0
    hash_qual1a_1 = 0; hash_qual2a_1 = 0
    hash_qual1a_2 = 0; hash_qual2a_2 = 0
    hash_qual1a_3 = 0; hash_qual2a_3 = 0
    hash_qual1a_4 = 0; hash_qual2a_4 = 0
    hash_qual1a_5 = 0; hash_qual2a_5 = 0
    hash_qual1a_6 = 0; hash_qual2a_6 = 0
    hash_qual1a_7 = 0; hash_qual2a_7 = 0
    hash_qual1a_8 = 0; hash_qual2a_8 = 0
    hash_qual1a_9 = 0; hash_qual2a_9 = 0
    hash_qual1a_10 = 0; hash_qual2a_10 = 0
    hash_qual1a_11 = 0; hash_qual2a_11 = 0
    hash_qual1a_12 = 0; hash_qual2a_12 = 0
    hash_qual1a_13 = 0; hash_qual2a_13 = 0
    hash_qual1a_14 = 0; hash_qual2a_14 = 0
    hash_qual1a_15 = 0; hash_qual2a_15 = 0
    hash_qual1a_16 = 0; hash_qual2a_16 = 0
    hash_qual1a_17 = 0; hash_qual2a_17 = 0
    hash_qual1a_18 = 0; hash_qual2a_18 = 0
    hash_qual1a_19 = 0; hash_qual2a_19 = 0
    hash_qual1a_20 = 0; hash_qual2a_20 = 0
    hash_qual1a_21 = 0; hash_qual2a_21 = 0
    TotQual1Sum = 0; TotQual2Sum = 0
    TotQual1aSum = 0; TotQual2aSum = 0
    NcontReads1 = 0; NcontReads2 = 0
    NcontReads1a = 0; NcontReads2a = 0
    TrimQualCt1 = 0; TrimQualCt2 = 0
    LenList1 = []; LenList2 = []
    LenList1a = []; LenList2a = []

    def __init__(self, InFile):
        InFile = open(InFile, 'r')
        InFile = csv.reader(InFile, delimiter='\t')
        for rec in InFile:
            self.RawReadCount += int(rec[0])
            self.TrimReadCount1 += int(rec[1])
            self.TrimReadCount2 += int(rec[2])
            self.ShortReadCount1 += int(rec[3])
            self.ShortReadCount2 += int(rec[4])
            self.TotGCSum1 += float(rec[6])
            self.TotGCSum2 += float(rec[7])
            self.TotLenSum1 += float(rec[8])
            self.TotLenSum2 += float(rec[9])
            self.TotGCSuma1 += float(rec[11])
            self.TotGCSuma2 += float(rec[12])
            self.TotLenSuma1 += float(rec[13])
            self.TotLenSuma2 += float(rec[14])
            self.CleanReadCount += int(rec[15])
            self.NreadCount1 += int(rec[16])
            self.NreadCount2 += int(rec[17])
            self.AdtrimCount1 += int(rec[18])
            self.AdtrimCount2 += int(rec[19])
            self.FailQualCount1 += int(rec[21])
            self.FailQualCount2 += int(rec[22])
            self.TotACount1 += int(rec[23])
            self.TotACount2 += int(rec[24])
            self.TotTCount1 += int(rec[25])
            self.TotTCount2 += int(rec[26])
            self.TotGCount1 += int(rec[27])
            self.TotGCount2 += int(rec[28])
            self.TotCCount1 += int(rec[29])
            self.TotCCount2 += int(rec[30])
            self.TotNCount1 += int(rec[31])
            self.TotNCount2 += int(rec[32])
            self.TotAaCount1 += int(rec[33])
            self.TotAaCount2 += int(rec[34])
            self.TotTaCount1 += int(rec[35])
            self.TotTaCount2 += int(rec[36])
            self.TotGaCount1 += int(rec[37])
            self.TotGaCount2 += int(rec[38])
            self.TotCaCount1 += int(rec[39])
            self.TotCaCount2 += int(rec[40])
            self.TotNaCount1 += int(rec[41])
            self.TotNaCount2 += int(rec[42])
            self.TotQual1_1 += int(rec[43])
            self.TotQual1_2 += int(rec[44])
            self.TotQual1_3 += int(rec[45])
            self.TotQual1_4 += int(rec[46])
            self.TotQual2_1 += int(rec[47])
            self.TotQual2_2 += int(rec[48])
            self.TotQual2_3 += int(rec[49])
            self.TotQual2_4 += int(rec[50])
            self.TotTrimLen1 += int(rec[59])
            self.TotTrimLen2 += int(rec[60])
            self.hash_gcList1_1 += int(rec[61]); self.hash_gcList2_1 += int(rec[71])
            self.hash_gcList1_2 += int(rec[62]); self.hash_gcList2_2 += int(rec[72])
            self.hash_gcList1_3 += int(rec[63]); self.hash_gcList2_3 += int(rec[73])
            self.hash_gcList1_4 += int(rec[64]); self.hash_gcList2_4 += int(rec[74])
            self.hash_gcList1_5 += int(rec[65]); self.hash_gcList2_5 += int(rec[75])
            self.hash_gcList1_6 += int(rec[66]); self.hash_gcList2_6 += int(rec[76])
            self.hash_gcList1_7 += int(rec[67]); self.hash_gcList2_7 += int(rec[77])
            self.hash_gcList1_8 += int(rec[68]); self.hash_gcList2_8 += int(rec[78])
            self.hash_gcList1_9 += int(rec[69]); self.hash_gcList2_9 += int(rec[79])
            self.hash_gcList1_10 += int(rec[70]); self.hash_gcList2_10 += int(rec[80])
            self.hash_gcList1a_1 += int(rec[81]); self.hash_gcList2a_1 += int(rec[91])
            self.hash_gcList1a_2 += int(rec[82]); self.hash_gcList2a_2 += int(rec[92])
            self.hash_gcList1a_3 += int(rec[83]); self.hash_gcList2a_3 += int(rec[93])
            self.hash_gcList1a_4 += int(rec[84]); self.hash_gcList2a_4 += int(rec[94])
            self.hash_gcList1a_5 += int(rec[85]); self.hash_gcList2a_5 += int(rec[95])
            self.hash_gcList1a_6 += int(rec[86]); self.hash_gcList2a_6 += int(rec[96])
            self.hash_gcList1a_7 += int(rec[87]); self.hash_gcList2a_7 += int(rec[97])
            self.hash_gcList1a_8 += int(rec[88]); self.hash_gcList2a_8 += int(rec[98])
            self.hash_gcList1a_9 += int(rec[89]); self.hash_gcList2a_9 += int(rec[99])
            self.hash_gcList1a_10 += int(rec[90]); self.hash_gcList2a_10 += int(rec[100])
            self.TotQual1Sum += float(rec[101]); self.TotQual2Sum += float(rec[103])
            self.TotQual1aSum += float(rec[102]); self.TotQual2aSum += float(rec[104])
            self.NcontReads1 += int(rec[105]); self.NcontReads2 += int(rec[107])
            self.NcontReads1a += int(rec[106]); self.NcontReads2a += int(rec[108])
            self.TrimQualCt1 += int(rec[109]); self.TrimQualCt2 += int(rec[110])
            self.hash_qual1_1 += int(rec[111]); self.hash_qual2_1 += int(rec[153])
            self.hash_qual1_2 += int(rec[112]); self.hash_qual2_2 += int(rec[154])
            self.hash_qual1_3 += int(rec[113]); self.hash_qual2_2 += int(rec[155])
            self.hash_qual1_4 += int(rec[114]); self.hash_qual2_4 += int(rec[156])
            self.hash_qual1_5 += int(rec[115]); self.hash_qual2_5 += int(rec[157])
            self.hash_qual1_6 += int(rec[116]); self.hash_qual2_6 += int(rec[158])
            self.hash_qual1_7 += int(rec[117]); self.hash_qual2_7 += int(rec[159])
            self.hash_qual1_8 += int(rec[118]); self.hash_qual2_8 += int(rec[160])
            self.hash_qual1_9 += int(rec[119]); self.hash_qual2_9 += int(rec[161])
            self.hash_qual1_10 += int(rec[120]); self.hash_qual2_10 += int(rec[162])
            self.hash_qual1_11 += int(rec[121]); self.hash_qual2_11 += int(rec[163])
            self.hash_qual1_12 += int(rec[122]); self.hash_qual2_12 += int(rec[164])
            self.hash_qual1_13 += int(rec[123]); self.hash_qual2_13 += int(rec[165])
            self.hash_qual1_14 += int(rec[124]); self.hash_qual2_14 += int(rec[166])
            self.hash_qual1_15 += int(rec[125]); self.hash_qual2_15 += int(rec[167])
            self.hash_qual1_16 += int(rec[126]); self.hash_qual2_16 += int(rec[168])
            self.hash_qual1_17 += int(rec[127]); self.hash_qual2_17 += int(rec[169])
            self.hash_qual1_18 += int(rec[128]); self.hash_qual2_18 += int(rec[170])
            self.hash_qual1_19 += int(rec[129]); self.hash_qual2_19 += int(rec[171])
            self.hash_qual1_20 += int(rec[130]); self.hash_qual2_20 += int(rec[172])
            self.hash_qual1_21 += int(rec[131]); self.hash_qual2_21 += int(rec[173])
            self.hash_qual1a_1 += int(rec[132]); self.hash_qual2a_1 += int(rec[174])
            self.hash_qual1a_2 += int(rec[133]); self.hash_qual2a_2 += int(rec[175])
            self.hash_qual1a_3 += int(rec[134]); self.hash_qual2a_3 += int(rec[176])
            self.hash_qual1a_4 += int(rec[135]); self.hash_qual2a_4 += int(rec[177])
            self.hash_qual1a_5 += int(rec[136]); self.hash_qual2a_5 += int(rec[178])
            self.hash_qual1a_6 += int(rec[137]); self.hash_qual2a_6 += int(rec[179])
            self.hash_qual1a_7 += int(rec[138]); self.hash_qual2a_7 += int(rec[180])
            self.hash_qual1a_8 += int(rec[139]); self.hash_qual2a_8 += int(rec[181])
            self.hash_qual1a_9 += int(rec[140]); self.hash_qual2a_9 += int(rec[182])
            self.hash_qual1a_10 += int(rec[141]); self.hash_qual2a_10 += int(rec[183])
            self.hash_qual1a_11 += int(rec[142]); self.hash_qual2a_11 += int(rec[184])
            self.hash_qual1a_12 += int(rec[143]); self.hash_qual2a_12 += int(rec[185])
            self.hash_qual1a_13 += int(rec[144]); self.hash_qual2a_13 += int(rec[186])
            self.hash_qual1a_14 += int(rec[145]); self.hash_qual2a_14 += int(rec[187])
            self.hash_qual1a_15 += int(rec[146]); self.hash_qual2a_15 += int(rec[188])
            self.hash_qual1a_16 += int(rec[147]); self.hash_qual2a_16 += int(rec[189])
            self.hash_qual1a_17 += int(rec[148]); self.hash_qual2a_17 += int(rec[190])
            self.hash_qual1a_18 += int(rec[149]); self.hash_qual2a_18 += int(rec[191])
            self.hash_qual1a_19 += int(rec[150]); self.hash_qual2a_19 += int(rec[192])
            self.hash_qual1a_20 += int(rec[151]); self.hash_qual2a_20 += int(rec[193])
            self.hash_qual1a_21 += int(rec[152]); self.hash_qual2a_21 += int(rec[194])
            self.LenList1.append(float(rec[51]))
            self.LenList1.append(float(rec[52]))
            self.LenList2.append(float(rec[53]))
            self.LenList2.append(float(rec[54]))
            self.LenList1a.append(float(rec[55]))
            self.LenList1a.append(float(rec[56]))
            self.LenList2a.append(float(rec[57]))
            self.LenList2a.append(float(rec[58]))
#       if no output obtained
        if int(self.CleanReadCount) == 0:
            print(colored("\nError: No Sequence reads reported as filtered\n", "red"))
            sys.exit()

    def stat_pair(self, Nbase, QualThresh, MinSize, Adapter, Trim, MinLen, File1, File2, QualFormat, OutputDir):
        os.chdir(OutputDir)
        TotBases1 = self.TotACount1+self.TotTCount1+self.TotGCount1+self.TotCCount1
        TotBasesWithN1 = self.TotACount1+self.TotTCount1+self.TotGCount1+self.TotCCount1+self.TotNCount1
        TotBases2 = self.TotACount2+self.TotTCount2+self.TotGCount2+self.TotCCount2
        TotBasesWithN2 = self.TotACount2+self.TotTCount2+self.TotGCount2+self.TotCCount2+self.TotNCount2
        TotBases1a = self.TotAaCount1+self.TotTaCount1+self.TotGaCount1+self.TotCaCount1
        TotBasesWithN1a = self.TotAaCount1+self.TotTaCount1+self.TotGaCount1+self.TotCaCount1+self.TotNaCount1
        TotBases2a = self.TotAaCount2+self.TotTaCount2+self.TotGaCount2+self.TotCaCount2
        TotBasesWithN2a = self.TotAaCount2+self.TotTaCount2+self.TotGaCount2+self.TotCaCount2+self.TotNaCount2
        MinLenRead1 = min(self.LenList1)
        MinLenRead2 = min(self.LenList2)
        MaxLenRead1 = max(self.LenList1)
        MaxLenRead2 = max(self.LenList2)
        try:
            self.LenList1a.remove(0)
            self.LenList1a.remove(1000000)
            self.LenList2a.remove(0)
            self.LenList2a.remove(1000000)
        except:
            pass
        MaxLenRead1a = max(self.LenList1a)
        MinLenRead1a = min(self.LenList1a)
        MinLenRead2a = min(self.LenList2a)
        MaxLenRead2a = max(self.LenList2a)
        Q301 = self.hash_qual1_16+self.hash_qual1_17+self.hash_qual1_18+self.hash_qual1_19+self.hash_qual1_20+self.hash_qual1_21
        Q302 = self.hash_qual2_16+self.hash_qual2_17+self.hash_qual2_18+self.hash_qual2_19+self.hash_qual2_20+self.hash_qual2_21
        Q301a = self.hash_qual1a_16+self.hash_qual1a_17+self.hash_qual1a_18+self.hash_qual1a_19+self.hash_qual1a_20+self.hash_qual1a_21
        Q302a = self.hash_qual2a_16+self.hash_qual2a_17+self.hash_qual2a_18+self.hash_qual2a_19+self.hash_qual2a_20+self.hash_qual2a_21
        StatFile = open("Statistics.txt", 'w')
        StatFile.write("%-s\n" % "Parameters specified for filtering")
        StatFile.write("===============================================\n\n")
        if Nbase < 101:
            StatFile.write("%-50s\t%s\n" % ("% of Uncalled bases(N)", Nbase))
        StatFile.write("%-50s\t%d\n" % ("Mean quality value threshold", QualThresh))
        StatFile.write("%-50s\t%d\n" % ("Minimum size of reads", int(MinSize)))
        if Adapter:
            StatFile.write("%-50s\t%-s\n" % ("Adapter sequences", Adapter))
        if Trim:
            StatFile.write("%-50s\t%-s\n" % ("Filtering mode for quality value", "Trim"))
            StatFile.write("%-50s\t%d\n" % ("Read with minimum length to keep after trimming", MinLen))
        else:
            StatFile.write("%-50s\t%-s\n" % ("Filtering mode for quality value", "Filter"))
        StatFile.write("\n\n")
        StatFile.write("%-s\n" % "Filtering Statistics for given Paired end files")
        StatFile.write("==================================================================================================================================\n\n")
        StatFile.write("%-65s\t%-20s\t%-20s\n" % ("Input file", File1, File2))
        if QualFormat == 1:
            StatFile.write("%-65s\t%-20s\t%-20s\n" % ("Quality format", "illumina 1.8+", "illumina 1.8+"))
        elif QualFormat == 2:
            StatFile.write("%-65s\t%-20s\t%-20s\n" % ("Quality format", "illumina 1.6+", "illumina 1.6+"))
        elif QualFormat == 3:
            StatFile.write("%-65s\t%-20s\t%-20s\n" % ("Quality format", "Sanger", "Sanger"))
        StatFile.write("%-65s\t%-20s\t%-20s\n" % ("Total number of reads analyzed", self.RawReadCount, self.RawReadCount))
        StatFile.write("%-65s\t%-20s\t%-20s\n" % ("Total Bases (ATGC) in unfiltered reads (bp)", TotBases1, TotBases2))
        StatFile.write("%-65s\t%-20s\t%-20s\n" % ("Total Bases (ATGC) in filtered reads (bp)", TotBases1a, TotBases2a))
        StatFile.write("%-65s\t%-20s\t%-20s\n" % ("Minimum size (bp) of unfiltered reads", MinLenRead1, MinLenRead2))
        StatFile.write("%-65s\t%-20s\t%-20s\n" % ("Maximum size (bp) of unfiltered reads", MaxLenRead1, MaxLenRead2))
        StatFile.write("%-65s\t%-20s\t%-20s\n" % ("Mean size (bp) of unfiltered reads", float(self.TotLenSum1/self.RawReadCount), float(self.TotLenSum2/self.RawReadCount)))
        StatFile.write("%-65s\t%-20s\t%-20s\n" % ("Average Quality value for unfiltered reads", float(self.TotQual1Sum/TotBasesWithN1), float(self.TotQual2Sum/TotBasesWithN2)))

        StatFile.write("%-65s\t%-20s\t%-20s\n" % ("Average Quality value for filtered reads", float(self.TotQual1aSum/TotBasesWithN1a), float(self.TotQual2aSum/TotBasesWithN2a)))
        if MinSize:
            StatFile.write("%-65s\t%-20s\t%-20s\n" % ("Total number of reads below minimum size", self.ShortReadCount1, self.ShortReadCount2))
        StatFile.write("%-65s\t%-20s\t%-20s\n" % ("Total unfiltered reads containing at least one uncalled base(N)", self.NcontReads1, self.NcontReads2))
        StatFile.write("%-65s\t%-20s\t%-20s\n" % ("Total filtered reads containing at least one uncalled base(N)", self.NcontReads1a, self.NcontReads2a))
        if Nbase:
            StatFile.write("%-65s\t%-20s\t%-20s\n" % ("Reads filtered out with more than given % of N", self.NreadCount1, self.NreadCount2))
        if Trim is None:
            StatFile.write("%-65s\t%-20s\t%-20s\n" % ("Reads filtered out for quality", self.FailQualCount1, self.FailQualCount2))
        if Adapter:
            StatFile.write("%-65s\t%-20s\t%-20s\n" % ("Reads trimmed for adapter", self.AdtrimCount1, self.AdtrimCount2))
        StatFile.write("%-65s\t%-20s\t%-20s\n" % ("Average %GC content in unfiltered reads", float((self.TotGCSum1*100)/self.TotLenSum1), float((self.TotGCSum2*100)/self.TotLenSum2)))
        StatFile.write("%-65s\t%-20s\t%-20s\n" % ("Average %GC content in filtered reads", float((self.TotGCSuma1*100)/self.TotLenSuma1), float((self.TotGCSuma2*100)/self.TotLenSuma2)))
        if Trim:
            StatFile.write("%-65s\t%-20s\t%-20s\n" % ("Reads trimmed for quality", self.TrimQualCt1, self.TrimQualCt2))
        StatFile.write("%-65s\t%-20s\t%-20s\n" % ("Total unfiltered reads with >Q30 mean quality value", Q301, Q302))
        StatFile.write("%-65s\t%-20s\t%-20s\n" % ("Total filtered reads with >Q30 mean quality value ", Q301a, Q302a))
        StatFile.write("%-65s\t%-20s\t%-20s\n" % ("Total number of reads kept", self.CleanReadCount, self.CleanReadCount))
        if Trim is None:
            StatFile.write("%-65s\t%-20s\t%-20s\n" % ("Minimum size (bp) of filtered reads", MinLenRead1a, MinLenRead2a))
            StatFile.write("%-65s\t%-20s\t%-20s\n" % ("Maximum size (bp) of filtered reads", MaxLenRead1a, MaxLenRead2a))
            StatFile.write("%-65s\t%-20s\t%-20s\n" % ("Mean size (bp) of filtered reads", float(self.TotLenSuma1/self.CleanReadCount), float(self.TotLenSuma2/self.CleanReadCount)))
        else:
            StatFile.write("%-65s\t%-20s\t%-20s\n" % ("Minimum size (bp) of filtered reads", MinLenRead1a, MinLenRead2a))
            StatFile.write("%-65s\t%-20s\t%-20s\n" % ("Maximum size (bp) of filtered reads", MaxLenRead1a, MaxLenRead2a))
            StatFile.write("%-65s\t%-20s\t%-20s\n" % ("Mean size (bp) of filtered reads", float(self.TotLenSuma1/self.CleanReadCount), float(self.TotLenSuma2/self.CleanReadCount)))

        Remain = self.RawReadCount - self.CleanReadCount
        StatFile.write("\n\n\nThe total number of sequences removed:"+str(Remain)+"\n")
        StatFile.close()
#        self.quality_vis()

    def stat_vis(self, FileP1_basname, FileP2_basname):

        Gcbin = range(10, 101, 10)
        Gcvalue1 = []
        Gcvalue1a = []
        Gcvalue2 = []
        Gcvalue2a = []

        Gcvalue1.extend([self.hash_gcList1_1, self.hash_gcList1_2, self.hash_gcList1_3, self.hash_gcList1_4, self.hash_gcList1_5, self.hash_gcList1_6])
        Gcvalue1.extend([self.hash_gcList1_7, self.hash_gcList1_8, self.hash_gcList1_9, self.hash_gcList1_10])

        Gcvalue1a.extend([self.hash_gcList1a_1, self.hash_gcList1a_2, self.hash_gcList1a_3, self.hash_gcList1a_4, self.hash_gcList1a_5, self.hash_gcList1a_6])
        Gcvalue1a.extend([self.hash_gcList1a_7, self.hash_gcList1a_8, self.hash_gcList1a_9, self.hash_gcList1a_10])

        Gcvalue2.extend([self.hash_gcList2_1, self.hash_gcList2_2, self.hash_gcList2_3, self.hash_gcList2_4, self.hash_gcList2_5, self.hash_gcList2_6])
        Gcvalue2.extend([self.hash_gcList2_7, self.hash_gcList2_8, self.hash_gcList2_9, self.hash_gcList2_10])

        Gcvalue2a.extend([self.hash_gcList2a_1, self.hash_gcList2a_2, self.hash_gcList2a_3, self.hash_gcList2a_4, self.hash_gcList2a_5, self.hash_gcList2a_6])
        Gcvalue2a.extend([self.hash_gcList2a_7, self.hash_gcList2a_8, self.hash_gcList2a_9, self.hash_gcList2a_10])

        plt.figure(1)
        plt.plot(Gcbin, Gcvalue1a, 'g', Gcbin, Gcvalue1, 'r')
        plt.xticks([10, 20, 30, 40, 50, 60, 70, 80, 90, 100], rotation='horizontal')
        plt.legend(["Filtered", "Unfiltered"], loc=2)
        plt.ylabel("# Reads", fontweight='bold')
        plt.xlabel("%GC Content", fontweight='bold')
        plt.savefig(FileP1_basname+'_GCdist.png')
        plt.clf()

        plt.figure(2)
        plt.plot(Gcbin, Gcvalue2a, 'g', Gcbin, Gcvalue2, 'r')
        plt.xticks([10, 20, 30, 40, 50, 60, 70, 80, 90, 100], rotation='horizontal')
        plt.legend(["Filtered", "Unfiltered"], loc=2)
        plt.ylabel("# Reads", fontweight='bold')
        plt.savefig(FileP2_basname+'_GCdist.png')
        plt.clf()

        #for base composition
        Bcbin = np.arange(5)
        Bcvalue1 = []
        Bcvalue1a = []
        Bcvalue2 = []
        Bcvalue2a = []

        Bcvalue1.extend([self.TotACount1, self.TotTCount1, self.TotGCount1, self.TotCCount1, self.TotNCount1])
        Bcvalue1a.extend([self.TotAaCount1, self.TotTaCount1, self.TotGaCount1, self.TotCaCount1, self.TotNaCount1])
        Bcvalue2.extend([self.TotACount2, self.TotTCount2, self.TotGCount2, self.TotCCount2, self.TotNCount2])
        Bcvalue2a.extend([self.TotAaCount2, self.TotTaCount2, self.TotGaCount2, self.TotCaCount2, self.TotNaCount2])

        wid = 0.25
        fig = plt.figure()
#        fig, ax = plt.subplots()
        plt.bar(Bcbin, Bcvalue1a, wid, color='#6A5ACD', label='Filter', align='center')
        plt.bar(Bcbin + wid, Bcvalue1, wid, color='#CD5C5C', label='Unfilter', align='center')
        plt.ylabel('# Bases',fontweight='bold')
        plt.xticks([0, 1, 2, 3, 4], ['A', 'T', 'G', 'C', 'N'], fontsize=12, ha='center')
        plt.legend(bbox_to_anchor=(0., 1., 1., .102), loc=3, prop={'size':15}, ncol=2, mode="expand")
        plt.margins(0.05, None)  #space in x axis
        plt.savefig(FileP1_basname+'_Basedist.png')
        plt.clf()

#        fig, ax = plt.subplots()
        fig = plt.figure()
        plt.bar(Bcbin, Bcvalue2a, wid, color='#6A5ACD', label='Filter', align='center')
        plt.bar(Bcbin + wid, Bcvalue2, wid, color='#CD5C5C', label='Unfilter', align='center')
        plt.ylabel('# Bases',fontweight='bold')
        plt.xticks([0, 1, 2, 3, 4], ['A', 'T', 'G', 'C', 'N'], rotation='0', fontsize=12, ha='center')
        plt.legend(bbox_to_anchor=(0., 1., 1., .102), loc=3, prop={'size':15}, ncol=2, mode="expand")
        plt.margins(0.05, None)  #space in x axis
        plt.savefig(FileP2_basname+'_Basedist.png')
        plt.clf()

        #for quality
        Qualbin = range(2, 43, 2)
        Qualvalue1 = []
        Qualvalue1a = []
        Qualvalue2 = []
        Qualvalue2a = []     #read count under each bin

        Qualvalue1.extend([self.hash_qual1_1, self.hash_qual1_2, self.hash_qual1_3, self.hash_qual1_4, self.hash_qual1_5, self.hash_qual1_6])
        Qualvalue1.extend([self.hash_qual1_7, self.hash_qual1_8, self.hash_qual1_9, self.hash_qual1_10, self.hash_qual1_11, self.hash_qual1_12])
        Qualvalue1.extend([self.hash_qual1_13, self.hash_qual1_14, self.hash_qual1_15, self.hash_qual1_16, self.hash_qual1_17, self.hash_qual1_18])
        Qualvalue1.extend([self.hash_qual1_19, self.hash_qual1_20, self.hash_qual1_21])

        Qualvalue1a.extend([self.hash_qual1a_1, self.hash_qual1a_2, self.hash_qual1a_3, self.hash_qual1a_4, self.hash_qual1a_5, self.hash_qual1a_6])
        Qualvalue1a.extend([self.hash_qual1a_7, self.hash_qual1a_8, self.hash_qual1a_9, self.hash_qual1a_10, self.hash_qual1a_11, self.hash_qual1a_12])
        Qualvalue1a.extend([self.hash_qual1a_13, self.hash_qual1a_14, self.hash_qual1a_15, self.hash_qual1a_16, self.hash_qual1a_17, self.hash_qual1a_18])
        Qualvalue1a.extend([self.hash_qual1a_19, self.hash_qual1a_20, self.hash_qual1a_21])

        Qualvalue2.extend([self.hash_qual2_1, self.hash_qual2_2, self.hash_qual2_3, self.hash_qual2_4, self.hash_qual2_5, self.hash_qual2_6])
        Qualvalue2.extend([self.hash_qual2_7, self.hash_qual2_8, self.hash_qual2_9, self.hash_qual2_10, self.hash_qual2_11, self.hash_qual2_12])
        Qualvalue2.extend([self.hash_qual2_13, self.hash_qual2_14, self.hash_qual2_15, self.hash_qual2_16, self.hash_qual2_17, self.hash_qual2_18])
        Qualvalue2.extend([self.hash_qual2_19, self.hash_qual2_20, self.hash_qual2_21])

        Qualvalue2a.extend([self.hash_qual2a_1, self.hash_qual2a_2, self.hash_qual2a_3, self.hash_qual2a_4, self.hash_qual2a_5, self.hash_qual2a_6])
        Qualvalue2a.extend([self.hash_qual2a_7, self.hash_qual2a_8, self.hash_qual2a_9, self.hash_qual2a_10, self.hash_qual2a_11, self.hash_qual2a_12])
        Qualvalue2a.extend([self.hash_qual2a_13, self.hash_qual2a_14, self.hash_qual2a_15, self.hash_qual2a_16, self.hash_qual2a_17, self.hash_qual2a_18])
        Qualvalue2a.extend([self.hash_qual2a_19, self.hash_qual2a_20, self.hash_qual2a_21])

        #for plotting quality comparison between between filtered and unfiltered data
        plt.figure(5)
        plt.plot(Qualbin, Qualvalue1a, 'g', Qualbin, Qualvalue1, 'r')
        plt.xticks([0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42], rotation='vertical')
        plt.legend(["Filtered", "Unfiltered"], loc=2)
        plt.ylabel("# Reads", fontweight='bold')
        plt.savefig(FileP1_basname+'_Qualdist.png')
        plt.clf()

        plt.figure(6)
        plt.plot(Qualbin, Qualvalue2a, 'g', Qualbin, Qualvalue2, 'r')
        plt.xticks([0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42], rotation='vertical')
        plt.legend(["Filtered", "Unfiltered"], loc=2)
        plt.ylabel("# Reads", fontweight='bold')
        plt.savefig(FileP2_basname+'_Qualdist.png')
        plt.clf()

        #pie chart for base composition
        q1_1 = self.hash_qual1_1+self.hash_qual1_2+self.hash_qual1_3+self.hash_qual1_4+self.hash_qual1_5  #0-10
        q1_2 = self.hash_qual1_6+self.hash_qual1_7+self.hash_qual1_8+self.hash_qual1_9+self.hash_qual1_10 #11-20
        q1_3 = self.hash_qual1_11+self.hash_qual1_12+self.hash_qual1_13+self.hash_qual1_14+self.hash_qual1_15 #21-30
        q1_4 = self.hash_qual1_16+self.hash_qual1_17+self.hash_qual1_18+self.hash_qual1_19+self.hash_qual1_20+self.hash_qual1_21 #31-42

        x_list = [q1_1, q1_2, q1_3, q1_4]
        label_list = ["0-10", "11-20", "21-30", "31-42"]
        plt.axis("equal")   #to make circle, default is oval
        plt.pie(x_list, labels=label_list, autopct="%1.1f%%")
        plt.title("Quality Score Distribution")
        plt.savefig(FileP1_basname+'_QualGroup.png')
        plt.clf()

        q2_1 = self.hash_qual2_1+self.hash_qual2_2+self.hash_qual2_3+self.hash_qual2_4+self.hash_qual2_5  #0-10
        q2_2 = self.hash_qual2_6+self.hash_qual2_7+self.hash_qual2_8+self.hash_qual2_9+self.hash_qual2_10 #11-20
        q2_3 = self.hash_qual2_11+self.hash_qual2_12+self.hash_qual2_13+self.hash_qual2_14+self.hash_qual2_15 #21-30
        q2_4 = self.hash_qual2_16+self.hash_qual2_17+self.hash_qual2_18+self.hash_qual2_19+self.hash_qual2_20+self.hash_qual2_21 #31-42

        x_list = [q2_1, q2_2, q2_3, q2_4]
        label_list = ["0-10", "11-20", "21-30", "31-42"]
        plt.axis("equal")   #to make circle, default is oval
        plt.pie(x_list, labels=label_list, autopct="%1.1f%%")
        plt.title("Quality Score Distribution")
        plt.savefig(FileP2_basname+'_QualGroup.png')
        plt.clf()
