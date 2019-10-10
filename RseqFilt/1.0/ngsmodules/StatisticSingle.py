#!/usr/bin/python

import csv
import os
from matplotlib import pyplot as plt
import numpy as np


class StatisticSingle:

    RawReadCount1 = 0
    TrimReadCount1 = 0
    ShortReadCountt1 = 0
    TotGCSum1 = 0
    TotLenSum1 = 0
    TotGCSuma1 = 0
    TotLenSuma1 = 0
    CleanReadCount1 = 0
    NreadCount1 = 0         ##Reads containing more than Nbase % (given)
    AdtrimCount1 = 0
    FailQualCount1 = 0
    TotACount1 = 0
    TotTCount1 = 0
    TotGCount1 = 0
    TotCCount1 = 0
    TotNCount1 = 0
    TotAaCount1 = 0
    TotTaCount1 = 0
    TotGaCount1 = 0
    TotCaCount1 = 0
    TotNaCount1 = 0
    TotQual1 = 0
    TotQual2 = 0
    TotQual3 = 0
    TotQual4 = 0
    TotTrimLen1 = 0
    hash_gcList1_1 = 0
    hash_gcList1_2 = 0
    hash_gcList1_3 = 0
    hash_gcList1_4 = 0
    hash_gcList1_5 = 0
    hash_gcList1_6 = 0
    hash_gcList1_7 = 0
    hash_gcList1_8 = 0
    hash_gcList1_9 = 0
    hash_gcList1_10 = 0
    hash_gcList1a_1 = 0
    hash_gcList1a_2 = 0
    hash_gcList1a_3 = 0
    hash_gcList1a_4 = 0
    hash_gcList1a_5 = 0
    hash_gcList1a_6 = 0
    hash_gcList1a_7 = 0
    hash_gcList1a_8 = 0
    hash_gcList1a_9 = 0
    hash_gcList1a_10 = 0
    hash_qual1_1 = 0
    hash_qual1_2 = 0
    hash_qual1_3 = 0
    hash_qual1_4 = 0
    hash_qual1_5 = 0
    hash_qual1_6 = 0
    hash_qual1_7 = 0
    hash_qual1_8 = 0
    hash_qual1_9 = 0
    hash_qual1_10 = 0
    hash_qual1_11 = 0
    hash_qual1_12 = 0
    hash_qual1_13 = 0
    hash_qual1_14 = 0
    hash_qual1_15 = 0
    hash_qual1_16 = 0
    hash_qual1_17 = 0
    hash_qual1_18 = 0
    hash_qual1_19 = 0
    hash_qual1_20 = 0
    hash_qual1_21 = 0
    hash_qual1a_1 = 0
    hash_qual1a_2 = 0
    hash_qual1a_3 = 0
    hash_qual1a_4 = 0
    hash_qual1a_5 = 0
    hash_qual1a_6 = 0
    hash_qual1a_7 = 0
    hash_qual1a_8 = 0
    hash_qual1a_9 = 0
    hash_qual1a_10 = 0
    hash_qual1a_11 = 0
    hash_qual1a_12 = 0
    hash_qual1a_13 = 0
    hash_qual1a_14 = 0
    hash_qual1a_15 = 0
    hash_qual1a_16 = 0
    hash_qual1a_17 = 0
    hash_qual1a_18 = 0
    hash_qual1a_19 = 0
    hash_qual1a_20 = 0
    hash_qual1a_21 = 0
    TotQual1Sum = 0
    TotQual1aSum = 0
    NcontReads1 = 0
    NcontReads1a = 0
    TrimQualCt = 0
    LowQualCount = 0
    LenList = []
    LenLista = []

    def __init__(self, InFile):
        InFile = open(InFile, 'r')
        InFile = csv.reader(InFile, delimiter='\t')
        for rec in InFile:
            self.RawReadCount1 += int(rec[0])
            self.TrimReadCount1 += int(rec[1])
            self.ShortReadCountt1 += int(rec[2])
            self.TotGCSum1 += float(rec[4])
            self.TotLenSum1 += float(rec[5])
            self.TotGCSuma1 += float(rec[7])
            self.TotLenSuma1 += float(rec[9])
            self.CleanReadCount1 += int(rec[11])
            self.NreadCount1 += int(rec[12])         ##Reads containing more than Nbase % (given)
            self.AdtrimCount1 += int(rec[14])
            self.FailQualCount1 += int(rec[17])
            self.TotACount1 += int(rec[19])
            self.TotTCount1 += int(rec[21])
            self.TotGCount1 += int(rec[23])
            self.TotCCount1 += int(rec[25])
            self.TotNCount1 += int(rec[27])
            self.TotAaCount1 += int(rec[29])
            self.TotTaCount1 += int(rec[31])
            self.TotGaCount1 += int(rec[33])
            self.TotCaCount1 += int(rec[35])
            self.TotNaCount1 += int(rec[37])
            # self.TotQual1 += int(rec[39])
            # self.TotQual2 += int(rec[40])
            # self.TotQual3 += int(rec[41])
            # self.TotQual4 += int(rec[42])
            self.TotTrimLen1 += int(rec[54])
            self.hash_gcList1_1 += int(rec[56])
            self.hash_gcList1_2 += int(rec[57])
            self.hash_gcList1_3 += int(rec[58])
            self.hash_gcList1_4 += int(rec[59])
            self.hash_gcList1_5 += int(rec[60])
            self.hash_gcList1_6 += int(rec[61])
            self.hash_gcList1_7 += int(rec[62])
            self.hash_gcList1_8 += int(rec[63])
            self.hash_gcList1_9 += int(rec[64])
            self.hash_gcList1_10 += int(rec[65])
            self.hash_gcList1a_1 += int(rec[77])
            self.hash_gcList1a_2 += int(rec[78])
            self.hash_gcList1a_3 += int(rec[79])
            self.hash_gcList1a_4 += int(rec[80])
            self.hash_gcList1a_5 += int(rec[81])
            self.hash_gcList1a_6 += int(rec[82])
            self.hash_gcList1a_7 += int(rec[83])
            self.hash_gcList1a_8 += int(rec[84])
            self.hash_gcList1a_9 += int(rec[85])
            self.hash_gcList1a_10 += int(rec[86])
            self.TotQual1Sum += float(rec[98])
            self.TotQual1aSum += float(rec[99])
            self.NcontReads1 += int(rec[102])
            self.NcontReads1a += int(rec[103])
            self.TrimQualCt += int(rec[104])
            self.hash_qual1_1 += int(rec[105])
            self.hash_qual1_2 += int(rec[106])
            self.hash_qual1_3 += int(rec[107])
            self.hash_qual1_4 += int(rec[108])
            self.hash_qual1_5 += int(rec[109])
            self.hash_qual1_6 += int(rec[110])
            self.hash_qual1_7 += int(rec[111])
            self.hash_qual1_8 += int(rec[112])
            self.hash_qual1_9 += int(rec[113])
            self.hash_qual1_10 += int(rec[114])
            self.hash_qual1_11 += int(rec[115])
            self.hash_qual1_12 += int(rec[116])
            self.hash_qual1_13 += int(rec[117])
            self.hash_qual1_14 += int(rec[118])
            self.hash_qual1_15 += int(rec[119])
            self.hash_qual1_16 += int(rec[120])
            self.hash_qual1_17 += int(rec[121])
            self.hash_qual1_18 += int(rec[122])
            self.hash_qual1_19 += int(rec[123])
            self.hash_qual1_20 += int(rec[124])
            self.hash_qual1_21 += int(rec[125])
            self.hash_qual1a_1 += int(rec[126])
            self.hash_qual1a_2 += int(rec[127])
            self.hash_qual1a_3 += int(rec[128])
            self.hash_qual1a_4 += int(rec[129])
            self.hash_qual1a_5 += int(rec[130])
            self.hash_qual1a_6 += int(rec[131])
            self.hash_qual1a_7 += int(rec[132])
            self.hash_qual1a_8 += int(rec[133])
            self.hash_qual1a_9 += int(rec[134])
            self.hash_qual1a_10 += int(rec[135])
            self.hash_qual1a_11 += int(rec[136])
            self.hash_qual1a_12 += int(rec[137])
            self.hash_qual1a_13 += int(rec[138])
            self.hash_qual1a_14 += int(rec[139])
            self.hash_qual1a_15 += int(rec[140])
            self.hash_qual1a_16 += int(rec[141])
            self.hash_qual1a_17 += int(rec[142])
            self.hash_qual1a_18 += int(rec[143])
            self.hash_qual1a_19 += int(rec[144])
            self.hash_qual1a_20 += int(rec[145])
            self.hash_qual1a_21 += int(rec[146])
            self.LenList.append(float(rec[43]))
            self.LenList.append(float(rec[44]))
            self.LenLista.append(float(rec[51]))
            self.LenLista.append(float(rec[52]))
#        print self.RawReadCount1

    def stat_single(self, Nbase, QualThresh, MinSize, Adapter, Trim, MinLen, File1, QualFormat, OutputDir):
        os.chdir(OutputDir)
        TotBases = self.TotACount1+self.TotTCount1+self.TotGCount1+self.TotCCount1
        TotBasesWithN = self.TotACount1+self.TotTCount1+self.TotGCount1+self.TotCCount1+self.TotNCount1
        TotBasesa = self.TotAaCount1+self.TotTaCount1+self.TotGaCount1+self.TotCaCount1
        TotBasesaWithN = self.TotAaCount1+self.TotTaCount1+self.TotGaCount1+self.TotCaCount1+self.TotNaCount1
#        print self.TotAaCount1, self.TotTaCount1, self.TotGaCount1, self.TotCaCount1,  self.TotNaCount1
        MinLenRead = min(self.LenList)
        MaxLenRead = max(self.LenList)
        try:
            self.LenLista.remove(0)
            self.LenLista.remove(1000000)
        except:
            pass
        MinLenReada = min(self.LenLista)
        MaxLenReada = max(self.LenLista)
        Q30 = self.hash_qual1_16+self.hash_qual1_17+self.hash_qual1_18+self.hash_qual1_19+self.hash_qual1_20+ \
              self.hash_qual1_21
        Q30a = self.hash_qual1a_16+self.hash_qual1a_17+self.hash_qual1a_18+self.hash_qual1a_19+self.hash_qual1a_20+ \
                self.hash_qual1a_21
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
            StatFile.write("%-50s\t%d\n" % ("Read with minimum length to keep after trimming", int(MinLen)))
        else:
            StatFile.write("%-50s\t%-s\n" % ("Filtering mode for quality value", "Filter"))
        StatFile.write("\n\n")
        StatFile.write("%-s\n" % "Filtering Statistics for given Single end files")
        StatFile.write("=============================================================================================="
                       "====================================\n\n")
        StatFile.write("%-65s\t%-20s\n" % ("Input file", File1))
        if QualFormat == 1:
            StatFile.write("%-65s\t%-20s\n" % ("Quality format", "illumina 1.8+"))
        elif QualFormat == 2:
            StatFile.write("%-65s\t%s\n" % ("Quality format", "illumina 1.6+"))
        elif QualFormat == 3:
            StatFile.write("%-65s\t%s\n" % ("Quality format", "Sanger"))
        StatFile.write("%-65s\t%-20s\n" % ("Total number of reads analyzed", self.RawReadCount1))
        StatFile.write("%-65s\t%-20s\n" % ("Total Bases (ATGC) in unfiltered reads (bp)", TotBases))
        StatFile.write("%-65s\t%-20s\n" % ("Total Bases (ATGC) in filtered reads (bp)", TotBasesa))
        StatFile.write("%-65s\t%-20s\n" % ("Minimum size (bp) of unfiltered reads", MinLenRead))
        StatFile.write("%-65s\t%-20s\n" % ("Maximum size (bp) of unfiltered reads", MaxLenRead))
        StatFile.write("%-65s\t%-20s\n" % ("Mean size (bp) of unfiltered reads",
                                           float(self.TotLenSum1/self.RawReadCount1)))
        StatFile.write("%-65s\t%.2f\n" % ("Average Quality value for unfiltered reads",
                                           float(self.TotQual1Sum/TotBasesWithN)))
        StatFile.write("%-65s\t%.2f\n" % ("Average Quality value for filtered reads",
                                           float(self.TotQual1aSum/TotBasesaWithN)))
        if MinSize:
            StatFile.write("%-65s\t%-20s\n" % ("Total number of reads below minimum size"+"("+MinSize+")",
                                               self.ShortReadCountt1))
        StatFile.write("%-65s\t%-20s\n" % ("Total unfiltered reads containing at least one uncalled base(N)",
                                           self.NcontReads1))
        StatFile.write("%-65s\t%-20s\n" % ("Total filtered reads containing at least one uncalled base(N)",
                                           self.NcontReads1a))
        if Nbase:
            StatFile.write("%-65s\t%-20s\n" % ("Reads filtered out with more than given % of N", self.NreadCount1))
        if Trim is None:
            StatFile.write("%-65s\t%-20s\n" % ("Reads filtered out for quality", self.FailQualCount1))
        if Adapter:
            StatFile.write("%-65s\t%-20s\n" % ("Reads trimmed for adapter", self.AdtrimCount1))
        StatFile.write("%-65s\t%.2f\n" % ("Average %GC content in unfiltered reads",
                                          float((self.TotGCSum1*100)/self.TotLenSum1)))
        StatFile.write("%-65s\t%.2f\n" % ("Average %GC content in filtered reads",
                                          float((self.TotGCSuma1*100)/self.TotLenSuma1)))
        if Trim:
            StatFile.write("%-65s\t%-20s\n" % ("Reads trimmed for quality", self.TrimQualCt))
        StatFile.write("%-65s\t%-20s\n" % ("Total unfiltered reads with >Q30 mean quality value", Q30))
        StatFile.write("%-65s\t%-20s\n" % ("Total filtered reads with >Q30 mean quality value ", Q30a))
        StatFile.write("%-65s\t%-20s\n" % ("Total number of reads kept", self.CleanReadCount1))
        if Trim is None:
            # StatFile.write("%-65s\t%-20s\n" % ("Minimum size of filtered reads", MinLenRead))
            # StatFile.write("%-65s\t%-20s\n" % ("Maximum size of filtered reads", MaxLenRead))
            StatFile.write("%-65s\t%-20s\n" % ("Mean size (bp) of filtered reads",
                                               float(self.TotLenSuma1/self.CleanReadCount1)))
        else:
            StatFile.write("%-65s\t%-20s\n" % ("Minimum size (bp) of filtered reads", MinLenReada))
            StatFile.write("%-65s\t%-20s\n" % ("Maximum size (bp) of filtered reads", MaxLenReada))
            StatFile.write("%-65s\t%-20s\n" % ("Mean size (bp) of filtered reads",
                                               float(self.TotLenSuma1/self.CleanReadCount1)))

        Remain = self.RawReadCount1 - self.CleanReadCount1
        StatFile.write("\n\n\nThe total number of sequences removed:"+str(Remain)+"\n")
        StatFile.close()

    def stat_vis(self, FileS_basname):

        Gcbin = range(10, 101, 10)
        Gcvalue1 = []
        Gcvalue1a = []

        Gcvalue1.extend([self.hash_gcList1_1, self.hash_gcList1_2, self.hash_gcList1_3, self.hash_gcList1_4,
                         self.hash_gcList1_5, self.hash_gcList1_6])
        Gcvalue1.extend([self.hash_gcList1_7, self.hash_gcList1_8, self.hash_gcList1_9, self.hash_gcList1_10])

        Gcvalue1a.extend([self.hash_gcList1a_1, self.hash_gcList1a_2, self.hash_gcList1a_3, self.hash_gcList1a_4,
                          self.hash_gcList1a_5, self.hash_gcList1a_6])
        Gcvalue1a.extend([self.hash_gcList1a_7, self.hash_gcList1a_8, self.hash_gcList1a_9, self.hash_gcList1a_10])

        plt.figure(1)
        plt.plot(Gcbin, Gcvalue1a, 'g', Gcbin, Gcvalue1, 'r')
        plt.xticks([10, 20, 30, 40, 50, 60, 70, 80, 90, 100], rotation='horizontal')
        plt.legend(["Filtered", "Unfiltered"], loc=2)
        plt.ylabel("# Reads", fontweight='bold')
        plt.xlabel("%GC Content", fontweight='bold')
        plt.savefig(FileS_basname+'_GCdist.png')
        plt.clf()

        # for base composition
        Bcbin = np.arange(5)
        Bcvalue1 = []
        Bcvalue1a = []

        Bcvalue1.extend([self.TotACount1, self.TotTCount1, self.TotGCount1, self.TotCCount1, self.TotNCount1])
        Bcvalue1a.extend([self.TotAaCount1, self.TotTaCount1, self.TotGaCount1, self.TotCaCount1, self.TotNaCount1])
#       print self.TotNCount1, self.TotNaCount1

        wid = 0.25
        fig = plt.figure()
        plt.bar(Bcbin, Bcvalue1a, wid, color='#6A5ACD', label='Filter', align='center')
        plt.bar(Bcbin + wid, Bcvalue1, wid, color='#CD5C5C', label='Unfilter', align='center')
        plt.ylabel('# Bases',fontweight='bold')
        plt.xticks([0, 1, 2, 3, 4], ['A', 'T', 'G', 'C', 'N'], fontsize=12, ha='center')
        plt.legend(bbox_to_anchor=(0., 1., 1., .102), loc=3, prop={'size':15}, ncol=2, mode="expand")
        #space in x axis
        plt.margins(0.05, None)
        plt.savefig(FileS_basname+'_Basedist.png')
        plt.clf()

#       #for quality
        Qualbin = range(2, 43, 2)
        Qualvalue1 = []
        Qualvalue1a = []

        Qualvalue1.extend([self.hash_qual1_1, self.hash_qual1_2, self.hash_qual1_3, self.hash_qual1_4,
                           self.hash_qual1_5, self.hash_qual1_6])
        Qualvalue1.extend([self.hash_qual1_7, self.hash_qual1_8, self.hash_qual1_9, self.hash_qual1_10,
                           self.hash_qual1_11, self.hash_qual1_12])
        Qualvalue1.extend([self.hash_qual1_13, self.hash_qual1_14, self.hash_qual1_15, self.hash_qual1_16,
                           self.hash_qual1_17, self.hash_qual1_18])
        Qualvalue1.extend([self.hash_qual1_19, self.hash_qual1_20, self.hash_qual1_21])

        Qualvalue1a.extend([self.hash_qual1a_1, self.hash_qual1a_2, self.hash_qual1a_3, self.hash_qual1a_4,
                            self.hash_qual1a_5, self.hash_qual1a_6])
        Qualvalue1a.extend([self.hash_qual1a_7, self.hash_qual1a_8, self.hash_qual1a_9, self.hash_qual1a_10,
                            self.hash_qual1a_11, self.hash_qual1a_12])
        Qualvalue1a.extend([self.hash_qual1a_13, self.hash_qual1a_14, self.hash_qual1a_15, self.hash_qual1a_16,
                            self.hash_qual1a_17, self.hash_qual1a_18])
        Qualvalue1a.extend([self.hash_qual1a_19, self.hash_qual1a_20, self.hash_qual1a_21])

        #for plotting quality comparison between between filtered and unfiltered data
        plt.figure(5)
        plt.plot(Qualbin, Qualvalue1a, 'g', Qualbin, Qualvalue1, 'r')
        plt.xticks([0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42], rotation='vertical')
        plt.legend(["Filtered", "Unfiltered"], loc=2)
        plt.ylabel("# Reads", fontweight='bold')
        plt.savefig(FileS_basname+'_Qualdist.png')
        plt.clf()

        #pie chart for base composition
        #0-10
        q1_1 = self.hash_qual1_1+self.hash_qual1_2+self.hash_qual1_3+self.hash_qual1_4+self.hash_qual1_5
        #11-20
        q1_2 = self.hash_qual1_6+self.hash_qual1_7+self.hash_qual1_8+self.hash_qual1_9+self.hash_qual1_10
        #21-30
        q1_3 = self.hash_qual1_11+self.hash_qual1_12+self.hash_qual1_13+self.hash_qual1_14+self.hash_qual1_15
        #31-42
        q1_4 = self.hash_qual1_16+self.hash_qual1_17+self.hash_qual1_18+self.hash_qual1_19+self.hash_qual1_20 \
                +self.hash_qual1_21

        x_list = [q1_1, q1_2, q1_3, q1_4]
        label_list = ["0-10", "11-20", "21-30", "31-42"]
        #to make circle, default is oval
        plt.axis("equal")
        plt.pie(x_list, labels=label_list, autopct="%1.1f%%")
        plt.title("Quality Score Distribution")
        plt.savefig(FileS_basname+'_QualGroup.png')
        plt.clf()

