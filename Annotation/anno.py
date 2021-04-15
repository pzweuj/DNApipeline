#! /usr/bin/python3
# -*- coding: utf-8 -*-

__Version__ = "0.19"
__Author__ = "pzweuj"
__Date__ = "20210415"

import os
import sys

FUN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(FUN_DIR)
from Other.function import mkdir


class Annotation(object):
    """
    注释模块
    后续结果将全部改放到中间文件夹中，再使用自建脚本处理为excel表格后
    再放置到annotation文件夹中，表格格式保持一致
    注释流程不再可选比对工具，而是改用先使用snpEff注释，再使用annovar注释，最后对结果进行整理的模式
    目前未完成，更新中
    """
    def __init__(self, runningInfo):
        self.runningInfo = runningInfo
        self.sample = runningInfo["sample"]
        self.pair = runningInfo["pair"]
        self.rawdata = runningInfo["rawdata"]
        self.output = runningInfo["output"]

        self.threads = str(runningInfo["process"]["threads"])
        self.runApp = runningInfo["process"]["Annotation"]

        mkdir(self.output)
        mkdir(self.output + "/tempFile")
        mkdir(self.output + "/annotation")


    # snpEff
    # https://pcingola.github.io/SnpEff/
    def snpeff(self):
        humandb = self.runningInfo["setting"]["Annotation"]["humandb"]
        buildver = self.runningInfo["setting"]["Annotation"]["buildver"]
        resultsDir = self.output
        sampleID = self.sample
        pairID = self.pair

        if pairID == None:
            vcfFile = resultsDir + "/vcf/" + sampleID + ".vcf"
        else:
            vcfFile = resultsDir + "/vcf/" + sampleID + ".filter.vcf"        

        tmpDir = resultsDir + "/tempFile/snpeff_" + sampleID
        cmd = """
            java -jar /home/bioinfo/ubuntu/software/snpEff/snpEff.jar \\
                -c /home/bioinfo/ubuntu/software/snpEff/snpEff.config \\
                {buildver} {vcfFile} > {tmpDir}/{sampleID}.snpeff.vcf
            cp {tmpDir}/{sampleID}.snpeff.vcf {resultsDir}/annotation/
        """.format(buildver=buildver, vcfFile=vcfFile, tmpDir=tmpDir, sampleID=sampleID, resultsDir=resultsDir)
        print(cmd)
        os.system(cmd)
    
    # annovar
    # https://annovar.openbioinformatics.org/en/latest/
    def annovar(self):
        humandb = self.runningInfo["setting"]["Annotation"]["humandb"]
        buildver = self.runningInfo["setting"]["Annotation"]["buildver"]
        resultsDir = self.output
        sampleID = self.sample
        pairID = self.pair
        threads = self.threads

        tmpDir = resultsDir + "/tempFile/annovar_" + sampleID
        mkdir(tmpDir)

        self.snpeff()

        cmd = """
            convert2annovar.pl -format vcf4 \\
                {resultsDir}/annotation/{sampleID}.snpeff.vcf \\
                --includeinfo > {tmpDir}/{sampleID}.avinput
            table_annovar.pl {tmpDir}/{sampleID}.avinput \\
                {humandb} -buildver {buildver} \\
                -out {resultsDir}/annotation/{sampleID} -remove \\
                -protocol refGene,avsnp150,gnomad211_genome,clinvar_20210308,JaxCkb,Civic,OncoKB,dbnsfp41a,cosmic92_coding \\
                -operation g,f,f,f,f,f,f,f,f \\
                -nastring - -thread {threads} -otherinfo
        """.format(tmpDir=tmpDir, resultsDir=resultsDir, sampleID=sampleID, humandb=humandb, threads=threads, buildver=buildver)
        print(cmd)
        os.system(cmd)
    
    # annovar结果过滤与标准化
    def annovarResultsFilter(self):
        buildver = self.runningInfo["setting"]["Annotation"]["buildver"]
        resultsDir = self.output
        sampleID = self.sample

        annvarResultsFile = """{resultsDir}/annotation/{sampleID}.{buildver}_multianno.txt""".format(resultsDir=resultsDir, sampleID=sampleID, buildver=buildver)
        SnpOutputFile = annvarResultsFile.replace(".txt", ".snp.txt")
        IndelOutputFile = annvarResultsFile.replace(".txt", ".indel.txt")

        annvarResults = open(annvarResultsFile, "r")
        SnpOutput = open(SnpOutputFile, "w")
        IndelOutput = open(IndelOutputFile, "w")

        for line in annvarResults:
            if line.startswith("Chr"):
                lineAfterSplit = line.split("\t")
                outputList = []
                for i in range(89):
                    need = lineAfterSplit[i]
                    outputList.append(need)
                outputList.append("GT")
                outputList.append("DP")
                outputList.append("Ref_AD")
                outputList.append("Alt_AD")
                outputList.append("AF")
                outputString = "\t".join(outputList) + "\n"
                SnpOutput.write(outputString)
                IndelOutput.write(outputString)        

            else:
                lineAfterSplit = line.split("\t")
                Func = lineAfterSplit[5]
                Ref = lineAfterSplit[3]
                Alt = lineAfterSplit[4]
                ExonicFunc = lineAfterSplit[8]
                FORMAT = lineAfterSplit[97].split(":")
                FORMAT_results = lineAfterSplit[98].split(":")

                Format_list_zipped = zip(FORMAT, FORMAT_results)
                Format_list = list(Format_list_zipped)
                Format_dict = {}
                for formats in Format_list:
                    Format_dict[formats[0]] = formats[1]

                GT = Format_dict["GT"]
                DP = Format_dict["DP"]
                AD = Format_dict["AD"].split(",")
                Ref_AD = AD[0]
                Alt_AD = AD[1]
                try:
                    AF = "%.2f" % ((float(Alt_AD) / float(DP)) * 100) + "%"
                except Exception:
                    AF = "-"

                outputList = []
                for i in range(89):
                    need = lineAfterSplit[i]
                    outputList.append(need)

                outputList.append(GT)
                outputList.append(DP)
                outputList.append(Ref_AD)
                outputList.append(Alt_AD)
                outputList.append(AF)

                if ("splicing" in Func) or ("exonic" in Func):
                    if "ncRNA" not in Func:
                        if ("-" in Ref) or ("-" in Alt):
                            indel = "\t".join(outputList) + "\n"
                            IndelOutput.write(indel)                    
                        else:
                            if ("nonsynonymous" in ExonicFunc) or ("stopgain" in ExonicFunc):
                                snp = "\t".join(outputList) + "\n"
                                SnpOutput.write(snp)

        annvarResults.close()
        SnpOutput.close()
        IndelOutput.close()
        print("完成annovar结果过滤与格式调整")

