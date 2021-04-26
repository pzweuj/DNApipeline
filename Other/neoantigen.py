#! /usr/bin/python3
# -*- coding: utf-8 -*-

__Version__ = "0.01"
__Author__ = "pzweuj"
__Date__ = "20210425"

import os
import sys

FUN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(FUN_DIR)
from Other.function import mkdir

class Neoantigen(object):
    """
    Neoantigen新抗原预测模块，进行MHC I新抗原预测
    """
    def __init__(self, runningInfo):
        self.runningInfo = runningInfo
        self.sample = runningInfo["sample"]
        self.pair = runningInfo["pair"]
        self.output = runningInfo["output"]

        self.threads = str(runningInfo["process"]["threads"])
        self.buildver = runningInfo["setting"]["Annotation"]["buildver"]

        mkdir(self.output)
        mkdir(self.output + "/tempFile")
        mkdir(self.output + "/Neoantigen")

    # 提取HLA结果
    def extractHLAResults(self):
        sampleID = self.sample
        pairID = self.pair
        resultsDir = self.output
        HLApath = resultsDir + "/HLA"
        hlaFiles = os.listdir(HLApath)

        tmpDir = resultsDir + "/tempFile/Neoantigen_" + sampleID
        mkdir(tmpDir)

        hla = "-"

        f_list = []
        for f in hlaFiles:
            if pairID in f:
                f_list.append(f)

        hlaCheck = "".join(f_list)
        if "optitype" in hlaCheck:
            hla = HLApath + "/" + pairID + ".optitype.txt"
        elif "seq2HLA" in hlaCheck:
            hla = HLApath + "/" + pairID + ".seq2HLA.txt"
        elif "hlahd" in hlaCheck:
            hla = HLApath + "/" + pairID + ".hlahd.txt"
        elif "hlascan" in hlaCheck:
            hla = HLApath + "/" + pairID + ".hlascan.txt"
        else:
            print("未找到正常样本 " + pairID + " HLA结果，退出")
            exit()

        print("正在分析 " + hla + " 文件")
        HLAClassI = {}
        if "optitype" in hla:
            optitype = open(hla, "r")
            for line in optitype:
                if not "Objective" in line:
                    if line != "\n":
                        lines = line.split("\t")
                        HLAClassI["A1"] = lines[1]
                        HLAClassI["A2"] = lines[2]
                        HLAClassI["B1"] = lines[3]
                        HLAClassI["B2"] = lines[4]
                        HLAClassI["C1"] = lines[5]
                        HLAClassI["C2"] = lines[6]
            optitype.close()

        elif "seq2HLA" in hla:
            seq2hla = open(hla, "r")
            for line in seq2hla:
                if not line.startswith("#"):
                    lines = line.split("\t")
                    if lines[0] == "A":
                        HLAClassI["A1"] = lines[1].replace("'", "")
                        HLAClassI["A2"] = lines[3].replace("'", "")
                    elif lines[0] == "B":
                        HLAClassI["B1"] = lines[1].replace("'", "")
                        HLAClassI["B2"] = lines[3].replace("'", "")
                    elif lines[0] == "C":
                        HLAClassI["C1"] = lines[1].replace("'", "")
                        HLAClassI["C2"] = lines[3].replace("'", "")
                    else:
                        continue
            seq2hla.close()

        elif "hlahd" in hla:
            hlahd = open(hla, "r")
            for line in hlahd:
                lines = line.replace("\n", "").split("\t")
                if lines[0] == "A":
                    HLAClassI["A1"] = lines[1].split("-")[1][:-3]
                    HLAClassI["A2"] = lines[2].split("-")[1][:-3]
                elif lines[0] == "B":
                    HLAClassI["B1"] = lines[1].split("-")[1][:-3]
                    HLAClassI["B2"] = lines[2].split("-")[1][:-3]
                elif lines[0] == "C":
                    HLAClassI["C1"] = lines[1].split("-")[1][:-3]
                    HLAClassI["C2"] = lines[2].split("-")[1][:-3]
                else:
                    continue
            hlahd.close()

        elif "hlascan" in hla:
            hlascan = open(hla, "r")
            l = []
            for line in hlascan:
                if line.startswith("["):
                    t = line.split("\t")[1][0:5]
                    l.append(t)
            HLAClassI["A1"] = "A*" + l[0]
            HLAClassI["A2"] = "A*" + l[1]
            HLAClassI["B1"] = "B*" + l[2]
            HLAClassI["B2"] = "B*" + l[3]
            HLAClassI["C1"] = "C*" + l[4]
            HLAClassI["C2"] = "C*" + l[5]
            hlascan.close()

        else:
            print("未找到正常样本 " + pairID + " HLA结果，退出")
            exit()

        print(HLAClassI)

        ## 格式处理例子
        """
        Patient HLA-A_1 HLA-A_2 HLA-B_1 HLA-B_2 HLA-C_1 HLA-C_2
        'NA' is used when the HLA typing predicts the same HLA subtype for A, B, or C.
        """
        if HLAClassI["A1"] == HLAClassI["A2"]:
            HLAClassI["A2"] == "NA"
        if HLAClassI["B1"] == HLAClassI["B2"]:
            HLAClassI["B2"] == "NA"
        if HLAClassI["C1"] == HLAClassI["C2"]:
            HLAClassI["C2"] == "NA"

        for k in HLAClassI.keys():
            HLAClassI[k] = HLAClassI[k].replace("A*", "hla_a_").replace("B*", "hla_b_").replace("C*", "hla_c_").replace(":", "_")
        output = "\t".join([sampleID, HLAClassI["A1"], HLAClassI["A2"], HLAClassI["B1"], HLAClassI["B2"], HLAClassI["C1"], HLAClassI["C2"]])
        print(output)

        results = open(tmpDir + "/" + sampleID + ".hlas", "w")
        results.write(output + "\n")
        results.close()

    # NeoPredPipe
    # https://github.com/MathOnco/NeoPredPipe
    def neopredpipe(self):
        sampleID = self.sample
        pairID = self.pair
        resultsDir = self.output

        threads = self.threads
        buildver = self.buildver

        tmpDir = resultsDir + "/tempFile/Neoantigen_" + sampleID
        mkdir(tmpDir)
        vcf = resultsDir + "/vcf/" + sampleID + ".filter.vcf"
        vcfDir = tmpDir + "/vcf"
        mkdir(vcfDir)
        shutil.copy(vcf, vcfDir + "/" + sampleID + ".vcf")
        NeoPred = "/home/bioinfo/ubuntu/software/NeoPredPipe-1.1/NeoPredPipe.py"

        cmd = """
            python {NeoPred} \\
                -I {tmpDir}/vcf \\
                -H {tmpDir}/{sampleID}.hlas \\
                -o {tmpDir} \\
                -n {sampleID} \\
                -c 0
        """.format(NeoPred=NeoPred, vcf=vcf, tmpDir=tmpDir, sampleID=sampleID)
        print(cmd)
        os.system(cmd)

    # 结果过滤
    def neopredpipe_filter(self):
        sampleID = self.sample
        pairID = self.pair
        resultsDir = self.output

        threads = self.threads
        buildver = self.buildver

        

