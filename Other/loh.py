#! /usr/bin/python3
# -*- coding: utf-8 -*-

__Version__ = "0.02"
__Author__ = "pzweuj"
__Date__ = "20210421"

import os
import sys

FUN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(FUN_DIR)
from Other.function import mkdir


class LOH(object):
    """
    LOH检测模块，用于检测HLA区域杂合性缺失
    """

    def __init__(self, runningInfo):
        self.runningInfo = runningInfo
        self.sample = runningInfo["sample"]
        self.rawdata = runningInfo["rawdata"]
        self.output = runningInfo["output"]

        self.threads = str(runningInfo["process"]["threads"])
        self.buildver = runningInfo["setting"]["Annotation"]["buildver"]

        mkdir(self.output)
        mkdir(self.output + "/tempFile")
        mkdir(self.output + "/HLA")


    # 提取HLA结果
    def extractHLAResults(self):
        sampleID = self.sample
        pairID = self.pair
        resultsDir = self.output
        HLApath = resultsDir + "/HLA"
        hlaFiles = os.listdir(HLApath)

        tmpDir = resultsDir + "/tempFile/LOH_" + sampleID
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
        HLA-A    hla_a_24_02    hla_a_69_01
        HLA-B    hla_b_15_01    hla_b_52_01
        HLA-C    hla_c_01_02    hla_c_12_02
        """

        output1 = "HLA-A\thla_a_" + HLAClassI["A1"].split("*")[1].replace(":", "_") + "\thla_a_" + HLAClassI["A2"].split("*")[1].replace(":", "_")
        output2 = "HLA-B\thla_b_" + HLAClassI["B1"].split("*")[1].replace(":", "_") + "\thla_b_" + HLAClassI["B2"].split("*")[1].replace(":", "_")
        output3 = "HLA-C\thla_c_" + HLAClassI["C1"].split("*")[1].replace(":", "_") + "\thla_c_" + HLAClassI["C2"].split("*")[1].replace(":", "_")

        print(output1)
        print(output2)
        print(output3)

        results = open(tmpDir + "/" + sampleID + ".hlas", "w")
        results.write(output1 + "\n")
        results.write(output2 + "\n")
        results.write(output3 + "\n")
        results.close()


    # LOHHLA
    # https://github.com/slagtermaarten/LOHHLA
    # https://bitbucket.org/mcgranahanlab/lohhla/src/master/
    # 与原版对比增加了hg38支持以及修正部分bug ---> slagtermaarten
    # 将picard部分更新为新picard ---> pzw
    # 参考https://github.com/pyc1216/hlaloh-pipeline
    def lohhla(self):
        self.extractHLAResults()
        pass



    # PyLOH
    # https://github.com/uci-cbcl/PyLOH
    def pyloh(self):
        self.extractHLAResults()
        pass

