#! /usr/bin/python3
# -*- coding: utf-8 -*-

__Version__ = "0.03"
__Author__ = "pzweuj"
__Date__ = "20210425"

import os
import sys
import math
import shutil

FUN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(FUN_DIR)
from Other.function import mkdir


class LOH(object):
    """
    LOH检测模块，用于检测HLA区域杂合性缺失
    需求HLA-A、HLA-B、HLA-C三个class I基因的分型结果，因此仅能在大panel或WES等项目中使用
    ！！！！未验证！！！！
    """

    def __init__(self, runningInfo):
        self.runningInfo = runningInfo
        self.sample = runningInfo["sample"]
        self.pair = runningInfo["pair"]
        self.rawdata = runningInfo["rawdata"]
        self.output = runningInfo["output"]

        self.threads = str(runningInfo["process"]["threads"])
        self.buildver = runningInfo["setting"]["Annotation"]["buildver"]

        mkdir(self.output)
        mkdir(self.output + "/tempFile")
        mkdir(self.output + "/HLA")
        mkdir(self.output + "/LOH")


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
        try:
            output1 = "hla_a_" + HLAClassI["A1"].split("*")[1].replace(":", "_") + "\nhla_a_" + HLAClassI["A2"].split("*")[1].replace(":", "_")
        except:
            output1 = "\n"
            print("无法分型HLA-A")
        
        try:
            output2 = "hla_b_" + HLAClassI["B1"].split("*")[1].replace(":", "_") + "\nhla_b_" + HLAClassI["B2"].split("*")[1].replace(":", "_")
        except:
            output2 = "\n"
            print("无法分型HLA-B")

        try:
            output3 = "hla_c_" + HLAClassI["C1"].split("*")[1].replace(":", "_") + "\nhla_c_" + HLAClassI["C2"].split("*")[1].replace(":", "_")
        except:
            output3 = "\n"
            print("无法分型HLA-C")

        print(output1)
        print(output2)
        print(output3)

        results = open(tmpDir + "/" + sampleID + ".hlas", "w")
        results.write(output1 + "\n")
        results.write(output2 + "\n")
        results.write(output3 + "\n")
        results.close()
        cmd = """
            sed -i '/^$/d' {tmpDir}/{sampleID}.hlas
        """.format(tmpDir=tmpDir, sampleID=sampleID)
        os.system(cmd)

    # 评估肿瘤倍性与肿瘤纯度
    # PyLOH
    # https://github.com/uci-cbcl/PyLOH
    def pyloh(self):
        pass

    # THetA
    # https://github.com/raphael-group/THetA
    def theta(self):
        pass

    # PureCN
    # https://github.com/lima1/PureCN
    def purecn(self):
        pass


    # 暂时先模拟一个用着，肿瘤纯度应该可以直接从病理获得，比通过软件评估准确
    def TempPurityPloidy(self):
        sampleID = self.sample
        pairID = self.pair
        resultsDir = self.output
        
        tmpDir = resultsDir + "/tempFile/LOH_" + sampleID
        mkdir(tmpDir)

        solution = open(tmpDir + "/" + sampleID + ".solutions.txt", "w")
        solution.write("Ploidy\ttumorPurity\ttumorPloidy\n")
        solution.write(sampleID + "\t2\t0.8\t2\n")
        solution.close()


    # LOHHLA
    # https://github.com/slagtermaarten/LOHHLA
    # https://bitbucket.org/mcgranahanlab/lohhla
    # 与原版对比增加了hg38支持以及修正部分bug ---> slagtermaarten
    # 将picard部分更新为picard 2.25.2 ---> pzw
    def lohhla(self):
        sampleID = self.sample
        pairID = self.pair
        resultsDir = self.output

        tmpDir = resultsDir + "/tempFile/LOH_" + sampleID
        mkdir(tmpDir)

        # 对HLA结果进行提取
        self.extractHLAResults()
        # 肿瘤纯度与倍性评估
        self.TempPurityPloidy()

        # HLAfasta来源
        # https://github.com/jason-weirather/hla-polysolver/blob/master/data/abc_complete.fasta
        HLAloc = "/home/bioinfo/ubuntu/software/LOHHLA"
        HLAfasta = "/home/bioinfo/ubuntu/software/LOHHLA/data/abc_complete.fasta"
        HLAexon = "/home/bioinfo/ubuntu/software/LOHHLA/data/hla.dat"
        cmd = """
            Rscript {HLAloc}/LOHHLAscript.R \\
                --patientId {sampleID} \\
                --outputDir {tmpDir} \\
                --normalBAMfile {resultsDir}/bam/{pairID}.bam \\
                --tumorBAMfile {resultsDir}/bam/{sampleID}.bam \\
                --BAMDir {resultsDir}/bam \\
                --hlaPath {tmpDir}/{sampleID}.hlas \\
                --HLAfastaLoc {HLAfasta} \\
                --CopyNumLoc {tmpDir}/{sampleID}.solutions.txt \\
                --mappingStep TRUE \\
                --minCoverageFilter 10 \\
                --fishingStep TRUE \\
                --cleanUp FALSE \\
                --gatkDir /home/bioinfo/ubuntu/software/picard \\
                --novoDir /home/bioinfo/ubuntu/software/novocraft \\
                --LOHHLA_loc {HLAloc} \\
                --HLAexonLoc {HLAexon}
        """.format(HLAloc=HLAloc, sampleID=sampleID, tmpDir=tmpDir, resultsDir=resultsDir, pairID=pairID, HLAfasta=HLAfasta, HLAexon=HLAexon)
        print(cmd)
        os.system(cmd)

        # results
        f_list = os.listdir(tmpDir)
        lohfile = "-"
        for f in f_list:
            if "DNA.HLAlossPrediction_" in f:
                lohfile = f
        if lohfile == "-":
            print("未找到结果")
            exit()
        else:
            lohf = open(tmpDir + "/" + lohfile, "r")
            loh_results = open(tmpDir + "/" + sampleID + ".loh.txt", "w")
            loh_results.write("HLAType\tHLACopyNumWithBAFBin\tpval\tLOHstat\n")
            for line in lohf:
                if not line.startswith("message"):
                    lines = line.split("\t")
                    HLAtype = lines[1][0:5]
                    if lines[0].startswith("homozygous_alleles"):
                        HLAtype2copyNumWithBAFBin = "-"
                        pVal = "-"
                        LOHstat = "FALSE"
                    else:
                        HLAtype2copyNumWithBAFBin = lines[28]
                        pVal = lines[33]
                        if pVal == "NA":
                            LOHstat = "-"
                        # 参考 https://github.com/pyc1216/hlaloh-pipeline/blob/master/scripts/get_result.py 获得回报结果
                        elif float(pVal) < 0.01 and float(HLAtype2copyNumWithBAFBin) < math.log2(0.5):
                            LOHstat = "TRUE"
                        else:
                            LOHstat = "FALSE"
                    loh_results.write(HLAtype + "\t" + HLAtype2copyNumWithBAFBin + "\t" + pVal + "\t" + LOHstat + "\n")
            loh_results.close()
            lohf.close()
        shutil.copy(tmpDir + "/" + sampleID + ".loh.txt", resultsDir + "/LOH/" + sampleID + ".loh.txt")
        print("LOH分析完成")




