#! /usr/bin/python3
# -*- coding: utf-8 -*-

__Version__ = "0.02"
__Author__ = "pzweuj"
__Date__ = "20210428"

import os
import sys
import shutil
import pandas as pd

FUN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(FUN_DIR)
from Other.function import mkdir

class Neoantigen(object):
    """
    Neoantigen新抗原预测模块，进行MHC I新抗原预测
    对应MHC I类结果，以%Rank_EL值为准
    MHC I  强： <0.5%； 弱： <2%
    MHC II 强： <2%  ； 弱： <10%
    
    IC50：
        IC50： <50nM
        被测量的拮抗剂的半抑制浓度。它能指示某一药物或者物质，
        在抑制某些生物程的半量。在凋亡方面，可理解为一定浓度的某种药物诱导肿瘤细胞凋亡50%，
        该浓度称为50%抑制浓度，即凋亡细胞与全部细胞数之比等于50%时所对应的浓度，
        IC50值可以用来衡量药物诱导凋亡的能力，即诱导能力越强，该数值越低。

    此功能速度很慢，大panel样本大概需要4小时

    多个预测算法对比文章：
    http://cancerimmunolres.aacrjournals.org/content/7/5/719.abstract
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
            HLAClassI["A2"] = "NA"
        if HLAClassI["B1"] == HLAClassI["B2"]:
            HLAClassI["B2"] = "NA"
        if HLAClassI["C1"] == HLAClassI["C2"]:
            HLAClassI["C2"] = "NA"

        for k in HLAClassI.keys():
            HLAClassI[k] = HLAClassI[k].replace("A*", "hla_a_").replace("B*", "hla_b_").replace("C*", "hla_c_").replace(":", "_")
        output = "\t".join([sampleID, HLAClassI["A1"], HLAClassI["A2"], HLAClassI["B1"], HLAClassI["B2"], HLAClassI["C1"], HLAClassI["C2"]])
        print(output)

        results = open(tmpDir + "/" + sampleID + ".hlas", "w")
        results.write(output + "\n")
        results.close()

    # NeoPredPipe
    # https://github.com/MathOnco/NeoPredPipe
    # NetMHCpan
    # http://www.cbs.dtu.dk/services/NetMHCpan/
    def neopredpipe(self):
        sampleID = self.sample
        pairID = self.pair
        resultsDir = self.output

        threads = self.threads
        buildver = self.buildver
        NeoPred = "/home/bioinfo/ubuntu/software/NeoPredPipe-1.1/NeoPredPipe.py"

        tmpDir = resultsDir + "/tempFile/Neoantigen_" + sampleID
        mkdir(tmpDir)
        vcf = resultsDir + "/vcf/" + sampleID + ".filter.vcf"
        vcfDir = tmpDir + "/vcf"
        mkdir(vcfDir)
        shutil.copy(vcf, vcfDir + "/" + sampleID + ".vcf")
        self.extractHLAResults()

        cmd = """
            python {NeoPred} \\
                -I {tmpDir}/vcf \\
                -H {tmpDir}/{sampleID}.hlas \\
                -o {tmpDir} \\
                -n {sampleID} \\
                -c 0
            sed '1iSample\\tR\\tLine\\tchrom\\tallelepos\\tref\\talt\\tGeneName\\tpos\\thla\\tpeptide\\tcore\\tOf\\tGp\\tGl\\tIp\\tIl\\tIcore\\tIdentity\\tScore_EL\\t%Rank_EL\\tScore_BA\\t%Rank_BA\\tAff(nM)\\tCandidate\\tBindLevel\\tNovelty' \\
                {tmpDir}/{sampleID}.neoantigens.txt > {tmpDir}/{sampleID}.neoantigens.tsv
            cp {tmpDir}/{sampleID}.neoantigens.tsv {resultsDir}/Neoantigen/
        """.format(NeoPred=NeoPred, vcf=vcf, tmpDir=tmpDir, sampleID=sampleID, resultsDir=resultsDir)
        print(cmd)
        os.system(cmd)

    # 结果过滤
    def neopredpipe_filter(self):
        sampleID = self.sample
        pairID = self.pair
        resultsDir = self.output

        threads = self.threads
        buildver = self.buildver

        df = pd.read_csv("Neoantigen/T22738.neoantigens.tsv", sep="\t", header=0)
        print("使用NetMHCpan发现新抗原：" + str(len(df)) + "个")

        # 过滤%Rank_EL <= 0.5, Aff(nM) <= 50, %Rank_EL升序排列
        df2 = df[(df["%Rank_EL"] <= 0.5) & (df["Aff(nM)"] <= 50)]
        df3 = df2.sort_values("%Rank_EL", ascending=True)
        df4 = df3.round({"Score_EL": 6, "Score_BA": 6, "%Rank_EL": 3, "%Rank_BA": 3})
        df4.to_csv("Neoantigen/T22738.neoantigens.filter.txt", header=True, index=None, sep="\t")

        print("过滤后：" + str(len(df4)) + "个")

# end