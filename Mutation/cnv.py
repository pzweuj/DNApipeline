#! /usr/bin/python3
# -*- coding: utf-8 -*-

__Version__ = "0.13"
__Author__ = "pzweuj"
__Date__ = "20210224"

import os
import sys

FUN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(FUN_DIR)
from Other.function import mkdir


class CNV(object):
    """
    拷贝数变异检测模块
    """
    def __init__(self, runningInfo):
        self.runningInfo = runningInfo
        self.sample = runningInfo["sample"]
        self.rawdata = runningInfo["rawdata"]
        self.output = runningInfo["output"]

        self.threads = str(runningInfo["process"]["threads"])
        self.runApp = runningInfo["process"]["Mutation"]["CNV"]

        mkdir(self.output)
        mkdir(self.output + "/tempFile")
        mkdir(self.output + "/cnv")


    # cnvkit
    # https://cnvkit.readthedocs.io/en/stable/pipeline.html
    def cnvkit(self):
        resultsDir = self.output
        sampleID = self.sample
        referenceCnn = self.runningInfo["setting"]["Mutation"]["CNN"]
        threads = self.threads

        tmpDir = resultsDir + "/tempFile/cnvkit_" + sampleID
        mkdir(tmpDir)

        cmd = """
            cnvkit.py batch \\
                {resultsDir}/bam/{sampleID}.bam \\
                -r {referenceCnn} \\
                -d {tmpDir} \\
                -m hybrid \\
                --diagram --scatter -p {threads}
            cp {tmpDir}/{sampleID}.call.cns {resultsDir}/cnv/{sampleID}.cnvkit.cns
        """.format(tmpDir=tmpDir, resultsDir=resultsDir, sampleID=sampleID, referenceCnn=referenceCnn, threads=threads)
        print(cmd)
        os.system(cmd)

    def cnvkit_filter(self):
        # 需要建立hg19坐标与外显子编号对应数据库，落实是哪个外显子的扩增
        resultsDir = self.output
        sampleID = self.sample
        tmpDir = resultsDir + "/tempFile/cnvkit_" + sampleID

        cnvkitResultsFile = open(tmpDir + "/" + sampleID + ".call.cns", "r")
        cnvkitResults = open(resultsDir + "/cnv/" + sampleID + ".cnvkit.txt", "w")
        cnvkitResults.write("chromsome\tstart\tend\tgene\tVAF\tlog2\tdepth\tp_ttest\tprobes\tweight\n")
        for line in cnvkitResultsFile:
            if line.startswith("chromosome"):
                continue
            else:
                lineAfterSplit = line.split("\t")
                chrom = lineAfterSplit[0]
                start = lineAfterSplit[1]
                end = lineAfterSplit[2]
                gene = lineAfterSplit[3]
                log2 = lineAfterSplit[4]
                VAF = lineAfterSplit[5]
                depth = lineAfterSplit[6]
                if "p_ttest" in line:
                    p_ttest = lineAfterSplit[7]
                    probes = lineAfterSplit[8]
                    weight = lineAfterSplit[9]
                else:
                    p_ttest = "-"
                    probes = lineAfterSplit[7]
                    weight = lineAfterSplit[8]

                if gene == "-":
                    continue
                elif VAF == "2":
                    continue
                else:
                    filter_output = [chrom, start, end, gene, VAF, log2, depth, p_ttest, probes, weight]
                    outputString = "\t".join(filter_output)
                    cnvkitResults.write(outputString)
        cnvkitResults.close()
        cnvkitResultsFile.close()


    def conifer(self):
        pass

    def freec(self):
        pass

    def cnvnator(self):
        pass