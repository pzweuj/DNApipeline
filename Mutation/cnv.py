#! /usr/bin/python3
# -*- coding: utf-8 -*-

__Version__ = "0.14"
__Author__ = "pzweuj"
__Date__ = "20210407"

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
        self.pair = runningInfo["pair"]
        self.rawdata = runningInfo["rawdata"]
        self.output = runningInfo["output"]

        self.reference = runningInfo["setting"]["Mapping"]["reference"]
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
        pairID = self.pair
        CNV = self.runningInfo["setting"]["Mutation"]["CNV"]
        baseline = CNV["baseline"]
        reference = self.reference
        target = CNV["target"]
        antitarget = CNV["antitarget"]
        refFlat = self.runningInfo["setting"]["Annotation"]["refFlat"]
        threads = self.threads

        mkdir(resultsDir + "/cnv")
        tmpDir = resultsDir + "/tempFile/cnvkit_" + sampleID
        mkdir(tmpDir)

        if pairID == None:
            cmd = """
                cnvkit.py coverage {resultsDir}/bam/{sampleID}.bam \\
                    {target} -o {tmpDir}/{sampleID}.targetcoverage.cnn -p {threads}
                cnvkit.py coverage {resultsDir}/bam/{sampleID}.bam \\
                    {antitarget} -o {tmpDir}/{sampleID}.antitargetcoverage.cnn -p {threads}
                cnvkit.py fix {tmpDir}/{sampleID}.targetcoverage.cnn \\
                    {tmpDir}/{sampleID}.antitargetcoverage.cnn \\
                    {baseline} -o {tmpDir}/{sampleID}.cnr
                cnvkit.py segment {tmpDir}/{sampleID}.cnr \\
                    -o {tmpDir}/{sampleID}.cns -p {threads}
                cnvkit.py call {tmpDir}/{sampleID}.cns \\
                    -o {tmpDir}/{sampleID}.call.cns
                cnvkit.py scatter {tmpDir}/{sampleID}.cnr \\
                    -s {tmpDir}/{sampleID}.cns -o {tmpDir}/{sampleID}.scatter.pdf
                cnvkit.py diagram {tmpDir}/{sampleID}.cnr \\
                    -s {tmpDir}/{sampleID}.cns -o {tmpDir}/{sampleID}.diagram.pdf
            """.format(resultsDir=resultsDir, sampleID=sampleID, target=target, tmpDir=tmpDir, antitarget=antitarget, baseline=baseline, threads=threads)
        else:
            cmd = """
                cnvkit.py coverage {resultsDir}/bam/{pairID}.bam \\
                    {target} -o {tmpDir}/{pairID}.targetcoverage.cnn -p {threads}
                cnvkit.py coverage {resultsDir}/bam/{sampleID}.bam \\
                    {target} -o {tmpDir}/{sampleID}.targetcoverage.cnn -p {threads}
                cnvkit.py coverage {resultsDir}/bam/{pairID}.bam \\
                    {antitarget} -o {tmpDir}/{pairID}.antitargetcoverage.cnn -p {threads}
                cnvkit.py coverage {resultsDir}/bam/{sampleID}.bam \\
                    {antitarget} -o {tmpDir}/{sampleID}.antitargetcoverage.cnn -p {threads}

                cnvkit.py reference {tmpDir}/{pairID}.antitargetcoverage.cnn \\
                    {tmpDir}/{pairID}.targetcoverage.cnn \\
                    --fasta {reference} -o {tmpDir}/{pairID}.cnn

                cnvkit.py fix {tmpDir}/{sampleID}.targetcoverage.cnn \\
                    {tmpDir}/{sampleID}.antitargetcoverage.cnn \\
                    {tmpDir}/{pairID}.cnn -o {tmpDir}/{sampleID}.cnr
                cnvkit.py segment {tmpDir}/{sampleID}.cnr \\
                    -o {tmpDir}/{sampleID}.cns -p {threads}
                cnvkit.py call {tmpDir}/{sampleID}.cns \\
                    -o {tmpDir}/{sampleID}.call.cns
                cnvkit.py scatter {tmpDir}/{sampleID}.cnr \\
                    -s {tmpDir}/{sampleID}.cns -o {tmpDir}/{sampleID}.scatter.pdf
                cnvkit.py diagram {tmpDir}/{sampleID}.cnr \\
                    -s {tmpDir}/{sampleID}.cns -o {tmpDir}/{sampleID}.diagram.pdf
            """.format(resultsDir=resultsDir, sampleID=sampleID, pairID=pairID, target=target, antitarget=antitarget, tmpDir=tmpDir, reference=reference, threads=threads)

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
                lineAfterSplit = line.split("\n")[0].split("\t")
                chrom = lineAfterSplit[0]
                start = lineAfterSplit[1]
                end = lineAfterSplit[2]
                gene = lineAfterSplit[3]
                log2 = lineAfterSplit[4]
                VAF = lineAfterSplit[5]
                depth = lineAfterSplit[6]
                if len(lineAfterSplit) == 10:
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
                    cnvkitResults.write(outputString + "\n")
        cnvkitResults.close()
        cnvkitResultsFile.close()

    # conifer
    # http://conifer.sourceforge.net/
    def conifer(self):
        pass

    # freec
    # http://boevalab.inf.ethz.ch/FREEC/
    def freec(self):
        pass

    # cnvnator
    # https://github.com/abyzovlab/CNVnator
    def cnvnator(self):
        pass