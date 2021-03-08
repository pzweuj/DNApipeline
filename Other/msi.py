#! /usr/bin/python3
# -*- coding: utf-8 -*-

__Version__ = "0.19"
__Author__ = "pzweuj"
__Date__ = "20210308"

import os
import sys

FUN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(FUN_DIR)
from Other.function import mkdir


class MSI(object):
    """
    MSI分析流程
    无配对样本使用msisensor-pro优先，使用若干正常样本建立基线
    """
    def __init__(self, runningInfo):
        self.runningInfo = runningInfo
        self.sample = runningInfo["sample"]
        self.rawdata = runningInfo["rawdata"]
        self.output = runningInfo["output"]

        self.threads = str(runningInfo["process"]["threads"])
        self.runApp = runningInfo["process"]["Other"]["MSI"]

        mkdir(self.output)
        mkdir(self.output + "/tempFile")
        mkdir(self.output + "/msi")

    # msisensor-pro
    # https://github.com/xjtu-omics/msisensor-pro
    def msisensor_pro(self):
        resultsDir = self.output
        sampleID = self.sample
        msi_baseline = self.runningInfo["setting"]["Other"]["msisensorpro_baseline"]
        msi_list = self.runningInfo["setting"]["Other"]["msi_list"]

        tmpDir = resultsDir + "/tempFile/msisensorpro_" + sampleID
        mkdir(tmpDir)
        cmd = """
            msisensor-pro pro -d {msi_list} \\
                -t {resultsDir}/bam/{sampleID}.bam \\
                -o {tmpDir}/{sampleID}
            mv {tmpDir}/{sampleID} {tmpDir}/{sampleID}.txt

        """.format(msi_list=msi_list, resultsDir=resultsDir, sampleID=sampleID, tmpDir=tmpDir)
        print(cmd)
        os.system(cmd)

    # msisensor2
    # https://github.com/niu-lab/msisensor2
    def msisensor2(self):
        resultsDir = self.output
        sampleID = self.sample
        threads = self.threads
        msi_model = self.runningInfo["setting"]["Other"]["MSIsensor2"]

        tmpDir = resultsDir + "/tempFile/msisensor2_" + sampleID
        mkdir(tmpDir)

        cmd = """
            msisensor2 msi -M {msi_model} -t {resultsDir}/bam/{sampleID}.bam \\
                -o {tmpDir} -b {threads}
        """.format(msi_model=msi_model, resultsDir=resultsDir, sampleID=sampleID, tmpDir=tmpDir, threads=threads)
        print(cmd)
        os.system(cmd)

    # msisensor-ct
    # https://github.com/niu-lab/msisensor-ct
    def msisensor_ct(self):
        resultsDir = self.output
        sampleID = self.sample
        threads = self.threads
        msi_model = self.runningInfo["setting"]["Other"]["MSIsensor_ct"]

        tmpDir = resultsDir + "/tempFile/msisensorct_" + sampleID
        mkdir(tmpDir)

        cmd = """
            msisensor-ct msi -D -M {msi_model} -t {resultsDir}/bam/{sampleID}.bam \\
                -o {tmpDir} -b {threads}
        """.format(msi_model=msi_model, resultsDir=resultsDir, sampleID=sampleID, tmpDir=tmpDir, threads=threads)
        print(cmd)
        os.system(cmd)        

    def visualMSI(self):
        pass

    def msing(self):
        pass

    def mantis(self):
        pass

