#! /usr/bin/python3
# -*- coding: utf-8 -*-

__Version__ = "0.1"
__Author__ = "pzweuj"
__Date__ = "20210223"

import os
import sys

FUN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(FUN_DIR)
from Other.function import mkdir


class MSI(object):
    """
    MSI分析流程
    无配对样本使用msisensor2优先
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

    def msisensor(self):
        pass

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

    def visualMSI(self):
        pass

    def msing(self):
        pass

    def mantis(self):
        pass

