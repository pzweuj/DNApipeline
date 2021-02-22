#! /usr/bin/python3
# -*- coding: utf-8 -*-

__Version__ = "0.1"
__Author__ = "pzweuj"
__Date__ = "20210222"

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

        cmd = """
            cnvkit.py batch \\
                {resultsDir}/bam/{sampleID}.bam \\
                -r {referenceCnn} \\
                -d {resultsDir}/tempFile/{sampleID}_cnvkit \\
                -m hybrid \\
                --diagram --scatter
            cp {resultsDir}/tempFile/{sampleID}_cnvkit/{sampleID}.cnr {resultsDir}/cnv/{sampleID}.cnr
        """.format(resultsDir=resultsDir, sampleID=sampleID, referenceCnn=referenceCnn)
        print(cmd)
        os.system(cmd)



    def conifer(self):
        pass

    def freec(self):
        pass

    def cnvnator(self):
        pass