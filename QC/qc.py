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


class QC(object):
    """
    质控模块
    """
    def __init__(self, runningInfo):
        self.runningInfo = runningInfo
        self.sample = runningInfo["sample"]
        self.rawdata = runningInfo["rawdata"]
        self.output = runningInfo["output"]
        
        self.threads = str(runningInfo["process"]["threads"])
        self.runApp = runningInfo["process"]["QC"]

        mkdir(self.output)
        mkdir(self.output + "/QC")
        mkdir(self.output + "/cleandata")

    # fastp
    # https://github.com/OpenGene/fastp
    def fastp(self):
        rawdataDir = self.rawdata
        sampleID = self.sample
        resultsDir = self.output
        threads = self.threads
        cmd = """
            fastp -i {rawdataDir}/{sampleID}_R1.fastq.gz \\
                -I {rawdataDir}/{sampleID}_R2.fastq.gz \\
                -o {resultsDir}/cleandata/{sampleID}.clean_R1.fastq.gz \\
                -O {resultsDir}/cleandata/{sampleID}.clean_R2.fastq.gz \\
                -j {resultsDir}/QC/{sampleID}.json \\
                -h {resultsDir}/QC/{sampleID}.html \\
                -w {threads}
        """.format(rawdataDir=rawdataDir, sampleID=sampleID, resultsDir=resultsDir, threads=threads)
        print(cmd)
        os.system(cmd)

    def cutadapt(self):
        pass

    def fastqc(self):
        pass