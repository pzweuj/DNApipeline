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

class SV(object):
    """
    结构变异与融合检测模块
    """
    def __init__(self, runningInfo):
        self.runningInfo = runningInfo
        self.sample = runningInfo["sample"]
        self.rawdata = runningInfo["rawdata"]
        self.output = runningInfo["output"]

        self.threads = str(runningInfo["process"]["threads"])
        self.runApp = runningInfo["process"]["Mutation"]["SV"]

        mkdir(self.output)
        mkdir(self.output + "/tempFile")
        mkdir(self.output + "/sv")

    # lumpy
    # https://github.com/arq5x/lumpy-sv
    # svtyper
    # https://github.com/hall-lab/svtyper
    def lumpy(self):
        resultsDir = self.output
        sampleID = self.sample
        threads = self.threads

        cmd = """
            samtools view -bh -F 1294 {resultsDir}/bam/{sampleID}.bam \\
                | samtools sort -@ {threads} - \\
                -o {resultsDir}/tempFile/{sampleID}.discordants.bam
            samtools index {resultsDir}/tempFile/{sampleID}.discordants.bam
            samtools view -h {resultsDir}/bam/{sampleID}.bam \\
                | extractSplitReads_BwaMem \\
                -i stdin \\
                | samtools view -bSh - \\
                | samtools sort -@ {threads} - \\
                -o {resultsDir}/tempFile/{sampleID}.splitters.bam
            samtools index {resultsDir}/tempFile/{sampleID}.splitters.bam
            lumpyexpress -B {resultsDir}/bam/{sampleID}.bam \\
                -D {resultsDir}/tempFile/{sampleID}.discordants.bam \\
                -S {resultsDir}/tempFile/{sampleID}.splitters.bam \\
                -O {resultsDir}/tempFile/{sampleID}.lumpy.vcf
            svtyper-sso \\
                -i {resultsDir}/tempFile/{sampleID}.lumpy.vcf \\
                -B {resultsDir}/bam/{sampleID}.bam \\
                --cores {threads} \\
                -o {resultsDir}/sv/{sampleID}.vcf
        """.format(resultsDir=resultsDir, sampleID=sampleID, threads=threads)
        print(cmd)
        os.system(cmd)

    def star_fusion(self):
        pass

    def subread_junction(self):
        pass

    def breakdancer(self):
        pass

    def crest(self):
        pass