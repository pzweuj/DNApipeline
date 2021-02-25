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
        self.reference = runningInfo["setting"]["Mapping"]["reference"]
        self.bed = runningInfo["setting"]["Mutation"]["Bed"]

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

        tmpDir = resultsDir + "/tempFile/lumpy_" + sampleID
        mkdir(tmpDir)

        cmd = """
            samtools view -bh -F 1294 {resultsDir}/bam/{sampleID}.bam \\
                | samtools sort -@ {threads} - \\
                -o {tmpDir}/{sampleID}.discordants.bam
            samtools index {tmpDir}/{sampleID}.discordants.bam
            samtools view -h {resultsDir}/bam/{sampleID}.bam \\
                | extractSplitReads_BwaMem \\
                -i stdin \\
                | samtools view -bSh - \\
                | samtools sort -@ {threads} - \\
                -o {tmpDir}/{sampleID}.splitters.bam
            samtools index {tmpDir}/{sampleID}.splitters.bam
            lumpyexpress -B {resultsDir}/bam/{sampleID}.bam \\
                -D {tmpDir}/{sampleID}.discordants.bam \\
                -S {tmpDir}/{sampleID}.splitters.bam \\
                -o {tmpDir}/{sampleID}.lumpy.vcf
            svtyper-sso \\
                -i {tmpDir}/{sampleID}.lumpy.vcf \\
                -B {resultsDir}/bam/{sampleID}.bam \\
                --cores {threads} \\
                -o {tmpDir}/{sampleID}.gt.vcf
            cp {tmpDir}/{sampleID}.gt.vcf {resultsDir}/sv/{sampleID}.vcf
        """.format(tmpDir=tmpDir, resultsDir=resultsDir, sampleID=sampleID, threads=threads)
        print(cmd)
        os.system(cmd)

    # manta
    # https://github.com/Illumina/manta
    def manta(self):
        manta = "/mnt/d/ubuntu/software/manta-1.6.0.centos6_x86_64/bin/configManta.py"
        reference = self.reference
        tumorBam = self.output + "/bam/" + self.sample + ".bam"
        bedFile = self.bed

        tmpDir = self.output + "/tempFile/manta_" + self.sample
        mkdir(tmpDir) 

        cmd = """
            {manta} \\
                --tumorBam {tumorBam} \\
                --referenceFasta {reference} \\
                --exome \\
                --callRegions {bedFile} \\
                --runDir {tmpDir}
        """.format(manta=manta, tumorBam=tumorBam, reference=reference, bedFile=bedFile, tmpDir=tmpDir)
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