#! /usr/bin/python3
# -*- coding: utf-8 -*-

__Version__ = "0.1"
__Author__ = "pzweuj"
__Date__ = "20210220"


import os
import sys

FUN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(FUN_DIR)
from Other.function import mkdir

class Mapping(object):
    """
    比对模块
    由于比对模块的目的是生成bam文件
    因此把同样是生成bam文件的标记重复以及校对均归类到比对模块中
    """
    def __init__(self, runningInfo):
        self.runningInfo = runningInfo
        self.sample = runningInfo["sample"]
        self.rawdata = runningInfo["rawdata"]
        self.output = runningInfo["output"]
        
        self.threads = str(runningInfo["process"]["threads"])
        self.runApp = runningInfo["process"]["Mapping"]
        
        self.databases = runningInfo["setting"]["Mapping"]["databases"]
        self.reference = runningInfo["setting"]["Mapping"]["reference"]
        self.removeDups = runningInfo["setting"]["Mapping"]["removeDups"]
        self.recalibrate = runningInfo["setting"]["Mapping"]["recalibrate"]

        mkdir(self.output)
        mkdir(self.output + "/tempFile")
        mkdir(self.output + "/bam")

    # bwa
    # https://github.com/lh3/bwa
    # samtools
    # https://github.com/samtools/samtools
    def bwa_mem(self):
        reference = self.reference
        threads = self.threads
        resultsDir = self.output
        sampleID = self.sample
        cmd = """
            bwa mem -t {threads} \\
                -M \\
                -R "@RG\\tID:{sampleID}\\tLB:{sampleID}\\tPL:illumina\\tPU:Hiseq\\tSM:{sampleID}" \\
                {reference} \\
                {resultsDir}/cleandata/{sampleID}.clean_R1.fastq.gz \\
                {resultsDir}/cleandata/{sampleID}.clean_R2.fastq.gz \\
                | samtools view -bSh --threads {threads} - > {resultsDir}/tempFile/{sampleID}.bam
            samtools sort {resultsDir}/tempFile/{sampleID}.bam -@ {threads} -o {resultsDir}/bam/{sampleID}.bam
            samtools index {resultsDir}/bam/{sampleID}.bam
        """.format(threads=threads, sampleID=sampleID, reference=reference, resultsDir=resultsDir)
        print(cmd)
        os.system(cmd)

    def subread(self):
        pass

    def bowtie2(self):
        pass

    # GATK4
    # https://github.com/broadinstitute/gatk
    def markDuplicates(self):
        resultsDir = self.output
        sampleID = self.sample
        cmd = """
            gatk MarkDuplicates \\
                -I {resultsDir}/bam/{sampleID}.bam \\
                -O {resultsDir}/bam/{sampleID}.marked.bam \\
                -M {resultsDir}/tempFile/{sampleID}.dups.txt \\
                --REMOVE_DUPLICATES true
            rm {resultsDir}/bam/{sampleID}.bam*
            mv {resultsDir}/bam/{sampleID}.marked.bam {resultsDir}/bam/{sampleID}.bam
            samtools index {resultsDir}/bam/{sampleID}.bam
        """.format(resultsDir=resultsDir, sampleID=sampleID)
        print(cmd)
        os.system(cmd)

    def recalibrator(self):
        reference = self.reference
        databases = self.databases
        snp = databases + "/" + self.runningInfo["setting"]["Mapping"]["snp_1000g"]
        mills = databases + "/" + self.runningInfo["setting"]["Mapping"]["mills"]
        indel_1000g = databases + "/" + self.runningInfo["setting"]["Mapping"]["indel_1000g"]
        resultsDir = self.output
        sampleID = self.sample
        cmd = """
            gatk BaseRecalibrator \\
                --known-sites {snp} \\
                --known-sites {mills} \\
                --known-sites {indel_1000g} \\
                -R {reference} \\
                -I {resultsDir}/bam/{sampleID}.bam \\
                -O {resultsDir}/tempFile/{sampleID}.recal.table

            gatk ApplyBQSR \\
                -R {reference} \\
                --bqsr-recal-file {resultsDir}/tempFile/{sampleID}.recal.table \\
                -I {resultsDir}/bam/{sampleID}.bam \\
                -O {resultsDir}/bam/{sampleID}.BQSR.bam

            mv {resultsDir}/bam/{sampleID}.BQSR.bai {resultsDir}/bam/{sampleID}.BQSR.bam.bai
        """.format(snp=snp, mills=mills, indel_1000g=indel_1000g, reference=reference, resultsDir=resultsDir, sampleID=sampleID)
        print(cmd)
        os.system(cmd)        