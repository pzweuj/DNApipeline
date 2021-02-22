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

class SNV_Indel(object):
    """
    SNV与indel检测模块，最终生成vcf文件
    为了流程的一致，每个caller最终生成的结果在vcf文件夹中命名为sample.vcf
    """
    def __init__(self, runningInfo):
        self.runningInfo = runningInfo
        self.sample = runningInfo["sample"]
        self.rawdata = runningInfo["rawdata"]
        self.output = runningInfo["output"]
        
        self.threads = str(runningInfo["process"]["threads"])
        self.runApp = runningInfo["process"]["Mutation"]["SNV_indel"]

        self.databases = runningInfo["setting"]["Mapping"]["databases"]
        self.reference = runningInfo["setting"]["Mapping"]["reference"]
        self.bed = runningInfo["setting"]["Mutation"]["Bed"]
        self.filtDP = str(runningInfo["setting"]["Mutation"]["filter"]["DP"])
        self.filtQUAL = str(runningInfo["setting"]["Mutation"]["filter"]["QUAL"])

        mkdir(self.output)
        mkdir(self.output + "/tempFile")
        mkdir(self.output + "/vcf")

    # GATK4
    # https://github.com/broadinstitute/gatk
    def gatk_m2(self):
        reference = self.reference
        gnomad = self.databases + "/" + self.runningInfo["setting"]["Mutation"]["gnomad"]
        resultsDir = self.output
        sampleID = self.sample
        pon = self.runningInfo["setting"]["Mutation"]["pon"]
        bedFile = self.bed
        threads = self.threads

        checkBQSR = sampleID + ".BQSR.bam"
        if checkBQSR in os.listdir(resultsDir + "/bam"):
            bamFile = checkBQSR
        else:
            bamFile = sampleID + ".bam"

        # 区分有无pon的情况，应可用配对样本或多个正常样本建立pon
        if pon == None:
            cmd = """
                gatk Mutect2 \\
                    -R {reference} \\
                    -I {resultsDir}/bam/{bamFile} \\
                    -O {resultsDir}/tempFile/{sampleID}.m2.vcf \\
                    -tumor {sampleID} \\
                    --af-of-alleles-not-in-resource 0.0000025 \\
                    --germline-resource {gnomad} \\
                    --native-pair-hmm-threads {threads} \\
                    -L {bedFile}
                cp {resultsDir}/tempFile/{sampleID}.m2.vcf {resultsDir}/vcf/{sampleID}.vcf
            """.format(bedFile=bedFile, reference=reference, resultsDir=resultsDir, sampleID=sampleID, gnomad=gnomad, threads=threads, bamFile=bamFile)
        else:
            cmd = """
                gatk Mutect2 \\
                    -R {reference} \\
                    -I {resultsDir}/bam/{bamFile} \\
                    -O {resultsDir}/tempFile/{sampleID}.m2.vcf \\
                    -tumor {sampleID} \\
                    --af-of-alleles-not-in-resource 0.0000025 \\
                    --germline-resource {gnomad} \\
                    -pon {pon} \\
                    --native-pair-hmm-threads {threads} \\
                    -L {bedFile}
                cp {resultsDir}/tempFile/{sampleID}.m2.vcf {resultsDir}/vcf/{sampleID}.vcf
            """.format(bedFile=bedFile, pon=pon, reference=reference, resultsDir=resultsDir, sampleID=sampleID, gnomad=gnomad, threads=threads, bamFile=bamFile)            
        print(cmd)
        os.system(cmd)

    def gatk_haplotypecaller(self):
        pass

    def mutscan(self):
        pass

    # bcftools
    # http://samtools.github.io/bcftools/bcftools.html
    def bcftools(self):
        pass

    # freebayes
    # https://github.com/freebayes/freebayes
    def freebayes(self):
        reference = self.reference
        resultsDir = self.output
        sampleID = self.sample
        bedFile = self.bed

        checkBQSR = sampleID + ".BQSR.bam"
        if checkBQSR in os.listdir(resultsDir + "/bam"):
            bamFile = checkBQSR
        else:
            bamFile = sampleID + ".bam"

        cmd = """
            freebayes -f {reference} \\
                {resultsDir}/bam/{bamFile} \\
                -t {bedFile} \\
                > {resultsDir}/tempFile/{sampleID}.freebayes.vcf
            cp {resultsDir}/tempFile/{sampleID}.freebayes.vcf {resultsDir}/vcf/{sampleID}.vcf
        """.format(bedFile=bedFile, bamFile=bamFile, reference=reference, resultsDir=resultsDir, sampleID=sampleID)
        print(cmd)
        os.system(cmd)

    # GATK4 mutect2 过滤
    def gatk_filter(self):
        small_exac = self.databases + "/" + self.runningInfo["setting"]["Mutation"]["gatk_filter"]["small_exac"]
        bedFile = self.bed
        reference = self.reference
        resultsDir = self.output
        sampleID = self.sample

        checkBQSR = sampleID + ".BQSR.bam"
        if checkBQSR in os.listdir(resultsDir + "/bam"):
            bamFile = checkBQSR
        else:
            bamFile = sampleID + ".bam"

        cmd = """
            gatk GetPileupSummaries \\
                -I {resultsDir}/bam/{bamFile} \\
                -O {resultsDir}/tempFile/{sampleID}.pileups.table \\
                -V {small_exac} \\
                -L {bedFile} \\
                -R {reference}

            gatk CalculateContamination \\
                -I {resultsDir}/tempFile/{sampleID}.pileups.table \\
                -O {resultsDir}/tempFile/{sampleID}.contamination.table

            gatk FilterMutectCalls \\
                -R {reference} \\
                -V {resultsDir}/tempFile/{sampleID}.m2.vcf \\
                -O {resultsDir}/tempFile/{sampleID}.m2.contFiltered.vcf \\
                --contamination-table {resultsDir}/tempFile/{sampleID}.contamination.table

            bcftools view \\
                {resultsDir}/tempFile/{sampleID}.m2.contFiltered.vcf \\
                -f PASS,clustered_events,slippage \\
                > {resultsDir}/vcf/{sampleID}.vcf
        """.format(bamFile=bamFile, resultsDir=resultsDir, sampleID=sampleID, small_exac=small_exac, bedFile=bedFile, reference=reference)
        print(cmd)
        os.system(cmd)

    # 过滤
    # bcftools
    # http://samtools.github.io/bcftools/bcftools.html
    def filter(self):
        resultsDir = self.output
        sampleID = self.sample
        DP = self.filtDP
        QUAL = self.filtQUAL

        cmd = """
            bcftools view \\
                -i 'MIN(FORMAT/DP)>={DP}' \\
                {resultsDir}/vcf/{sampleID}.vcf \\
                > {resultsDir}/vcf/{sampleID}.filter.vcf
        """.format(QUAL=QUAL, DP=DP, resultsDir=resultsDir, sampleID=sampleID)
        print(cmd)
        os.system(cmd)