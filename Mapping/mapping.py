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
        self.markDups = runningInfo["setting"]["Mapping"]["markDups"]
        self.removeDups = runningInfo["setting"]["Mapping"]["removeDups"]
        self.recalibrate = runningInfo["setting"]["Mapping"]["recalibrate"]
        self.bed = runningInfo["setting"]["Mutation"]["Bed"]
        self.MAF = runningInfo["setting"]["Mutation"]["filter"]["MAF"]
        self.umi = runningInfo["setting"]["QC"]["UMI_loc"]

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

        tmpDir = resultsDir + "/tempFile/bwa_" + sampleID
        mkdir(tmpDir)
        cmd = """
            bwa mem -t {threads} \\
                -M \\
                -R "@RG\\tID:{sampleID}\\tLB:{sampleID}\\tPL:illumina\\tPU:Hiseq\\tSM:{sampleID}" \\
                {reference} \\
                {resultsDir}/cleandata/{sampleID}.clean_R1.fastq.gz \\
                {resultsDir}/cleandata/{sampleID}.clean_R2.fastq.gz \\
                | samtools view -bSh --threads {threads} - > {tmpDir}/{sampleID}.bam
            samtools sort {tmpDir}/{sampleID}.bam -@ {threads} -o {tmpDir}/{sampleID}.sort.bam
            samtools index {tmpDir}/{sampleID}.sort.bam
            cp {tmpDir}/{sampleID}.sort.bam {resultsDir}/bam/{sampleID}.bam
            cp {tmpDir}/{sampleID}.sort.bam.bai {resultsDir}/bam/{sampleID}.bam.bai
        """.format(tmpDir=tmpDir, threads=threads, sampleID=sampleID, reference=reference, resultsDir=resultsDir)
        print(cmd)
        os.system(cmd)

    def subread(self):
        pass

    def bowtie2(self):
        pass

    # gencore
    # https://github.com/OpenGene/gencore
    # 当使用fastp识别UMI时使用此方法进行去重
    def gencore(self):
        resultsDir = self.output
        sampleID = self.sample
        threads = self.threads
        reference = self.reference

        tmpDir = resultsDir + "/tempFile/gencore_" + sampleID
        mkdir(tmpDir)

        cmd = """
            gencore -i {resultsDir}/bam/{sampleID}.bam \\
                -r {reference} \\
                -o {tmpDir}/{sampleID}.umi.bam \\
                -u UMI -s 2 -d 1 \\
                -j {tmpDir}/{sampleID}.json -h {tmpDir}/{sampleID}.html
            samtools sort -@ {threads} {tmpDir}/{sampleID}.umi.bam -o {tmpDir}/{sampleID}.umi.sort.bam
            samtools index {tmpDir}/{sampleID}.umi.sort.bam
            cp {tmpDir}/{sampleID}.umi.sort.bam {resultsDir}/bam/{sampleID}.bam
            cp {tmpDir}/{sampleID}.umi.sort.bam.bai {resultsDir}/bam/{sampleID}.bam.bai
        """.format(reference=reference, resultsDir=resultsDir, sampleID=sampleID, tmpDir=tmpDir, threads=threads)
        print(cmd)
        os.system(cmd)

    # GATK4
    # https://github.com/broadinstitute/gatk
    def markDuplicates(self):
        resultsDir = self.output
        sampleID = self.sample
        if self.removeDups:
            remove = "true"
        else:
            remove = "false"

        tmpDir = resultsDir + "/tempFile/markDups_" + sampleID
        mkdir(tmpDir)
        cmd = """
            gatk MarkDuplicates \\
                -I {resultsDir}/bam/{sampleID}.bam \\
                -O {tmpDir}/{sampleID}.marked.bam \\
                -M {tmpDir}/{sampleID}.dups.txt \\
                --REMOVE_DUPLICATES {remove}
            cp {tmpDir}/{sampleID}.marked.bam {resultsDir}/bam/{sampleID}.bam
            samtools index {resultsDir}/bam/{sampleID}.bam
        """.format(tmpDir=tmpDir, resultsDir=resultsDir, sampleID=sampleID, remove=remove)
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

        tmpDir = resultsDir + "/tempFile/gatk_" + sampleID
        mkdir(tmpDir)

        cmd = """
            gatk BaseRecalibrator \\
                --known-sites {snp} \\
                --known-sites {mills} \\
                --known-sites {indel_1000g} \\
                -R {reference} \\
                -I {resultsDir}/bam/{sampleID}.bam \\
                -O {tmpDir}/{sampleID}.recal.table

            gatk ApplyBQSR \\
                -R {reference} \\
                --bqsr-recal-file {tmpDir}/{sampleID}.recal.table \\
                -I {resultsDir}/bam/{sampleID}.bam \\
                -O {tmpDir}/{sampleID}.BQSR.bam

            cp {tmpDir}/{sampleID}.BQSR.bam {resultsDir}/bam/{sampleID}.bam
            cp {tmpDir}/{sampleID}.BQSR.bai {resultsDir}/bam/{sampleID}.bam.bai
        """.format(tmpDir=tmpDir, snp=snp, mills=mills, indel_1000g=indel_1000g, reference=reference, resultsDir=resultsDir, sampleID=sampleID)
        print(cmd)
        os.system(cmd)

    # bamdst
    # https://github.com/shiquan/bamdst
    # 用于统计bam比对效果
    def bamdst(self):
        resultsDir = self.output
        sampleID = self.sample
        bedFile = self.bed
        
        tmpDir = resultsDir + "/tempFile/bamdst_" + sampleID
        mkdir(tmpDir)

        if bedFile == None:
            print("必须导入bed文件才能进行捕获分析！")
        else:
            cmd = """
                bamdst -p {bedFile} -o {tmpDir} {resultsDir}/bam/{sampleID}.bam
            """.format(bedFile=bedFile, tmpDir=tmpDir, resultsDir=resultsDir, sampleID=sampleID)
            print(cmd)
            os.system(cmd)

        bamdstReportFile = open(tmpDir + "/coverage.report", "r")
        bamdstReport = open(resultsDir + "/QC/" + sampleID + ".bamdst.txt", "w")
        bamdstReport.write(sampleID + " Bamdst QC Report\n")
        for line in bamdstReportFile:
            if line.startswith("#"):
                continue
            else:
                lines = line.lstrip()
                bamdstReport.write(lines)
        bamdstReport.close()
        bamdstReportFile.close()
        print("bamdst捕获分析完成！")

