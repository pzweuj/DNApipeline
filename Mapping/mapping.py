#! /usr/bin/python3
# -*- coding: utf-8 -*-

__Version__ = "0.31"
__Author__ = "pzweuj"
__Date__ = "20210416"


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
        self.pair = runningInfo["pair"]
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
    # sambamba
    # https://lomereiter.github.io/sambamba/
    def bwa_mem(self):
        reference = self.reference
        threads = self.threads
        resultsDir = self.output
        sampleID = self.sample
        pairID = self.pair

        tmpDir = resultsDir + "/tempFile/bwa_" + sampleID
        mkdir(tmpDir)
        tmp = tmpDir + "/tmp"
        mkdir(tmp)

        cmd = """
            bwa mem -t {threads} \\
                -M \\
                -R "@RG\\tID:{sampleID}\\tLB:{sampleID}\\tPL:illumina\\tPU:Hiseq\\tSM:{sampleID}" \\
                {reference} \\
                {resultsDir}/cleandata/{sampleID}.clean_R1.fastq.gz \\
                {resultsDir}/cleandata/{sampleID}.clean_R2.fastq.gz \\
                | sambamba view -f bam -t {threads} -S /dev/stdin > {tmpDir}/{sampleID}.bam
            sambamba sort {tmpDir}/{sampleID}.bam -t {threads} -o {tmpDir}/{sampleID}.sort.bam --tmpdir {tmp} -p
            rm {tmpDir}/{sampleID}.bam
            cp {tmpDir}/{sampleID}.sort.bam {resultsDir}/bam/{sampleID}.bam
            cp {tmpDir}/{sampleID}.sort.bam.bai {resultsDir}/bam/{sampleID}.bam.bai
            rm -rf {tmp}
        """.format(tmpDir=tmpDir, threads=threads, sampleID=sampleID, reference=reference, resultsDir=resultsDir, tmp=tmp)
        print(cmd)
        os.system(cmd)

        if pairID != None:
            pairDir = resultsDir + "/tempFile/bwa_" + pairID
            mkdir(pairDir)
            tmp = pairDir + "/tmp"
            mkdir(tmp)
            p = """
                bwa mem -t {threads} \\
                    -M \\
                    -R "@RG\\tID:{pairID}\\tLB:{pairID}\\tPL:illumina\\tPU:Hiseq\\tSM:{pairID}" \\
                    {reference} \\
                    {resultsDir}/cleandata/{pairID}.clean_R1.fastq.gz \\
                    {resultsDir}/cleandata/{pairID}.clean_R2.fastq.gz \\
                    | sambamba view -f bam -t {threads} -S /dev/stdin > {pairDir}/{pairID}.bam
                sambamba sort {pairDir}/{pairID}.bam -t {threads} -o {pairDir}/{pairID}.sort.bam --tmpdir {tmp} -p
                rm {pairDir}/{pairID}.bam
                cp {pairDir}/{pairID}.sort.bam {resultsDir}/bam/{pairID}.bam
                cp {pairDir}/{pairID}.sort.bam.bai {resultsDir}/bam/{pairID}.bam.bai
                rm -rf {tmp}
            """.format(pairDir=pairDir, threads=threads, pairID=pairID, reference=reference, resultsDir=resultsDir, tmp=tmp)
            print(p)
            os.system(p)

    # gencore
    # https://github.com/OpenGene/gencore
    # 当使用fastp识别UMI时使用此方法进行去重
    def gencore(self):
        resultsDir = self.output
        sampleID = self.sample
        pairID = self.pair
        threads = self.threads
        reference = self.reference

        tmpDir = resultsDir + "/tempFile/gencore_" + sampleID
        mkdir(tmpDir)
        tmp = tmpDir + "/tmp"
        mkdir(tmp)
        cmd = """
            gencore -i {resultsDir}/bam/{sampleID}.bam \\
                -r {reference} \\
                -o {tmpDir}/{sampleID}.umi.bam \\
                -u UMI -s 2 -d 1 \\
                -j {tmpDir}/{sampleID}.json -h {tmpDir}/{sampleID}.html
            sambamba sort -t {threads} {tmpDir}/{sampleID}.umi.bam -o {tmpDir}/{sampleID}.umi.sort.bam --tmpdir {tmp} -p
            cp {tmpDir}/{sampleID}.umi.sort.bam {resultsDir}/bam/{sampleID}.bam
            cp {tmpDir}/{sampleID}.umi.sort.bam.bai {resultsDir}/bam/{sampleID}.bam.bai
            rm -rf {tmp}
        """.format(reference=reference, resultsDir=resultsDir, sampleID=sampleID, tmpDir=tmpDir, threads=threads, tmp=tmp)
        print(cmd)
        os.system(cmd)

        if pairID != None:
            pairDir = resultsDir + "/tempFile/gencore_" + pairID
            mkdir(pairDir)
            p = """
                gencore -i {resultsDir}/bam/{pairID}.bam \\
                    -r {reference} \\
                    -o {pairDir}/{pairID}.umi.bam \\
                    -u UMI -s 2 -d 1 \\
                    -j {pairDir}/{pairID}.json -h {pairDir}/{pairID}.html
                sambamba sort -t {threads} {pairDir}/{pairID}.umi.bam -o {pairDir}/{pairID}.umi.sort.bam --tmpdir {tmp} -p
                cp {pairDir}/{pairID}.umi.sort.bam {resultsDir}/bam/{pairID}.bam
                cp {pairDir}/{pairID}.umi.sort.bam.bai {resultsDir}/bam/{pairID}.bam.bai
                rm -rf {tmp}
            """.format(reference=reference, resultsDir=resultsDir, pairID=pairID, pairDir=pairDir, threads=threads, tmp=tmp)
            print(p)
            os.system(p)

    # GATK4
    # https://github.com/broadinstitute/gatk
    def markDuplicates(self):
        resultsDir = self.output
        sampleID = self.sample
        pairID = self.pair
        threads = self.threads
        if self.removeDups:
            remove = "-r"
        else:
            remove = ""

        tmpDir = resultsDir + "/tempFile/markDups_" + sampleID
        mkdir(tmpDir)
        tmp = tmpDir + "/tmp"
        mkdir(tmp)
        
        cmd = """
            sambamba markdup \\
                {resultsDir}/bam/{sampleID}.bam \\
                {tmpDir}/{sampleID}.marked.bam \\
                -p --overflow-list-size 600000 \\
                --tmpdir {tmp} \\
                -t {threads} {remove}
            rm -rf {tmp}
            cp {tmpDir}/{sampleID}.marked.bam {resultsDir}/bam/{sampleID}.bam
            cp {tmpDir}/{sampleID}.marked.bam.bai {resultsDir}/bam/{sampleID}.bam.bai            
        """.format(threads=threads, tmpDir=tmpDir, resultsDir=resultsDir, sampleID=sampleID, remove=remove, tmp=tmp)
        print(cmd)
        os.system(cmd)

        if pairID != None:
            pairDir = resultsDir + "/tempFile/markDups_" + pairID
            mkdir(pairDir)
            tmp = pairDir + "/tmp"
            mkdir(tmp)
            p = """
                sambamba markdup \\
                    {resultsDir}/bam/{pairID}.bam \\
                    {pairDir}/{pairID}.marked.bam \\
                    -p --overflow-list-size 600000 \\
                    --tmpdir {tmp} \\
                    -t {threads} {remove}
                rm -rf {tmp}
                cp {pairDir}/{pairID}.marked.bam {resultsDir}/bam/{pairID}.bam
                cp {pairDir}/{pairID}.marked.bam.bai {resultsDir}/bam/{pairID}.bam.bai
            """.format(threads=threads, pairDir=pairDir, resultsDir=resultsDir, pairID=pairID, remove=remove, tmp=tmp)
            print(p)
            os.system(p)

    def recalibrator(self):
        reference = self.reference
        databases = self.databases
        snp = databases + "/" + self.runningInfo["setting"]["Mapping"]["snp_1000g"]
        mills = databases + "/" + self.runningInfo["setting"]["Mapping"]["mills"]
        indel_1000g = databases + "/" + self.runningInfo["setting"]["Mapping"]["indel_1000g"]
        resultsDir = self.output
        sampleID = self.sample
        pairID = self.pair

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

        if pairID != None:
            pairDir = resultsDir + "/tempFile/gatk_" + pairID
            mkdir(pairDir)
            p = """
                gatk BaseRecalibrator \\
                    --known-sites {snp} \\
                    --known-sites {mills} \\
                    --known-sites {indel_1000g} \\
                    -R {reference} \\
                    -I {resultsDir}/bam/{pairID}.bam \\
                    -O {pairDir}/{pairID}.recal.table

                gatk ApplyBQSR \\
                    -R {reference} \\
                    --bqsr-recal-file {pairDir}/{pairID}.recal.table \\
                    -I {resultsDir}/bam/{pairID}.bam \\
                    -O {pairDir}/{pairID}.BQSR.bam

                cp {pairDir}/{pairID}.BQSR.bam {resultsDir}/bam/{pairID}.bam
                cp {pairDir}/{pairID}.BQSR.bai {resultsDir}/bam/{pairID}.bam.bai
            """.format(pairDir=pairDir, snp=snp, mills=mills, indel_1000g=indel_1000g, reference=reference, resultsDir=resultsDir, pairID=pairID)
            print(p)
            os.system(p)


    # Gemini
    # https://github.com/Illumina/Pisces/wiki/Gemini-5.2.10-Design-Document
    # 未测试
    def gemini(self):
        gemini_multi = "/home/bioinfo/ubuntu/software/GeminiMulti_5.2.10.49/GeminiMulti.dll"

        databases = self.databases
        resultsDir = self.output
        sampleID = self.sample
        pairID = self.pair
        threads = self.threads

        tmpDir = resultsDir + "/tempFile/gemini_" + sampleID
        mkdir(tmpDir)
        cmd = """
            dotnet {gemini_multi} -bam {resultsDir}/bam/{sampleID}.bam \\
                -genome {databases} \\
                --outFolder {tmpDir} \\
                --numprocesses {threads} \\
                --samtools /home/bioinfo/ubuntu/software/samtools-1.11
        """.format(gemini_multi=gemini_multi, resultsDir=resultsDir, sampleID=sampleID, databases=databases, tmpDir=tmpDir, threads=threads)
        print(cmd)
        os.system(cmd)

        if pairID != None:
            pairDir = resultsDir + "/tempFile/gemini_" + pairID
            mkdir(pairDir)
            p = """
                dotnet {gemini_multi} -bam {resultsDir}/bam/{pairID}.bam \\
                    -genome {databases} \\
                    --outFolder {pairDir} \\
                    --numprocesses {threads} \\
                    --samtools /home/bioinfo/ubuntu/software/samtools-1.11
            """.format(gemini_multi=gemini_multi, resultsDir=resultsDir, pairID=pairID, databases=databases, pairDir=pairDir, threads=threads)
            print(p)
            os.system(p)

    # bamdst
    # https://github.com/shiquan/bamdst
    # 用于统计bam比对效果
    def bamdst(self):
        resultsDir = self.output
        sampleID = self.sample
        pairID = self.pair
        bedFile = self.bed
        
        tmpDir = resultsDir + "/tempFile/bamdst_" + sampleID
        mkdir(tmpDir)
        if pairID != None:
            pairDir = resultsDir + "/tempFile/bamdst_" + pairID
            mkdir(pairDir)

        if bedFile == None:
            print("必须导入bed文件才能进行捕获分析！")
        else:
            cmd = """
                bamdst -p {bedFile} -o {tmpDir} {resultsDir}/bam/{sampleID}.bam
            """.format(bedFile=bedFile, tmpDir=tmpDir, resultsDir=resultsDir, sampleID=sampleID)
            print(cmd)
            os.system(cmd)

            if pairID != None:
                p = """
                    bamdst -p {bedFile} -o {pairDir} {resultsDir}/bam/{pairID}.bam
                """.format(bedFile=bedFile, pairDir=pairDir, resultsDir=resultsDir, pairID=pairID)
                print(p)
                os.system(p)


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
        
        if pairID != None:
            bamdstReportFile = open(pairDir + "/coverage.report", "r")
            bamdstReport = open(resultsDir + "/QC/" + pairID + ".bamdst.txt", "w")
            bamdstReport.write(pairID + " Bamdst QC Report\n")
            for line in bamdstReportFile:
                if line.startswith("#"):
                    continue
                else:
                    lines = line.lstrip()
                    bamdstReport.write(lines)
            bamdstReport.close()
            bamdstReportFile.close()
        
        print("bamdst捕获分析完成！")

# end