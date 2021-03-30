#! /usr/bin/python3
# -*- coding: utf-8 -*-

__Version__ = "0.15"
__Author__ = "pzweuj"
__Date__ = "20210301"

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
        self.MAF = str(runningInfo["setting"]["Mutation"]["filter"]["MAF"])

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

        tmpDir = resultsDir + "/tempFile/gatk_" + sampleID
        mkdir(tmpDir)

        # 区分有无pon的情况，应可用配对样本或多个正常样本建立pon
        if pon == None:
            pon = "null"

        if bedFile == None:
            bedFile = "null"
        
        cmd = """
            gatk Mutect2 \\
                -R {reference} \\
                -I {resultsDir}/bam/{bamFile} \\
                -O {tmpDir}/{sampleID}.m2.vcf \\
                -tumor {sampleID} \\
                --germline-resource {gnomad} \\
                -pon {pon} \\
                --native-pair-hmm-threads {threads} \\
                -L {bedFile} \\
                -A Coverage -A GenotypeSummaries \\
                --genotype-germline-sites true \\
                --max-reads-per-alignment-start 0
            cp {tmpDir}/{sampleID}.m2.vcf {resultsDir}/vcf/{sampleID}.vcf
        """.format(tmpDir=tmpDir, bedFile=bedFile, pon=pon, reference=reference, resultsDir=resultsDir, sampleID=sampleID, gnomad=gnomad, threads=threads, bamFile=bamFile)
        print(cmd)
        os.system(cmd)

    # GATK4
    def gatk_haplotypecaller(self):
        reference = self.reference
        resultsDir = self.output
        sampleID = self.sample
        bedFile = self.bed
        threads = self.threads
        bedFile = self.bed

        tmpDir = resultsDir + "/tempFile/HaplotypeCaller_" + sampleID
        mkdir(tmpDir)

        if bedFile == None:
            bedFile = "null"
        cmd = """
            gatk HaplotypeCaller \\
                -R {reference} \\
                -I {resultsDir}/bam/{sampleID}.bam \\
                -O {tmpDir}/{sampleID}.htc.vcf \\
                -L {bedFile} \\
                --native-pair-hmm-threads {threads}
            cp {tmpDir}/{sampleID}.htc.vcf {resultsDir}/vcf/{sampleID}.vcf
        """.format(reference=reference, resultsDir=resultsDir, tmpDir=tmpDir, sampleID=sampleID, bedFile=bedFile, threads=threads)
        print(cmd)
        os.system(cmd)

    def mutscan(self):
        pass

    def vardict(self):
        pass

    # pisces
    # https://github.com/Illumina/Pisces
    def pisces(self):
        """
        需建立索引
        dotnet CreateGenomeSizeFile.dll \
            -g hg19/ \
            -s "Homo sapiens (UCSC hg19)" \
            -o hg19/
        """
        database = self.databases
        resultsDir = self.output
        sampleID = self.sample
        bedFile = self.bed
        threads = self.threads
        minDP = self.filtDP
        minMAF = self.MAF

        piscesBin = "/home/bioinfo/ubuntu/software/Pisces_5.2.10.49/Pisces.dll"

        tmpDir = resultsDir + "/tempFile/pisces_" + sampleID
        mkdir(tmpDir)

        if bedFile != None:
            cmd = """
                dotnet {piscesBin} -b {resultsDir}/bam/{sampleID}.bam \\
                    -g {database} \\
                    -o {tmpDir} \\
                    -t {threads} \\
                    -i {bedFile} \\
                    --mindp {minDP} \\
                    --minvf {minMAF} \\
                    --minvq 0 --threadbychr true
            """.format(bedFile=bedFile, minDP=minDP, minMAF=minMAF, piscesBin=piscesBin, resultsDir=resultsDir, sampleID=sampleID, database=database, tmpDir=tmpDir, threads=threads)
        else:
            cmd = """
                dotnet {piscesBin} -b {resultsDir}/bam/{sampleID}.bam \\
                    -g {database} \\
                    -o {tmpDir} \\
                    -t {threads} \\
                    --mindp {minDP} \\
                    --minvf {minMAF} \\
                    --minvq 0 --threadbychr true
            """.format(minDP=minDP, minMAF=minMAF, piscesBin=piscesBin, resultsDir=resultsDir, sampleID=sampleID, database=database, tmpDir=tmpDir, threads=threads)
        print(cmd)
        os.system(cmd)

        filt = """
            bcftools view \\
                -e "GT='0/0' | GT='./.' | GT='0/.'" \\
                {tmpDir}/{sampleID}.genome.vcf > {tmpDir}/{sampleID}.muts.vcf
            bcftools view \\
                -e "FILTER='LowDP'" \\
                {tmpDir}/{sampleID}.muts.vcf > {tmpDir}/{sampleID}.pisces.vcf
            cp {tmpDir}/{sampleID}.pisces.vcf {resultsDir}/vcf/{sampleID}.vcf
        """.format(tmpDir=tmpDir, sampleID=sampleID, resultsDir=resultsDir)
        print(filt)
        os.system(filt)

    # bcftools
    # http://samtools.github.io/bcftools/bcftools.html
    def bcftools(self):
        pass

    # freebayes
    # https://github.com/freebayes/freebayes
    # 为了避免后续annovar注释的bug，此处设定genotyping-max-banddepth
    # 同一位点只输出最多突变条数的突变方向
    def freebayes(self):
        reference = self.reference
        resultsDir = self.output
        sampleID = self.sample
        bedFile = self.bed
        MAF = self.MAF
        DP = self.filtDP

        checkBQSR = sampleID + ".BQSR.bam"
        if checkBQSR in os.listdir(resultsDir + "/bam"):
            bamFile = checkBQSR
        else:
            bamFile = sampleID + ".bam"

        tmpDir = resultsDir + "/tempFile/freebayes_" + sampleID
        mkdir(tmpDir)

        if bedFile != None:
            cmd = """
                freebayes -f {reference} \\
                    {resultsDir}/bam/{bamFile} \\
                    -t {bedFile} \\
                    -F {MAF} -C 5 --min-coverage {DP} \\
                    --genotyping-max-banddepth 1 \\
                    > {tmpDir}/{sampleID}.freebayes.vcf
                sed -i "s/0\\/0/0\\/1/g" {tmpDir}/{sampleID}.freebayes.vcf
                cp {tmpDir}/{sampleID}.freebayes.vcf {resultsDir}/vcf/{sampleID}.vcf
            """.format(DP=DP, MAF=MAF, tmpDir=tmpDir, bedFile=bedFile, bamFile=bamFile, reference=reference, resultsDir=resultsDir, sampleID=sampleID)
        else:
            cmd = """
                freebayes -f {reference} \\
                    {resultsDir}/bam/{bamFile} \\
                    -F {MAF} -C 5 --min-coverage {DP} \\
                    --genotyping-max-banddepth 1 \\
                    > {tmpDir}/{sampleID}.freebayes.vcf
                sed -i "s/0\\/0/0\\/1/g" {tmpDir}/{sampleID}.freebayes.vcf
                cp {tmpDir}/{sampleID}.freebayes.vcf {resultsDir}/vcf/{sampleID}.vcf
            """.format(DP=DP, MAF=MAF, tmpDir=tmpDir, bedFile=bedFile, bamFile=bamFile, reference=reference, resultsDir=resultsDir, sampleID=sampleID)

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

        tmpDir = resultsDir + "/tempFile/gatk_" + sampleID
        mkdir(tmpDir)
        
        cmd = """
            gatk GetPileupSummaries \\
                -I {resultsDir}/bam/{bamFile} \\
                -O {tmpDir}/{sampleID}.pileups.table \\
                -V {small_exac} \\
                -L {bedFile} \\
                -R {reference}

            gatk CalculateContamination \\
                -I {tmpDir}/{sampleID}.pileups.table \\
                -O {tmpDir}/{sampleID}.contamination.table

            gatk FilterMutectCalls \\
                -R {reference} \\
                -V {tmpDir}/{sampleID}.m2.vcf \\
                -O {tmpDir}/{sampleID}.m2.contFiltered.vcf \\
                --contamination-table {tmpDir}/{sampleID}.contamination.table

            bcftools view \\
                {tmpDir}/{sampleID}.m2.contFiltered.vcf \\
                -f PASS,clustered_events,slippage \\
                > {tmpDir}/{sampleID}.filter.vcf
            cp {tmpDir}/{sampleID}.filter.vcf {resultsDir}/vcf/{sampleID}.vcf
        """.format(tmpDir=tmpDir, bamFile=bamFile, resultsDir=resultsDir, sampleID=sampleID, small_exac=small_exac, bedFile=bedFile, reference=reference)
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