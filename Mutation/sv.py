#! /usr/bin/python3
# -*- coding: utf-8 -*-

__Version__ = "0.16"
__Author__ = "pzweuj"
__Date__ = "20210303"


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
    # 当前建议使用manta，因为manta可同时输出read depth
    def manta(self):
        manta = "/mnt/d/ubuntu/software/manta-1.6.0.centos6_x86_64/bin/configManta.py"
        reference = self.reference
        tumorBam = self.output + "/bam/" + self.sample + ".bam"
        resultsDir = self.output
        sampleID = self.sample

        tmpDir = resultsDir + "/tempFile/manta_" + sampleID
        mkdir(tmpDir) 

        cmd = """
            rm -rf {tmpDir}/*
            {manta} \\
                --tumorBam {tumorBam} \\
                --referenceFasta {reference} \\
                --exome \\
                --generateEvidenceBam \\
                --runDir {tmpDir}
            {tmpDir}/runWorkflow.py
            zcat {tmpDir}/results/variants/tumorSV.vcf.gz > {tmpDir}/{sampleID}.manta.vcf
        """.format(sampleID=sampleID, manta=manta, tumorBam=tumorBam, reference=reference, tmpDir=tmpDir)
        print(cmd)
        os.system(cmd)

    # manta结果过滤
    def manta_filter(self):
        resultsDir = self.output
        sampleID = self.sample

        tmpDir = resultsDir + "/tempFile/manta_" + sampleID
        mantaResultsFile = tmpDir + "/" + sampleID + ".manta.vcf"

        filterDP = self.runningInfo["setting"]["Mutation"]["filter"]["DP"]
        filterVAF = self.runningInfo["setting"]["Mutation"]["filter"]["MAF"]
        filterAlt = int(filterDP * filterVAF)

        mantaResults = open(mantaResultsFile, "r")
        mantaResultsFilter = open(mantaResultsFile.replace(".vcf", ".filter.vcf"), "w")
        for line in mantaResults:
            if line.startswith("#"):
                mantaResultsFilter.write(line)
            else:
                lineAfterSplit = line.split("\n")[0].split("\t")
                FORMAT = lineAfterSplit[8]
                FORMAT_info = lineAfterSplit[9]

                if "PR" not in FORMAT:
                    SR = FORMAT_info.split(",")
                    PR_ref = "0"
                    PR_alt = "0"
                    SR_ref = SR[0]
                    SR_alt = SR[1]
                elif "SR" not in FORMAT:
                    PR = FORMAT_info.split(",")
                    SR_ref = "0"
                    SR_alt = "0"
                    PR_ref = PR[0]
                    PR_alt = PR[1]
                else:
                    FORMAT_infos = FORMAT_info.split(":")
                    PR = FORMAT_infos[0].split(",")
                    SR = FORMAT_infos[1].split(",")
                    PR_ref = PR[0]
                    PR_alt = PR[1]
                    SR_ref = SR[0]
                    SR_alt = SR[1]
                
                Ref_total = int(PR_ref) + int(SR_ref)
                Alt_total = int(PR_alt) + int(SR_alt)
                if Ref_total < filterDP:
                    continue
                if Alt_total < filterAlt:
                    continue

                mantaResultsFilter.write(line)
        mantaResultsFilter.close()
        mantaResults.close()
        os.system("cp " + tmpDir + "/" + sampleID + ".manta.filter.vcf " + resultsDir + "/sv/" + sampleID + ".sv.vcf")

    # factera
    # https://factera.stanford.edu/
    def factera(self):
        factera = "/home/bioinfo/ubuntu/software/factera/factera.pl"
        resultsDir = self.output
        sampleID = self.sample
        threads = self.threads

        tmpDir = resultsDir + "/tempFile/factera_" + sampleID
        mkdir(tmpDir)

        cmd = """
            {factera} -F -p {threads} \\
                -o {tmpDir} \\
                {resultsDir}/bam/{sampleID}.bam \\
                {exonBed} \\
                {referenceTwoBit}
        """.format(factera=factera, resultsDir=resultsDir, sampleID=sampleID, exonBed=exonBed, referenceTwoBit=referenceTwoBit, tmpDir=tmpDir, threads=threads)
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

    def svict(self):
        pass