#! /usr/bin/python3
# -*- coding: utf-8 -*-

__Version__ = "0.3"
__Author__ = "pzweuj"
__Date__ = "20210402"

import os
import sys
import json

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
        self.pair = runningInfo["pair"]
        self.rawdata = runningInfo["rawdata"]
        self.output = runningInfo["output"]
        
        self.threads = str(runningInfo["process"]["threads"])
        self.runApp = runningInfo["process"]["QC"]

        self.UMI_loc = runningInfo["setting"]["QC"]["UMI_loc"]
        self.UMI_len = runningInfo["setting"]["QC"]["UMI_len"]

        mkdir(self.output)
        mkdir(self.output + "/QC")
        mkdir(self.output + "/cleandata")
        mkdir(self.output + "/tempFile")

    # fastp
    # https://github.com/OpenGene/fastp
    def fastp(self):
        rawdataDir = self.rawdata
        sampleID = self.sample
        pairID = self.pair
        resultsDir = self.output
        threads = self.threads
        UMI_loc = self.UMI_loc
        UMI_len = self.UMI_len

        tmpDir = resultsDir + "/tempFile/fastp_" + sampleID
        mkdir(tmpDir)

        if pairID != None:
            pairDir = resultsDir + "/tempFile/fastp_" + pairID
            mkdir(pairDir)

        if UMI_loc == None:
            cmd = """
                fastp -i {rawdataDir}/{sampleID}_R1.fastq.gz \\
                    -I {rawdataDir}/{sampleID}_R2.fastq.gz \\
                    -o {resultsDir}/cleandata/{sampleID}.clean_R1.fastq.gz \\
                    -O {resultsDir}/cleandata/{sampleID}.clean_R2.fastq.gz \\
                    -j {tmpDir}/{sampleID}.json \\
                    -h {tmpDir}/{sampleID}.html \\
                    -w {threads}
            """.format(tmpDir=tmpDir, rawdataDir=rawdataDir, sampleID=sampleID, resultsDir=resultsDir, threads=threads)
        
            if pairID != None:
                p = """
                    fastp -i {rawdataDir}/{pairID}_R1.fastq.gz \\
                        -I {rawdataDir}/{pairID}_R2.fastq.gz \\
                        -o {resultsDir}/cleandata/{pairID}.clean_R1.fastq.gz \\
                        -O {resultsDir}/cleandata/{pairID}.clean_R2.fastq.gz \\
                        -j {pairDir}/{pairID}.json \\
                        -h {pairDir}/{pairID}.html \\
                        -w {threads}
                """.format(pairDir=pairDir, rawdataDir=rawdataDir, pairID=pairID, resultsDir=resultsDir, threads=threads)

        else:
            cmd = """
                fastp -i {rawdataDir}/{sampleID}_R1.fastq.gz \\
                    -I {rawdataDir}/{sampleID}_R2.fastq.gz \\
                    -o {resultsDir}/cleandata/{sampleID}.clean_R1.fastq.gz \\
                    -O {resultsDir}/cleandata/{sampleID}.clean_R2.fastq.gz \\
                    -j {tmpDir}/{sampleID}.json \\
                    -h {tmpDir}/{sampleID}.html \\
                    -w {threads} -A \\
                    -U --umi_prefix=UMI --umi_loc={UMI_loc} --umi_len={UMI_len} --umi_skip=2
            """.format(UMI_loc=UMI_loc, UMI_len=UMI_len, tmpDir=tmpDir, rawdataDir=rawdataDir, sampleID=sampleID, resultsDir=resultsDir, threads=threads)
            
            if pairID != None:
                cmd = """
                    fastp -i {rawdataDir}/{pairID}_R1.fastq.gz \\
                        -I {rawdataDir}/{pairID}_R2.fastq.gz \\
                        -o {resultsDir}/cleandata/{pairID}.clean_R1.fastq.gz \\
                        -O {resultsDir}/cleandata/{pairID}.clean_R2.fastq.gz \\
                        -j {pairDir}/{pairID}.json \\
                        -h {pairDir}/{pairID}.html \\
                        -w {threads} -A \\
                        -U --umi_prefix=UMI --umi_loc={UMI_loc} --umi_len={UMI_len} --umi_skip=2
                """.format(UMI_loc=UMI_loc, UMI_len=UMI_len, pairDir=pairDir, rawdataDir=rawdataDir, pairID=pairID, resultsDir=resultsDir, threads=threads)

        print(cmd)
        os.system(cmd)

        if pairID != None:
            print(p)
            os.system(p)

    # fastp结果整理
    def fastp_filter(self):
        rawdataDir = self.rawdata
        sampleID = self.sample
        pairID = self.pair
        resultsDir = self.output

        fastpResultsTmpDir = resultsDir + "/tempFile/fastp_" + sampleID
        fastpJsonReport = fastpResultsTmpDir + "/" + sampleID + ".json"
        jsonFile = json.load(open(fastpJsonReport, "r"))

        rawReads = jsonFile["summary"]["before_filtering"]["total_reads"]
        rawBases = jsonFile["summary"]["before_filtering"]["total_bases"]
        rawQ20Bases = jsonFile["summary"]["before_filtering"]["q20_bases"]
        rawQ30Bases = jsonFile["summary"]["before_filtering"]["q30_bases"]
        rawQ20Rate = "%.2f" % (float(rawQ20Bases / rawBases) * 100) + "%"
        rawQ30Rate = "%.2f" % (float(rawQ30Bases / rawBases) * 100) + "%"
        rawGC = "%.2f" % (jsonFile["summary"]["before_filtering"]["gc_content"] * 100) + "%"
        duplicationRate = "%.2f" % (jsonFile["duplication"]["rate"] * 100) + "%"

        cleanReads = jsonFile["summary"]["after_filtering"]["total_reads"]
        cleanBases = jsonFile["summary"]["after_filtering"]["total_bases"]
        cleanQ20Bases = jsonFile["summary"]["after_filtering"]["q20_bases"]
        cleanQ30Bases = jsonFile["summary"]["after_filtering"]["q30_bases"]
        cleanQ20Rate = "%.2f" % (float(cleanQ20Bases / cleanBases) * 100) + "%"
        cleanQ30Rate = "%.2f" % (float(cleanQ30Bases / cleanBases) * 100) + "%"
        cleanGC = "%.2f" % (jsonFile["summary"]["after_filtering"]["gc_content"] * 100) + "%"

        fastpReport = open(resultsDir + "/QC/" + sampleID + ".fastp.txt", "w")
        fastpReport.write(sampleID + " fastp QC Report\n")
        fastpReport.write("rawReads\t" + str(rawReads) + "\n")
        fastpReport.write("rawBases\t" + str(rawBases) + "\n")
        fastpReport.write("raw_q20\t" + str(rawQ20Bases) + "\n")
        fastpReport.write("raw_q20_rate\t" + rawQ20Rate + "\n")
        fastpReport.write("raw_q30\t" + str(rawQ30Bases) + "\n")
        fastpReport.write("raw_q30_rate\t" + rawQ30Rate + "\n")
        fastpReport.write("raw_GC_content\t" + rawGC + "\n")
        fastpReport.write("duplicationRate\t" + duplicationRate + "\n")
        fastpReport.write("cleanReads\t" + str(cleanReads) + "\n")
        fastpReport.write("cleanBases\t" + str(cleanBases) + "\n")
        fastpReport.write("clean_q20\t" + str(cleanQ20Bases) + "\n")
        fastpReport.write("clean_q20_rate\t" + cleanQ20Rate + "\n")
        fastpReport.write("clean_q30\t" + str(cleanQ30Bases) + "\n")
        fastpReport.write("clean_q30_rate\t" + cleanQ30Rate + "\n")
        fastpReport.write("clean_GC_content\t" + cleanGC + "\n")
        fastpReport.close()

        if pairID != None:
            fastpResultsTmpDir = resultsDir + "/tempFile/fastp_" + pairID
            fastpJsonReport = fastpResultsTmpDir + "/" + pairID + ".json"
            jsonFile = json.load(open(fastpJsonReport, "r"))

            rawReads = jsonFile["summary"]["before_filtering"]["total_reads"]
            rawBases = jsonFile["summary"]["before_filtering"]["total_bases"]
            rawQ20Bases = jsonFile["summary"]["before_filtering"]["q20_bases"]
            rawQ30Bases = jsonFile["summary"]["before_filtering"]["q30_bases"]
            rawQ20Rate = "%.2f" % (float(rawQ20Bases / rawBases) * 100) + "%"
            rawQ30Rate = "%.2f" % (float(rawQ30Bases / rawBases) * 100) + "%"
            rawGC = "%.2f" % (jsonFile["summary"]["before_filtering"]["gc_content"] * 100) + "%"
            duplicationRate = "%.2f" % (jsonFile["duplication"]["rate"] * 100) + "%"

            cleanReads = jsonFile["summary"]["after_filtering"]["total_reads"]
            cleanBases = jsonFile["summary"]["after_filtering"]["total_bases"]
            cleanQ20Bases = jsonFile["summary"]["after_filtering"]["q20_bases"]
            cleanQ30Bases = jsonFile["summary"]["after_filtering"]["q30_bases"]
            cleanQ20Rate = "%.2f" % (float(cleanQ20Bases / cleanBases) * 100) + "%"
            cleanQ30Rate = "%.2f" % (float(cleanQ30Bases / cleanBases) * 100) + "%"
            cleanGC = "%.2f" % (jsonFile["summary"]["after_filtering"]["gc_content"] * 100) + "%"

            fastpReport = open(resultsDir + "/QC/" + pairID + ".fastp.txt", "w")
            fastpReport.write(pairID + " fastp QC Report\n")
            fastpReport.write("rawReads\t" + str(rawReads) + "\n")
            fastpReport.write("rawBases\t" + str(rawBases) + "\n")
            fastpReport.write("raw_q20\t" + str(rawQ20Bases) + "\n")
            fastpReport.write("raw_q20_rate\t" + rawQ20Rate + "\n")
            fastpReport.write("raw_q30\t" + str(rawQ30Bases) + "\n")
            fastpReport.write("raw_q30_rate\t" + rawQ30Rate + "\n")
            fastpReport.write("raw_GC_content\t" + rawGC + "\n")
            fastpReport.write("duplicationRate\t" + duplicationRate + "\n")
            fastpReport.write("cleanReads\t" + str(cleanReads) + "\n")
            fastpReport.write("cleanBases\t" + str(cleanBases) + "\n")
            fastpReport.write("clean_q20\t" + str(cleanQ20Bases) + "\n")
            fastpReport.write("clean_q20_rate\t" + cleanQ20Rate + "\n")
            fastpReport.write("clean_q30\t" + str(cleanQ30Bases) + "\n")
            fastpReport.write("clean_q30_rate\t" + cleanQ30Rate + "\n")
            fastpReport.write("clean_GC_content\t" + cleanGC + "\n")
            fastpReport.close()


    # fastqc
    # https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
    def fastqc(self):
        rawdataDir = self.rawdata
        sampleID = self.sample
        pairID = self.pair
        resultsDir = self.output
        threads = self.threads

        tmpDir = resultsDir + "/tempFile/fastqc_" + sampleID
        mkdir(tmpDir)

        if pairID != None:
            pairDir = resultsDir + "/tempFile/fastqc_" + pairID
            mkdir(pairDir)

        cmd = """
            fastqc \\
                {rawdataDir}/{sampleID}_R1.fastq.gz \\
                {rawdataDir}/{sampleID}_R2.fastq.gz \\
                -t {threads} -o {tmpDir}
        
        """.format(rawdataDir=rawdataDir, sampleID=sampleID, threads=threads, tmpDir=tmpDir)
        print(cmd)
        os.system(cmd)

        if pairID != None:
            p = """
                fastqc \\
                    {rawdataDir}/{pairID}_R1.fastq.gz \\
                    {rawdataDir}/{pairID}_R2.fastq.gz \\
                    -t {threads} -o {pairDir}
            
            """.format(rawdataDir=rawdataDir, pairID=pairID, threads=threads, pairDir=pairDir)
            print(p)
            os.system(p)