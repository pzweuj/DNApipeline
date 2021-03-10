#! /usr/bin/python3
# -*- coding: utf-8 -*-

__Version__ = "0.12"
__Author__ = "pzweuj"
__Date__ = "20210309"

import os
import sys

FUN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(FUN_DIR)
from Other.function import mkdir


class HLA(object):
    """
    HLA分型分析模块
    HLA 区域为
    https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37
    chr6:28477797-33448354(hg19)
    
    https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38
    chr6:28510120-33480577(hg38)
    """
    def __init__(self, runningInfo):
        self.runningInfo = runningInfo
        self.sample = runningInfo["sample"]
        self.rawdata = runningInfo["rawdata"]
        self.output = runningInfo["output"]

        self.threads = str(runningInfo["process"]["threads"])
        self.runApp = runningInfo["process"]["Other"]["HLA"]
        self.buildver = runningInfo["setting"]["Annotation"]["buildver"]

        mkdir(self.output)
        mkdir(self.output + "/tempFile")
        mkdir(self.output + "/HLA")

    # HLA-HD
    # https://www.genome.med.kyoto-u.ac.jp/HLA-HD/
    def hlahd(self):
        resultsDir = self.output
        sampleID = self.output
        buildver = self.buildver
        threads = self.threads

        # 以下数据库无需指定参考基因坐标，为通用数据库，因此不写入配置文件中
        freq = "/home/bioinfo/ubuntu/software/hlahd.1.3.0/freq_data"
        dictionary = "/home/bioinfo/ubuntu/software/hlahd.1.3.0/dictionary"
        split_file = "/home/bioinfo/ubuntu/software/hlahd.1.3.0/HLA_gene.split.3.32.0.txt"

        if buildver == "hg19":
            extractRegion = "chr6:28477797-33448354"
        elif buildver == "b37":
            extractRegion = "6:28477797-33448354"
        elif buildver == "hg38":
            extractRegion = "chr6:28510120-33480577"
        else:
            print("Cannot get buildver")
            exit()

        tmpDir = resultsDir + "/tempFile/hlahd_" + sampleID
        mkdir(tmpDir)

        cmd = """
            samtools view {resultsDir}/bam/{sampleID}.bam {extractRegion} -b > {tmpDir}/{sampleID}.HLA.bam
            bedtools bamtofastq -i {tmpDir}/{sampleID}.HLA.bam \\
                -fq {tmpDir}/{sampleID}.HLA.R1.fastq -fq2 {tmpDir}/{sampleID}.HLA.R2.fastq
            hlahd.sh -t {threads} -m 100 -c 0.95 -f {freq} \\
                {tmpDir}/{sampleID}.HLA.R1.fastq {tmpDir}/{sampleID}.HLA.R2.fastq \\
                {split_file} {dictionary} {sampleID} {tmpDir}
            cp {tmpDir}/{sampleID}/result/G3761_final.result.txt {resultsDir}/HLA/
        """.format(threads=threads, freq=freq, split_file=split_file, dictionary=dictionary, resultsDir=resultsDir, sampleID=sampleID, extractRegion=extractRegion, tmpDir=tmpDir)
        print(cmd)
        os.system(cmd)


    def seq2hla(self):
        pass
    