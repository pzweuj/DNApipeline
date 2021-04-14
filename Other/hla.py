#! /usr/bin/python3
# -*- coding: utf-8 -*-

__Version__ = "1.00"
__Author__ = "pzweuj"
__Date__ = "20210414"

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
        self.pair = runningInfo["pair"]
        self.rawdata = runningInfo["rawdata"]
        self.output = runningInfo["output"]

        self.threads = str(runningInfo["process"]["threads"])
        self.runApp = runningInfo["process"]["Other"]["HLA"]
        self.buildver = runningInfo["setting"]["Annotation"]["buildver"]

        mkdir(self.output)
        mkdir(self.output + "/tempFile")
        mkdir(self.output + "/HLA")

    # 提取HLA区域
    def extractHLA(self):
        resultsDir = self.output
        sampleID = self.sample
        pairID = self.pair
        buildver = self.buildver
        threads = self.threads
        tmpDir = self.tmpDir

        if buildver == "hg19":
            HLA = [
                "chr6:28477797-33448354",
                "chr6_apd_hap1:1-4622290",
                "chr6_cox_hap2:1-4795371",
                "chr6_dbb_hap3:1-4610396",
                "chr6_mann_hap4:1-4683263",
                "chr6_mcf_hap5:1-4833398",
                "chr6_qbl_hap6:1-4611984",
                "chr6_ssto_hap7:1-4928567"
            ]

        elif buildver == "hg38":
            HLA = [
                "chr6:28510120-33480577",
                "chr6_GL000250v2_alt:1066038-4433734",
                "chr6_GL000251v2_alt:1283988-4540572",
                "chr6_GL000252v2_alt:1063230-4372611",
                "chr6_GL000253v2_alt:1062914-4548533",
                "chr6_GL000254v2_alt:1062887-4416229",
                "chr6_GL000255v2_alt:1063190-4323464",
                "chr6_GL000256v2_alt:1106450-4577757"
            ]

        else:
            print("Cannot get buildver")
            exit()

        extractRegion = " ".join(HLA)
        
        # 当有配对样本，优先使用配对样本进行HLA分型
        if pairID != None:
            sampleID = pairID

        cmd = """
            samtools view {resultsDir}/bam/{sampleID}.bam \\
                {extractRegion} \\
                -b > {tmpDir}/{sampleID}.HLA.bam
            samtools view {resultsDir}/bam/{sampleID}.bam -bh -f 12 -@ {threads} > {tmpDir}/{sampleID}.unmapped.bam
            samtools merge {tmpDir}/{sampleID}.merge.bam {tmpDir}/{sampleID}.HLA.bam {tmpDir}/{sampleID}.unmapped.bam -@ {threads}
            samtools sort {tmpDir}/{sampleID}.merge.bam -n -@ {threads} -o {tmpDir}/{sampleID}.sort.bam
            samtools fastq {tmpDir}/{sampleID}.sort.bam \\
                -1 {tmpDir}/{sampleID}.HLA.R1.fastq \\
                -2 {tmpDir}/{sampleID}.HLA.R2.fastq -@ {threads} -s /dev/null
            rm {tmpDir}/*.bam
        """.format(threads=threads, resultsDir=resultsDir, sampleID=sampleID, extractRegion=extractRegion, tmpDir=tmpDir)
        print(cmd)
        os.system(cmd)


    # HLA-HD
    # https://www.genome.med.kyoto-u.ac.jp/HLA-HD/
    # 速度慢
    def hlahd(self):
        resultsDir = self.output
        sampleID = self.sample
        pairID = self.pair
        buildver = self.buildver
        threads = self.threads

        # 以下数据库无需指定参考基因坐标，为通用数据库，因此不写入配置文件中
        freq = "/home/bioinfo/ubuntu/software/hlahd.1.3.0/freq_data"
        dictionary = "/home/bioinfo/ubuntu/software/hlahd.1.3.0/dictionary"
        split_file = "/home/bioinfo/ubuntu/software/hlahd.1.3.0/HLA_gene.ABC.txt"

        tmpDir = resultsDir + "/tempFile/hlahd_" + sampleID
        self.tmpDir = tmpDir
        mkdir(tmpDir)
        self.extractHLA()

        if pairID != None:
            sampleID = pairID
        cmd = """
            hlahd.sh -t {threads} -m 100 -c 0.95 -f {freq} \\
                {tmpDir}/{sampleID}.HLA.R1.fastq {tmpDir}/{sampleID}.HLA.R2.fastq \\
                {split_file} {dictionary} {sampleID} {tmpDir}
            cp {tmpDir}/{sampleID}/result/{sampleID}_final.result.txt {resultsDir}/HLA/{sampleID}.hlahd.txt
        """.format(threads=threads, freq=freq, split_file=split_file, dictionary=dictionary, resultsDir=resultsDir, sampleID=sampleID, tmpDir=tmpDir)
        print(cmd)
        os.system(cmd)


    # HLAscan
    # https://github.com/SyntekabioTools/HLAscan
    # 用起来麻烦，还得跑仨次
    def hlascan(self):
        resultsDir = self.output
        sampleID = self.sample
        pairID = self.pair
        buildver = self.buildver
        threads = self.threads
        hla_scan = "/home/bioinfo/ubuntu/software/HLAscan/hla_scan_r_v2.1.4"
        hla_db = "/home/bioinfo/ubuntu/software/HLAscan/HLA-ALL.IMGT"

        tmpDir = resultsDir + "/tempFile/hlascan_" + sampleID
        self.tmpDir = tmpDir
        mkdir(tmpDir)
        self.extractHLA()

        if pairID != None:
            sampleID = pairID

        print("开始进行HLA分型")
        cmd = """
            {hla_scan} -l {tmpDir}/{sampleID}.HLA.R1.fastq \\
                -r {tmpDir}/{sampleID}.HLA.R2.fastq \\
                -t {threads} \\
                -d {hla_db} -g HLA-A > {tmpDir}/{sampleID}.HLA-A.txt
            {hla_scan} -l {tmpDir}/{sampleID}.HLA.R1.fastq \\
                -r {tmpDir}/{sampleID}.HLA.R2.fastq \\
                -t {threads} \\
                -d {hla_db} -g HLA-B > {tmpDir}/{sampleID}.HLA-B.txt
            {hla_scan} -l {tmpDir}/{sampleID}.HLA.R1.fastq \\
                -r {tmpDir}/{sampleID}.HLA.R2.fastq \\
                -t {threads} \\
                -d {hla_db} -g HLA-C > {tmpDir}/{sampleID}.HLA-C.txt
            cat {tmpDir}/{sampleID}.HLA-A.txt {tmpDir}/{sampleID}.HLA-B.txt {tmpDir}/{sampleID}.HLA-C.txt \\
                > {tmpDir}/{sampleID}.hlahd.txt
            cp {tmpDir}/{sampleID}.hlahd.txt {resultsDir}/HLA/
        """.format(hla_scan=hla_scan, hla_db=hla_db, tmpDir=tmpDir, sampleID=sampleID, threads=threads, resultsDir=resultsDir)
        print(cmd)
        os.system(cmd)

    # OptiType
    # https://github.com/FRED-2/OptiType
    # optitype线程数需要修改软件配置文件，目前默认为8
    # 速度慢，较准确，推荐使用
    def optitype(self):
        resultsDir = self.output
        sampleID = self.sample
        pairID = self.pair
        buildver = self.buildver
        threads = self.threads
        optipipe = "/home/bioinfo/ubuntu/software/OptiType-1.3.5/OptiTypePipeline.py"

        tmpDir = resultsDir + "/tempFile/optitype_" + sampleID
        self.tmpDir = tmpDir
        mkdir(tmpDir)
        self.extractHLA()

        if pairID != None:
            sampleID = pairID
        cmd = """
            python {optipipe} \\
                -i {tmpDir}/{sampleID}.HLA.R1.fastq {tmpDir}/{sampleID}.HLA.R2.fastq \\
                -d -o {tmpDir} -p {sampleID} -v
            cp {tmpDir}/{sampleID}_result.tsv {resultsDir}/HLA/{sampleID}.optitype.txt
        """.format(optipipe=optipipe, tmpDir=tmpDir, sampleID=sampleID, resultsDir=resultsDir)
        print(cmd)
        os.system(cmd)

    # seq2HLA
    # https://github.com/TRON-Bioinformatics/seq2HLA
    # 速度快，准确度稍差，建议在数据量非常大时使用
    def seq2hla(self):
        resultsDir = self.output
        sampleID = self.sample
        pairID = self.pair
        buildver = self.buildver
        threads = self.threads
        
        tmpDir = resultsDir + "/tempFile/seq2hla_" + sampleID
        self.tmpDir = tmpDir
        mkdir(tmpDir)
        self.extractHLA()

        if pairID != None:
            sampleID = pairID
        s2h = "/home/bioinfo/ubuntu/software/seq2HLA/seq2HLA.py"
        cmd = """
            python {s2h} \\
                -1 {tmpDir}/{sampleID}.HLA.R1.fastq \\
                -2 {tmpDir}/{sampleID}.HLA.R2.fastq \\
                -r {tmpDir}/{sampleID} -p {threads}
            cp {tmpDir}/{sampleID}-ClassI-class.HLAgenotype4digits {resultsDir}/HLA/{sampleID}.seq2HLA.txt
        """.format(s2h=s2h, tmpDir=tmpDir, sampleID=sampleID, threads=threads, resultsDir=resultsDir)
        print(cmd)
        os.system(cmd)

# end