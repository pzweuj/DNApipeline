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


class Annotation(object):
    """
    注释模块
    后续结果将全部改放到中间文件夹中，再使用自建脚本处理为excel表格后
    再放置到annotation文件夹中，表格格式保持一致
    """
    def __init__(self, runningInfo):
        self.runningInfo = runningInfo
        self.sample = runningInfo["sample"]
        self.rawdata = runningInfo["rawdata"]
        self.output = runningInfo["output"]

        self.threads = str(runningInfo["process"]["threads"])
        self.runApp = runningInfo["process"]["Annotation"]

        mkdir(self.output)
        mkdir(self.output + "/tempFile")
        mkdir(self.output + "/annotation")

    # annovar
    # https://annovar.openbioinformatics.org/en/latest/
    def annovar(self):
        humandb = self.runningInfo["setting"]["Annotation"]["humandb"]
        buildver = self.runningInfo["setting"]["Annotation"]["buildver"]
        resultsDir = self.output
        sampleID = self.sample
        threads = self.threads

        cmd = """
            convert2annovar.pl -format vcf4 {resultsDir}/vcf/{sampleID}.vcf \\
                --includeinfo > {resultsDir}/tempFile/{sampleID}.avinput
            table_annovar.pl {resultsDir}/tempFile/{sampleID}.avinput \\
                {humandb} -buildver {buildver} \\
                -out {resultsDir}/annotation/{sampleID} -remove \\
                -protocol refGene,cytoBand,avsnp150,gnomad211_genome,clinvar_20210131,cosmic70,dbnsfp41a \\
                -operation g,r,f,f,f,f,f \\
                -nastring . -thread {threads} -otherinfo
        """.format(resultsDir=resultsDir, sampleID=sampleID, humandb=humandb, threads=threads, buildver=buildver)
        print(cmd)
        os.system(cmd)        

    def snpeff(self):
        pass

    def vep(self):
        pass