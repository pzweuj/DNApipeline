#! /usr/bin/python3
# -*- coding: utf-8 -*-

__Version__ = "0.1"
__Author__ = "pzweuj"
__Date__ = "20210308"

import os
import sys

FUN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(FUN_DIR)
from Other.function import mkdir

class TMB(object):
    """
    TMB计算模块
    """
    def __init__(self):
        self.runningInfo = runningInfo
        self.sample = runningInfo["sample"]
        self.rawdata = runningInfo["rawdata"]
        self.output = runningInfo["output"]

        self.runApp = runningInfo["process"]["Other"]["TMB"]
        self.panelSize = runningInfo["setting"]["Other"]["PanelSize"]

        mkdir(self.output)
        mkdir(self.output + "/tempFile")
        mkdir(self.output + "/TMB")

    # 有待测试
    # https://github.com/bioinfo-pf-curie/TMB
    # https://github.com/Clinical-Genomics-Lund/TMB-Calc
    # https://github.com/fanyucai1/TMB
    def tmb_counter(self):
        resultsDir = self.output
        sampleID = self.sample

        tmpDir = resultsDir + "/tempFile/TMB_" + sampleID
        mkdir(tmpDir)

        annovarFile = open(resultsDir + "/annotation/" + sampleID + "anno.txt", "r")
        for line in annovarFile:
        	if not line.startswith("#"):
        		pass