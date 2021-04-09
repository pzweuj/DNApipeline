#! /usr/bin/python3
# -*- coding: utf-8 -*-

__Version__ = "0.01"
__Author__ = "pzweuj"
__Date__ = "20210402"

import os
import sys

FUN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(FUN_DIR)
from Other.function import mkdir


class LOH(object):
    """
    LOH检测模块，用于检测杂合性缺失
    """

    def __init__(self, runningInfo):
        self.runningInfo = runningInfo
        self.sample = runningInfo["sample"]
        self.rawdata = runningInfo["rawdata"]
        self.output = runningInfo["output"]

        self.threads = str(runningInfo["process"]["threads"])
        # self.runApp = runningInfo["process"]["Other"]["HLA"]
        self.buildver = runningInfo["setting"]["Annotation"]["buildver"]

        mkdir(self.output)
        mkdir(self.output + "/tempFile")
        mkdir(self.output + "/HLA")

    # LOHHLA
    # https://bitbucket.org/mcgranahanlab/lohhla/src/master/
    def lohhla(self):
        pass