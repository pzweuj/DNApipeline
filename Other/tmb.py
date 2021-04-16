#! /usr/bin/python3
# -*- coding: utf-8 -*-

__Version__ = "0.12"
__Author__ = "pzweuj"
__Date__ = "20210409"

import os
import sys

FUN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(FUN_DIR)
from Other.function import mkdir

class TMB(object):
    """
    TMB计算模块
    """
    def __init__(self, runningInfo):
        self.runningInfo = runningInfo
        self.sample = runningInfo["sample"]
        self.rawdata = runningInfo["rawdata"]
        self.output = runningInfo["output"]

        self.runApp = runningInfo["process"]["Other"]["TMB"]
        self.panelSize = runningInfo["setting"]["Other"]["PanelSize"]
        self.buildver = runningInfo["setting"]["Annotation"]["buildver"]

        mkdir(self.output)
        mkdir(self.output + "/tempFile")
        mkdir(self.output + "/TMB")

    # 参考
    # https://www.accessdata.fda.gov/cdrh_docs/pdf17/P170019B.pdf
    def tmb_counter(self):
        resultsDir = self.output
        sampleID = self.sample
        panelSize = self.panelSize

        tmpDir = resultsDir + "/tempFile/TMB_" + sampleID
        mkdir(tmpDir)

        polyDB_filter = 0.01
        cosmic_filter = 5
        AF_filter = 0.05

        annovarFile = open(resultsDir + "/annotation/" + sampleID + ".Anno.txt", "r")
        resultsFile = open(tmpDir + "/" + sampleID + ".tmb.txt", "w")
        resultsFile.write("# " + sampleID + "\n")
        resultsFile.write("panelSize\t" + str(panelSize) + "\n")
        resultsFile.write("过滤条件：\n")
        resultsFile.write("人群频率(需小于此值)\t" + str(polyDB_filter) + "\n")
        resultsFile.write("cosmic出现次数(需小于此值)\t" + str(cosmic_filter) + "\n")
        resultsFile.write("突变丰度(需大于此值)\t" + str(AF_filter) + "\n")


        # 计数
        n = 0
        for line in annovarFile:
            if not line.startswith("Chr\t"):
                lines = line.replace("\n", "").split("\t")
                func = lines[15]
                cosmic = lines[31]
                polyDB = lines[18]
                polyDB_eas = lines[22]
                Jax = lines[28]
                Civic = lines[29]
                Onco = lines[30]
                DP = int(lines[13])
                AD = int(lines[14])
                AF = float(lines[10].replace("%", "")) / 100

                # 处理内容
                if cosmic == "-":
                    cosmic_count = 0
                else:
                    cosmic_count = 0
                    cos_occur = cosmic.split(";")[1].split("=")[1].split(",")
                    for o in cos_occur:
                        occur = int(o.split("(")[0])
                        cosmic_count += occur

                if polyDB == "-" or polyDB == ".":
                    polyDB = 0
                else:
                    polyDB = float(polyDB)

                if polyDB_eas == "-" or polyDB_eas == ".":
                    polyDB_eas = 0
                else:
                    polyDB_eas = float(polyDB_eas)
                
                # 过滤
                # 滤去AF小于等于5%的点
                if AF <= AF_filter:
                    continue

                # 根据人群频率过滤
                if polyDB >= polyDB_filter:
                    continue
                if polyDB_eas >= polyDB_filter:
                    continue

                # 根据cosmic滤去肿瘤热点
                if cosmic_count >= cosmic_filter:
                    continue

                # 仅保留exonic与splicing
                if "ncRNA" in func:
                	continue
                elif "exonic" in func:
                    pass
                elif "splicing" in func:
                    pass
                else:
                    continue

                # 过滤热点
                if Jax != "-":
                	continue
                if Civic != "-":
                	continue
                if Onco != "-":
                	continue

                # 通过过滤
                print(line)
                n += 1
        annovarFile.close()
        panelSize_mb = float(panelSize) / 1000000
        TMB_out = n / panelSize_mb
        resultsFile.write("TMB\t" + "%.2f" % TMB_out + "\n")
        print("TMB: ", TMB_out)


        # 在相同癌种中的排位，数据缺乏，暂缺

        resultsFile.close()
        cmd = """
            cp {tmpDir}/{sampleID}.tmb.txt {resultsDir}/TMB/
        """.format(tmpDir=tmpDir, sampleID=sampleID, resultsDir=resultsDir)
        os.system(cmd)
