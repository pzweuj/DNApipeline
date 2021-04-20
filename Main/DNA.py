#! /usr/bin/python3
# -*- coding: utf-8 -*-

__Version__ = "1.01"
__Author__ = "pzweuj"
__Date__ = "20210420"


"""
本程序为DNA自动化分析主流程，采用配置文件作为输入的方式运行程序，其中配置文件模板位于
Config文件夹下。
"""

import os
import sys
import argparse
import shutil

FUN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(FUN_DIR)
from Other.function import *
from QC.qc import QC
from Mapping.mapping import Mapping
from Mutation.snv_indel import SNV_Indel
from Mutation.cnv import CNV
from Mutation.sv import SV
from Annotation.anno import Annotation
from Other.msi import MSI
from Other.hla import HLA
from Other.tmb import TMB

def main(runInfo):
    # 基本信息获取
    runningInformation = getRunningInfo(runInfo)
    project = runningInformation["project"]
    rawdataDir = runningInformation["rawdata"]
    output = runningInformation["output"]
    sample = runningInformation["sample"]
    pair = runningInformation["pair"]
    process = runningInformation["process"]
    setting = runningInformation["setting"]
    print("使用配置 " + runInfo)
    print("运行项目 " + project)

    # 检查原始数据及重命名
    for fileName in os.listdir(rawdataDir):
        if (sample + "_") in fileName:
            print("已成功找到 " + fileName)
            if "R1" in fileName:
                os.rename(rawdataDir + "/" + fileName, rawdataDir + "/" + sample + "_R1.fastq.gz")
            elif "R2" in fileName:
                os.rename(rawdataDir + "/" + fileName, rawdataDir + "/" + sample + "_R2.fastq.gz")
            else:
                continue
        else:
            continue

    if pair != None:
        print("使用配对样本模式！")
        for fileName in os.listdir(rawdataDir):
            if (pair + "_") in fileName:
                print("已成功找到配对样本 " + fileName)
                if "R1" in fileName:
                    os.rename(rawdataDir + "/" + fileName, rawdataDir + "/" + pair + "_R1.fastq.gz")
                elif "R2" in fileName:
                    os.rename(rawdataDir + "/" + fileName, rawdataDir + "/" + pair + "_R2.fastq.gz")
                else:
                    continue
            else:
                continue

    # 运行信息获取
    threads_ = process["threads"]
    QC_ = process["QC"]
    Mapping_ = process["Mapping"]
    SnvIndel_ = process["Mutation"]["SNV_indel"]
    SV_ = process["Mutation"]["SV"]
    CNV_ = process["Mutation"]["CNV"]
    Annotation_ = process["Annotation"]
    MSI_ = process["Other"]["MSI"]
    HLA_ = process["Other"]["HLA"]
    TMB_ = process["Other"]["TMB"]

    # 质控
    # [fastp]
    if QC_ == None:
        print("根据设定不进行质控")
    else:
        QC_process = QC(runningInformation)
        print("使用 " + QC_process.runApp + " 进行质控")
        print("质控后文件输出目录： " + QC_process.output + "/cleandata")
        print("使用线程数 " + QC_process.threads)
        if QC_process.runApp == "fastp":
            QC_process.fastp()
            QC_process.fastp_filter()
        else:
            print("未找到此质控方法")

    # 比对
    # [bwa_mem]
    if Mapping_ == None:
        print("根据设定不进行比对")
    else:
        Mapping_process = Mapping(runningInformation)
        print("使用 " + Mapping_process.runApp + " 进行比对")
        print("比对后文件输出目录： " + Mapping_process.output + "/bam")
        print("使用线程数 " + Mapping_process.threads)
        if Mapping_process.runApp == "bwa_mem":
            Mapping_process.bwa_mem()
        else:
            print("未找到此比对方法")

        # 只有进行了比对才会考虑进行去重和校对
        if Mapping_process.runApp:
            # 只要进行了umi的提取，就进行gencore过滤
            if Mapping_process.umi:
                Mapping_process.gencore()
            # 此步会将比对结果bam文件进行替换
            if Mapping_process.markDups:
                Mapping_process.markDuplicates()
            # 此步会在bam文件夹下新建BQSR.bam文件
            if Mapping_process.recalibrate:
                Mapping_process.recalibrator()

        # 使用bamdst进行捕获分析，仅当设定bed文件时进行
        if Mapping_process.bed != None:
            Mapping_process.bamdst()

    # 变异检测
    ## SNV indel
    ## [gatk_m2, varscan2, pisces, freebayes, gatk_haplotypercaller]
    if SnvIndel_ == None:
        print("根据设定不进行SNV/indel检测")
    else:
        SnvIndel_process = SNV_Indel(runningInformation)
        print("使用 " + SnvIndel_process.runApp + " 进行SNV/indel检测")
        print("检测后文件输出目录： " + SnvIndel_process.output + "/vcf")
        print("使用线程数 " + SnvIndel_process.threads)
        if SnvIndel_process.runApp == "gatk_m2":
            SnvIndel_process.gatk_m2()
            # 是否过滤，此步必须要bed文件
            if SnvIndel_process.runningInfo["setting"]["Mutation"]["gatk_filter"]["run"]:
                print("进行GATK4 Mutect2过滤")
                SnvIndel_process.gatk_filter()
        elif SnvIndel_process.runApp == "freebayes":
            SnvIndel_process.freebayes()
        elif SnvIndel_process.runApp == "gatk_haplotypecaller":
            SnvIndel_process.gatk_haplotypecaller()
        elif SnvIndel_process.runApp == "pisces":
            SnvIndel_process.pisces()
        elif SnvIndel_process.runApp == "varscan" or "varscan2":
            SnvIndel_process.varscan2()
            SnvIndel_process.varscan_filter()
        else:
            print("未找到此变异检测方法")
        print("全局过滤")
        print("最低深度 " + SnvIndel_process.filtDP + "X")
        print("QUAL " + SnvIndel_process.filtQUAL)
        SnvIndel_process.filter()

    ## CNV
    ## [cnvkit]
    if CNV_ == None:
        print("根据设定不进行CNV检测")
    else:
        CNV_process = CNV(runningInformation)
        print("使用 " + CNV_process.runApp + " 进行结构变异/融合变异检测")
        print("检测后文件输出目录： " + CNV_process.output + "/sv")
        print("使用线程数 " + CNV_process.threads)
        if CNV_process.runApp == "cnvkit":
            CNV_process.cnvkit()
            CNV_process.cnvkit_filter()
        else:
            print("未找到此CNV检测方法")

    ## SV
    ## [lumpy, manta]
    if SV_ == None:
        print("根据设定不进行SV检测")
    else:
        SV_process = SV(runningInformation)
        print("使用 " + SV_process.runApp + " 进行结构变异/融合变异检测")
        print("检测后文件输出目录： " + SV_process.output + "/Fusion")
        print("使用线程数 " + SV_process.threads)
        if SV_process.runApp == "lumpy":
            SV_process.lumpy()
            SV_process.lumpy_filter()
            SV_process.sv_anno()
        elif SV_process.runApp == "manta":
            SV_process.manta()
            SV_process.manta_filter()
            SV_process.sv_anno()
        else:
            print("未找到此SV检测方法")

    # 注释
    # [annovar]
    if Annotation_ == None:
        print("根据设定不进行注释")
    elif Annotation_ == False:
        print("根据设定不进行注释")
    else:
        Annotation_process = Annotation(runningInformation)
        print("使用多工具进行注释")
        print("检测后文件输出目录： " + Annotation_process.output + "/annotation")
        print("使用线程数 " + Annotation_process.threads)
        Annotation_process.annovar()
        Annotation_process.ResultsFilter()


    # Other
    ## MSI
    ## [msisensorpro, msisensor2, msisensor_ct]
    if MSI_ == None:
        print("根据设定不进行MSI检测")
    else:
        MSI_process = MSI(runningInformation)
        print("使用 " + MSI_process.runApp + " 进行微卫星不稳定检测")
        print("检测后文件输出目录： " + MSI_process.output + "/msi")
        print("使用线程数 " + MSI_process.threads)
        if MSI_process.runApp == "msisensor2":
            MSI_process.msisensor2()
        elif MSI_process.runApp == "msisensor_pro" or MSI_process.runApp == "msisensorpro":
            MSI_process.msisensor_pro()
        elif MSI_process.runApp == "msisensor_ct":
            MSI_process.msisensor_ct()
        else:
            print("未找到此MSI分析方法")

    ## HLA
    ## [optitype, seq2hla, hlascan, hlahd]
    if HLA_ == None:
        print("根据设定不进行HLA分析")
    else:
        HLA_process = HLA(runningInformation)
        print("使用 " + HLA_process.runApp + "  进行HLA分型分析")
        print("检测后文件输出目录： " + HLA_process.output + "/HLA")
        print("使用线程数 " + HLA_process.threads)
        if HLA_process.runApp == "hlahd":
            HLA_process.hlahd()
        elif HLA_process.runApp == "optitype":
            HLA_process.optitype()
        elif HLA_process.runApp == "hlascan":
            HLA_process.hlascan()
        elif HLA_process.runApp == "seq2hla":
            HLA_process.seq2hla()
        else:
            print("未找到此HLA分析方法")


    ## TMB
    if TMB_ == None:
        print("根据设定不进行TMB计算")
    elif TMB_ == False:
        print("根据设定不进行TMB计算")
    else:
        TMB_process = TMB(runningInformation)
        print("进行TMB计算 Panel大小设定："  + str(TMB_process.panelSize))
        print("检测后文件输出目录： " + TMB_process.output + "/TMB")
        TMB_process.tmb_counter()

    # 合并结果到excel表中
    mergeResultsToExcel(output, sample)
    if pair != None:
        mergeResultsToExcel(output, pair)

    # 删除中间文件
    if os.path.exists(output + "/tempFile"):
        if runningInformation["setting"]["REMOVE_TMP"]:
            tmpDirReady = os.listdir(output + "/tempFile")
            if len(tmpDirReady) != 0:
                for t in tmpDirReady:
                    if sample in t:
                        tmpDirReadyToRemove = output + "/tempFile/" + t
                        if os.path.isdir(tmpDirReadyToRemove):
                            shutil.rmtree(tmpDirReadyToRemove)        
            tmpDirReady = os.listdir(output + "/tempFile")
            if len(tmpDirReady) == 0:
                shutil.rmtree(output + "/tempFile")

    if os.path.exists(output + "/tempFile"):
        if pair != None:
            if runningInformation["setting"]["REMOVE_TMP"]:
                tmpDirReady = os.listdir(output + "/tempFile")
                if len(tmpDirReady) != 0:
                    for t in tmpDirReady:
                        if pair in t:
                            tmpDirReadyToRemove = output + "/tempFile/" + t
                            if os.path.isdir(tmpDirReadyToRemove):
                                shutil.rmtree(tmpDirReadyToRemove)
                tmpDirReady = os.listdir(output + "/tempFile")
                if len(tmpDirReady) != 0:
                    shutil.rmtree(output + "/tempFile")


    # 报告
    if runningInformation["process"]["Report"]:
        if project == "Chem":
            print("正在生成 " + project + " 项目报告")
            runningPath = getAbsPath()
            cmd = """
                python3 {runningPath}/../Report/Report_Chem.py {resultsDir} {sampleID} {runInfo}
            """.format(runningPath=runningPath, resultsDir=output, sampleID=sample, runInfo=runInfo)
            print(cmd)
            os.system(cmd)
        
        if project == "93Genes":
            pass

        if project == "Pancancer":
            pass

        if project == "Lung":
            print("正在生成 " + project + " 项目报告")
            runningPath = getAbsPath()
            cmd = """
                python3 {runningPath}/../Report/Report_Lung.py {resultsDir} {sampleID} {runInfo}
            """.format(runningPath=runningPath, resultsDir=output, sampleID=sample, runInfo=runInfo)
            print(cmd)
            os.system(cmd)

        if project == "Colon":
            pass

        if project == "BRCA":
            pass

    else:
        print("根据设定不进行报告输出")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="DNA Pipeline",
        prog="DNA Pipeline",
        usage="python3 DNA.py -i <Config Yaml File>")
    group = parser.add_mutually_exclusive_group()
    parser.add_argument("-v", "--version", action="version",
        version="Version 1.01 20210420")
    parser.add_argument("-i", "--input", type=str,
        help="输入配置文件")

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()
    main(runInfo=args.input)

# 测试
# NCCL = "/home/bioinfo/ubuntu/bin/DNApipeline-main/Config/NCCL.yaml"
# main(NCCL)
