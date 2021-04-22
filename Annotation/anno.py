#! /usr/bin/python3
# -*- coding: utf-8 -*-

__Version__ = "1.0"
__Author__ = "pzweuj"
__Date__ = "20210416"

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
    注释流程不再可选比对工具，而是改用先使用snpEff注释，再使用annovar注释，最后对结果进行整理的模式
    """
    def __init__(self, runningInfo):
        self.runningInfo = runningInfo
        self.sample = runningInfo["sample"]
        self.pair = runningInfo["pair"]
        self.rawdata = runningInfo["rawdata"]
        self.output = runningInfo["output"]

        self.threads = str(runningInfo["process"]["threads"])

        mkdir(self.output)
        mkdir(self.output + "/tempFile")
        mkdir(self.output + "/annotation")


    # snpEff
    # https://pcingola.github.io/SnpEff/
    def snpeff(self):
        humandb = self.runningInfo["setting"]["Annotation"]["humandb"]
        buildver = self.runningInfo["setting"]["Annotation"]["buildver"]
        resultsDir = self.output
        sampleID = self.sample
        pairID = self.pair

        if pairID == None:
            vcfFile = resultsDir + "/vcf/" + sampleID + ".vcf"
        else:
            vcfFile = resultsDir + "/vcf/" + sampleID + ".filter.vcf"        

        tmpDir = resultsDir + "/tempFile/snpeff_" + sampleID
        mkdir(tmpDir)
        cmd = """
            java -jar /home/bioinfo/ubuntu/software/snpEff/snpEff.jar \\
                -c /home/bioinfo/ubuntu/software/snpEff/snpEff.config \\
                -s {tmpDir}/{sampleID}.summary.html \\
                {buildver} {vcfFile} > {tmpDir}/{sampleID}.snpeff.vcf
            cp {tmpDir}/{sampleID}.snpeff.vcf {resultsDir}/annotation/
        """.format(buildver=buildver, vcfFile=vcfFile, tmpDir=tmpDir, sampleID=sampleID, resultsDir=resultsDir)
        print(cmd)
        os.system(cmd)
    
    # annovar
    # https://annovar.openbioinformatics.org/en/latest/
    def annovar(self):
        humandb = self.runningInfo["setting"]["Annotation"]["humandb"]
        buildver = self.runningInfo["setting"]["Annotation"]["buildver"]
        resultsDir = self.output
        sampleID = self.sample
        pairID = self.pair
        threads = self.threads

        tmpDir = resultsDir + "/tempFile/annovar_" + sampleID
        mkdir(tmpDir)

        self.snpeff()

        cmd = """
            convert2annovar.pl -format vcf4 \\
                {resultsDir}/annotation/{sampleID}.snpeff.vcf \\
                --includeinfo > {tmpDir}/{sampleID}.avinput
            table_annovar.pl {tmpDir}/{sampleID}.avinput \\
                {humandb} -buildver {buildver} \\
                -out {tmpDir}/{sampleID} -remove \\
                -protocol refGene,avsnp150,gnomad211_genome,clinvar_20210308,JaxCkb,Civic,OncoKB,dbnsfp41a,cosmic92_coding,intervar_20180118 \\
                -operation g,f,f,f,f,f,f,f,f,f \\
                -nastring - -thread {threads} -otherinfo
            cp {tmpDir}/{sampleID}.{buildver}_multianno.txt {resultsDir}/annotation/
        """.format(tmpDir=tmpDir, resultsDir=resultsDir, sampleID=sampleID, humandb=humandb, threads=threads, buildver=buildver)
        print(cmd)
        os.system(cmd)
    

    # 转录本字典
    # hg19与hg38均适用
    def refTranscript(self):
        refT = open("/home/bioinfo/ubuntu/database/hg19/LRG/refTransript.txt", "r")
        refDict = {}
        for r in refT:
            if not r.startswith("#"):
                gene = r.split("\t")[0]
                ts = r.split("\t")[1]
                refDict[gene] = ts
        refT.close()
        return refDict

    # annovar结果过滤与标准化
    def ResultsFilter(self):
        buildver = self.runningInfo["setting"]["Annotation"]["buildver"]
        resultsDir = self.output
        sampleID = self.sample

        annvarResultsFile = """{resultsDir}/annotation/{sampleID}.{buildver}_multianno.txt""".format(resultsDir=resultsDir, sampleID=sampleID, buildver=buildver)
        annvarFixFile = resultsDir + "/annotation/" + sampleID + ".Anno.txt"

        annvarResults = open(annvarResultsFile, "r")
        results = open(annvarFixFile, "w")
        results.write("Chr\tStart\tEnd\tRef\tAlt\tGene\tType\tTranscript\tcHGVS\tpHGVS\tVAF\tConsequence\t"\
            "AffectedExon\tDepth\tAlt_AD\tFunc.refGene\tAAChange.refGene\tavsnp150\tAF\tAF_popmax\tAF_male\t"\
            "AF_female\tAF_eas\tCLNDN\tCLNSIG\tInterVar_automated\tInterVar_sig\tJax_Ckb_Variants_Summary\tJax_Ckb_Drug_Summary\t"\
            "Civic\tOncoKB\tcosmic92_coding\tDamagePredCount\tSIFT_pred\tPolyphen2_HDIV_pred\tPolyphen2_HVAR_pred\t"\
            "LRT_pred\tMutationTaster_pred\tMutationAssessor_pred\tFATHMM_pred\tPROVEAN_pred\t"\
            "M-CAP_pred\tREVEL_score\n")

        # 生成字典keys
        lineDict = {}
        for line in annvarResults:
            if line.startswith("Chr\tStart"):
                lines = line.replace("\n", "").split("\t")
                for k in lines:
                    lineDict[k] = "-"
        annvarResults.close()

        # 再次打开
        annvarResults = open(annvarResultsFile, "r")
        for line in annvarResults:
            if not line.startswith("Chr\tStart"):
                lines = line.replace("\n", "").split("\t")
                k = list(lineDict.keys())
                
                # 生成字典
                for i in range(len(lines)):
                    lineDict[k[i]] = lines[i]
                
                # 开始处理
                # 处理Insertion状态下的End
                if lineDict["Ref"] == "-":
                    lineDict["Type"] = "Insertion"
                    lineDict["End"] = "-"
                elif lineDict["Alt"] == "-":
                    lineDict["Type"] = "Deletion"
                elif len(lineDict["Ref"]) == 1 and len(lineDict["Alt"]) == 1:
                    lineDict["Type"] = "SNV"
                else:
                    lineDict["Type"] = "Complex"

                # 处理基因，使用snpeff的注释信息
                otif8 = lineDict["Otherinfo8"].split(";")
                snpeff_anno = "ANN=?|N|N|N|N|transcript|N|N|N|N|N|N|N|N|N|"
                for o8 in otif8:
                    if o8.startswith("ANN="):
                        snpeff_anno = o8
                snpeff_annos = snpeff_anno.replace("ANN=", "").split(",")
                
                # 此处需使用固定转录本，未完成
                if len(snpeff_annos) > 1:
                    snpeff_select = snpeff_annos[0]
                    gene = snpeff_select.split("|")[3]
                    refTrans = self.refTranscript()
                    try:
                        ts = refTrans[gene]
                        for sa in snpeff_annos:
                            if ts.split(".")[0] in snpeff_annos:
                                if sa.split("|")[6] != "":
                                    snpeff_select = sa
                    except:
                        snpeff_select = snpeff_select
                else:
                    snpeff_select = snpeff_annos[0]
                s = snpeff_select.split("|")
                for ss in range(len(s)):
                    if s[ss] == "":
                        s[ss] = "-"

                lineDict["Gene"] = s[3]
                lineDict["Transcript"] = s[6]
                lineDict["cHGVS"] = s[9]
                lineDict["pHGVS"] = s[10]
                lineDict["AffectedExon"] = s[8]

                # 处理Consequence
                # print(s[1])
                if lineDict["Func.refGene"] == "exonic" and lineDict["Type"] == "Complex":
                    lineDict["Consequence"] = "Complex_mutation"
                elif "missense" in s[1]:
                    lineDict["Consequence"] = "Missense_substitution"
                elif "splice" in s[1]:
                    lineDict["Consequence"] = "Splice_Site_mutation"
                elif "synonymous" in s[1]:
                    lineDict["Consequence"] = "Synonymous_substitution"
                elif "inframe_deletion" in s[1]:
                    lineDict["Consequence"] = "Inframe_deletion"
                elif "inframe_insertion" in s[1]:
                    lineDict["Consequence"] = "Inframe_insertion"
                elif "frameshift" in s[1]:
                    if lineDict["Type"] == "Insertion":
                        lineDict["Consequence"] = "Frameshift_insertion"
                    elif lineDict["Type"] == "Deletion":
                        lineDict["Consequence"] = "Frameshift_deletion"
                    else:
                        lineDict["Consequence"] = "Frameshift_variant"
                elif "stop_gained" in s[1]:
                    lineDict["Consequence"] = "Nonsense_substitution"
                else:
                    lineDict["Consequence"] = "Other"

                # 处理VAF
                otif10 = lineDict["Otherinfo10"].split(":")
                otif9 = lineDict["Otherinfo9"].split(":")        
                format_zip = list(zip(otif9, otif10))
                format_dict = {}
                for f in format_zip:
                    format_dict[f[0]] = f[1]

                lineDict["DP"] = format_dict["DP"]
                lineDict["Alt_AD"] = format_dict["AD"]
                if not "," in lineDict["Alt_AD"]:
                    lineDict["Alt_AD"] = lineDict["Alt_AD"]
                else:
                    lineDict["Alt_AD"] = lineDict["Alt_AD"].split(",")[1]
                lineDict["Ref_AD"] = str(int(lineDict["DP"]) - int(lineDict["Alt_AD"]))

                try:
                    AF = "%.2f" % ((float(lineDict["Alt_AD"]) / float(lineDict["DP"])) * 100) + "%"
                except Exception:
                    AF = "-"
                lineDict["VAF"] = AF

                # 合并证据等级
                lineDict["InterVar_sig"] = "PVS1=" + lineDict["PVS1"] + ";" + \
                    "PS1=" + lineDict["PS1"] + ";" + \
                    "PS2=" + lineDict["PS2"] + ";" + \
                    "PS3=" + lineDict["PS3"] + ";" + \
                    "PS4=" + lineDict["PS4"] + ";" + \
                    "PM1=" + lineDict["PM1"] + ";" + \
                    "PM2=" + lineDict["PM2"] + ";" + \
                    "PM3=" + lineDict["PM3"] + ";" + \
                    "PM4=" + lineDict["PM4"] + ";" + \
                    "PM5=" + lineDict["PM5"] + ";" + \
                    "PM6=" + lineDict["PM6"] + ";" + \
                    "PP1=" + lineDict["PP1"] + ";" + \
                    "PP2=" + lineDict["PP2"] + ";" + \
                    "PP3=" + lineDict["PP3"] + ";" + \
                    "PP4=" + lineDict["PP4"] + ";" + \
                    "PP5=" + lineDict["PP5"] + ";" + \
                    "BA1=" + lineDict["BA1"] + ";" + \
                    "BS1=" + lineDict["BS1"] + ";" + \
                    "BS2=" + lineDict["BS2"] + ";" + \
                    "BS3=" + lineDict["BS3"] + ";" + \
                    "BS4=" + lineDict["BS4"] + ";" + \
                    "BP1=" + lineDict["BP1"] + ";" + \
                    "BP2=" + lineDict["BP2"] + ";" + \
                    "BP3=" + lineDict["BP3"] + ";" + \
                    "BP4=" + lineDict["BP4"] + ";" + \
                    "BP5=" + lineDict["BP5"] + ";" + \
                    "BP6=" + lineDict["BP6"] + ";" + \
                    "BP7=" + lineDict["BP7"]

                # 整理输出结果
                output = [lineDict["Chr"], lineDict["Start"], lineDict["End"], lineDict["Ref"],
                    lineDict["Alt"], lineDict["Gene"], lineDict["Type"], lineDict["Transcript"],
                    lineDict["cHGVS"], lineDict["pHGVS"], lineDict["VAF"], lineDict["Consequence"],
                    lineDict["AffectedExon"], lineDict["DP"], lineDict["Alt_AD"], lineDict["Func.refGene"],
                    lineDict["AAChange.refGene"], lineDict["avsnp150"], lineDict["AF"], lineDict["AF_popmax"], lineDict["AF_male"], 
                    lineDict["AF_female"], lineDict["AF_eas"], lineDict["CLNDN"], lineDict["CLNSIG"],
                    lineDict["InterVar_automated"], lineDict["InterVar_sig"], lineDict["Jax_Ckb_Variants_Summary"],
                    lineDict["Jax_Ckb_Drug_Summary"], lineDict["Civic"], lineDict["OncoKB"], lineDict["cosmic92_coding"],
                    lineDict["DamagePredCount"], lineDict["SIFT_pred"], lineDict["Polyphen2_HDIV_pred"],
                    lineDict["Polyphen2_HVAR_pred"], lineDict["LRT_pred"], lineDict["MutationTaster_pred"],
                    lineDict["MutationAssessor_pred"], lineDict["FATHMM_pred"], lineDict["PROVEAN_pred"],
                    lineDict["M-CAP_pred"], lineDict["REVEL_score"]
                ]

                print("\t".join(output))
                results.write("\t".join(output) + "\n")
        print("完成annovar结果过滤与格式调整")

# end