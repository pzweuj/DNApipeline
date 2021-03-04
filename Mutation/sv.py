#! /usr/bin/python3
# -*- coding: utf-8 -*-

__Version__ = "0.17"
__Author__ = "pzweuj"
__Date__ = "20210304"


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

    # 自编manta注释
    # 对数据库进行解析，数据库使用refFlat，可通过ucsc直接下载
    def checkFusionHotSpot(self):
        database = self.runningInfo["setting"]["Annotation"]["refFlat"]
        fusionStrand = self.svStrand
        chrom = self.svChrom
        breakPoint = self.breakPoint

        db = open(database, "r")
        output = []
        ts_output = []
        for line in db:
            l = line.split("\n")[0].split("\t")
            gene = l[0]
            ts = l[1]
            chrom_db = l[2]
            strand = l[3]
            gene_start = l[4]
            gene_end = l[5]
            cds_start = l[6]
            cds_end = l[7]
            exon_nums = l[8]
            exon_starts = l[9]
            exon_ends = l[10]

            # 寻找注释
            if chrom == chrom_db:
                # 找到目标基因
                if (int(breakPoint) >= int(gene_start)) and (int(breakPoint) <= int(gene_end)):

                    exon_starts_list = exon_starts.split(",")[0: -1]
                    exon_ends_list = exon_ends.split(",")[0: -1]

                    # 将外显子起止点形成列表
                    exon_checkpoint = []
                    i = 0
                    while i < int(exon_nums):
                        exon_checkpoint.append(int(exon_starts_list[i]))
                        exon_checkpoint.append(int(exon_ends_list[i]))
                        i += 1

                    # 判断目标基因的转录方向
                    if strand == "+":
                        if int(breakPoint) < exon_checkpoint[0]:
                            if fusionStrand == "+":
                                exon_out = gene + "_5UTR"
                            elif fusionStrand == "-":
                                exon_out = gene + "_exon1"
                            else:
                                exon_out = gene + "_?"

                        elif int(breakPoint) > exon_checkpoint[-1]:
                            if fusionStrand == "+":
                                exon_out = gene + "_exon" + exon_nums
                            elif fusionStrand == "-":
                                exon_out = gene + "_3UTR"
                            else:
                                exon_out = gene + "_?"

                        else:
                            n = 0
                            while n < (2 * int(exon_nums) - 1):
                                if (int(breakPoint) >= exon_checkpoint[n]) and (int(breakPoint) <= exon_checkpoint[n+1]):
                                    checkN = n
                                n += 1

                            if fusionStrand == "+":
                                exon_out = gene + "_exon" + str(checkN // 2 + 1)
                            elif fusionStrand == "-":
                                exon_out = gene + "_exon" + str((checkN + 1) // 2 + 1)
                            else:
                                exon_out = gene + "_?"


                    if strand == "-":
                        if int(breakPoint) < exon_checkpoint[0]:
                            if fusionStrand == "+":
                                exon_out = gene + "_3UTR"
                            elif fusionStrand == "-":
                                exon_out = gene + "_exon" + exon_nums
                            else:
                                exon_out = gene + "_?"


                        elif int(breakPoint) > exon_checkpoint[-1]:
                            if fusionStrand == "+":
                                exon_out = gene + "_exon1"
                            elif fusionStrand == "-":
                                exon_out = gene + "_5UTR"
                            else:
                                exon_out = gene + "_?"

                        else:
                            n = 0
                            while n < (2 * int(exon_nums) - 1):
                                if (int(breakPoint) >= exon_checkpoint[n]) and (int(breakPoint) <= exon_checkpoint[n+1]):
                                    checkN = n
                                n += 1

                            if fusionStrand == "+":
                                exon_out = gene + "_exon" + str(int(exon_nums) - checkN // 2)
                            elif fusionStrand == "-":
                                exon_out = gene + "_exon" + str(int(exon_nums) - (checkN + 1) // 2)
                            else:
                                exon_out = gene + "_?"

                    output.append(exon_out)
                    ts_output.append(ts)

        if not output:
            output.append("GeneUnknown_exon?")        
        if not ts_output:
            ts_output.append("Transcript?")
        db.close()
        return [output, ts_output]

    def manta_anno(self):
        sampleID = self.sample
        resultsDir = self.output
        refFlat = self.runningInfo["setting"]["Annotation"]["refFlat"]

        svFile = open(resultsDir + "/sv/" + sampleID + ".sv.vcf", "r")
        svAnno = open(resultsDir + "/sv/" + sampleID + ".ann.txt", "w")
        svAnno.write("chrom1\tbreakpoint1\tgene1\tchrom2\tbreakpoint2\tgene2\tfusionType\tAlt\tgeneSymbol\tPR\tSR\tDP\tVAF\tExon\tTranscript\n")
        for line in svFile:
            if line.startswith("#"):
                continue
            else:
                lineAfterSplit = line.split("\t")
                chrom = lineAfterSplit[0]
                Pos = lineAfterSplit[1]
                ID = lineAfterSplit[2]
                Ref = lineAfterSplit[3]
                Alt = lineAfterSplit[4]
                Qual = lineAfterSplit[5]
                Filter = lineAfterSplit[6]
                Info = lineAfterSplit[7]
                Format = lineAfterSplit[8]
                Sample = lineAfterSplit[9]

                if chrom == "chrM":
                    continue

                if ":" in Format:
                    S = Sample.split(":")
                    PR = S[0].split(",")
                    SR = S[1].split(",")
                    PR_ref = int(PR[0])
                    PR_alt = int(PR[1])
                    SR_ref = int(SR[0])
                    SR_alt = int(SR[1])
                else:
                    PR = Sample.split(",")
                    PR_ref = int(PR[0])
                    PR_alt = int(PR[1])
                    SR_ref = SR_alt = 0

                DP = PR_ref + SR_ref + PR_alt + SR_alt
                VAF = "%.2f" % (100 * float(PR_alt + SR_alt) / DP) + "%"

                if ("[" in Alt) or ("]" in Alt):
                    if chrom in Alt:
                        fusionType = "Inversion"
                    else:
                        fusionType = "Translocation"

                    # G    [chr1:1111111[G
                    if Alt[0] == "[":
                        chrom2 = Alt.split(":")[0].split("[")[1]
                        breakpoint2 = Alt.split(":")[1].split("[")[0]
                        strand = "--"
                    
                    # G    G[chr1:1111111[
                    elif Alt[-1] == "[":
                        chrom2 = Alt.split(":")[0].split("[")[1]
                        breakpoint2 = Alt.split(":")[1].split("[")[0]
                        strand = "+-"

                    # G    ]chr1:1111111]G
                    elif Alt[0] == "]":
                        chrom2 = Alt.split(":")[0].split("]")[1]
                        breakpoint2 = Alt.split(":")[1].split("]")[0]
                        strand = "-+"

                    # G    G]chr1:1111111]
                    elif Alt[-1] == "]":
                        chrom2 = Alt.split(":")[0].split("]")[1]
                        breakpoint2 = Alt.split(":")[1].split("]")[0]
                        strand = "++"

                    else:
                        chrom2 = "Unknown"
                        breakpoint2 = "Unknown"
                        strand = "Unknown"
                        print("can not analysis: " + ID)

                    self.database = refFlat
                    self.svStrand = strand[0]
                    self.svChrom = chrom
                    self.breakPoint = Pos
                    anno1 = self.checkFusionHotSpot()
                    if chrom2 != "Unknown":
                        self.svStrand = strand[1]
                        self.svChrom = chrom2
                        self.breakPoint = breakpoint2
                        anno2 = self.checkFusionHotSpot()
                    else:
                        anno2 = [["GeneUnknown_exon?"], ["Transcript?"]]

                    gene1 = anno1[0][0].split("_")[0]
                    gene2 = anno2[0][0].split("_")[0]
                    exon_output = ",".join(anno1[0]) + "|" + ",".join(anno2[0])
                    transcript_output = ",".join(anno1[1]) + "|" + ",".join(anno2[1])
                    outputStringList = [chrom, Pos, gene1, chrom2, breakpoint2, gene2, fusionType, Alt, gene1 + "-" + gene2, str(PR_alt), str(SR_alt), str(DP), VAF, exon_output, transcript_output]

                    outputString = "\t".join(outputStringList)
                    print(outputString)
                    svAnno.write(outputString + "\n")
        svFile.close()
        svAnno.close()



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