#! /usr/bin/python3
# -*- coding: utf-8 -*-

__Version__ = "1.00"
__Author__ = "pzweuj"
__Date__ = "20210420"

import os
import sys

FUN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(FUN_DIR)
from Other.function import mkdir

class SNV_Indel(object):
    """
    SNV与indel检测模块，最终生成vcf文件
    为了流程的一致，每个caller最终生成的结果在vcf文件夹中命名为sample.vcf
    """
    def __init__(self, runningInfo):
        self.runningInfo = runningInfo
        self.sample = runningInfo["sample"]
        self.pair = runningInfo["pair"]
        self.rawdata = runningInfo["rawdata"]
        self.output = runningInfo["output"]
        
        self.threads = str(runningInfo["process"]["threads"])
        self.runApp = runningInfo["process"]["Mutation"]["SNV_indel"]

        self.databases = runningInfo["setting"]["Mapping"]["databases"]
        self.reference = runningInfo["setting"]["Mapping"]["reference"]
        self.bed = runningInfo["setting"]["Mutation"]["Bed"]
        self.filtDP = str(runningInfo["setting"]["Mutation"]["filter"]["DP"])
        self.filtQUAL = str(runningInfo["setting"]["Mutation"]["filter"]["QUAL"])
        self.MAF = str(runningInfo["setting"]["Mutation"]["filter"]["MAF"])

        mkdir(self.output)
        mkdir(self.output + "/tempFile")
        mkdir(self.output + "/vcf")

    # GATK4
    # https://github.com/broadinstitute/gatk
    def gatk_m2(self):
        reference = self.reference
        gnomad = self.databases + "/" + self.runningInfo["setting"]["Mutation"]["gnomad"]
        resultsDir = self.output
        sampleID = self.sample
        pairID = self.pair
        pon = self.runningInfo["setting"]["Mutation"]["pon"]
        bedFile = self.bed
        threads = self.threads

        checkBQSR = sampleID + ".BQSR.bam"
        if checkBQSR in os.listdir(resultsDir + "/bam"):
            bamFile = checkBQSR
        else:
            bamFile = sampleID + ".bam"

        if pairID != None:
            pairBQSR = pairID + ".BQSR.bam"
            if pairBQSR in os.listdir(resultsDir + "/bam"):
                pairFile = pairBQSR
            else:
                pairFile = pairID + ".bam"

        tmpDir = resultsDir + "/tempFile/gatk_" + sampleID
        mkdir(tmpDir)

        # 区分有无pon的情况，应可用配对样本或多个正常样本建立pon
        if pon == None:
            pon = "null"

        if bedFile == None:
            bedFile = "null"
        
        if pairID == None:
            cmd = """
                gatk Mutect2 \\
                    -R {reference} \\
                    -I {resultsDir}/bam/{bamFile} \\
                    -O {tmpDir}/{sampleID}.m2.vcf \\
                    -tumor {sampleID} \\
                    --germline-resource {gnomad} \\
                    -pon {pon} \\
                    --native-pair-hmm-threads {threads} \\
                    -L {bedFile} \\
                    -A Coverage -A GenotypeSummaries \\
                    --genotype-germline-sites false \\
                    --max-reads-per-alignment-start 0
                cp {tmpDir}/{sampleID}.m2.vcf {resultsDir}/vcf/{sampleID}.vcf
            """.format(tmpDir=tmpDir, bedFile=bedFile, pon=pon, reference=reference, resultsDir=resultsDir, sampleID=sampleID, gnomad=gnomad, threads=threads, bamFile=bamFile)
            print(cmd)
            os.system(cmd)

        else:
            cmd = """
                gatk Mutect2 \\
                    -R {reference} \\
                    -I {resultsDir}/bam/{bamFile} -tumor {sampleID} \\
                    -I {resultsDir}/bam/{pairFile} -normal {pairID} \\
                    -O {tmpDir}/{sampleID}_{pairID}.m2.vcf \\
                    --germline-resource {gnomad} \\
                    -pon {pon} \\
                    --native-pair-hmm-threads {threads} \\
                    -L {bedFile} \\
                    -A Coverage -A GenotypeSummaries \\
                    --genotype-germline-sites false \\
                    --max-reads-per-alignment-start 0
                cp {tmpDir}/{sampleID}_{pairID}.m2.vcf {resultsDir}/vcf/{sampleID}_{pairID}.vcf
            """.format(tmpDir=tmpDir, bedFile=bedFile, pon=pon, reference=reference, resultsDir=resultsDir, sampleID=sampleID, gnomad=gnomad, threads=threads, bamFile=bamFile, pairFile=pairFile, pairID=pairID)
            print(cmd)
            os.system(cmd)

    # GATK4 single sample only
    def gatk_haplotypecaller(self):
        reference = self.reference
        resultsDir = self.output
        sampleID = self.sample
        bedFile = self.bed
        threads = self.threads
        bedFile = self.bed

        tmpDir = resultsDir + "/tempFile/HaplotypeCaller_" + sampleID
        mkdir(tmpDir)

        if bedFile == None:
            bedFile = "null"
        cmd = """
            gatk HaplotypeCaller \\
                -R {reference} \\
                -I {resultsDir}/bam/{sampleID}.bam \\
                -O {tmpDir}/{sampleID}.htc.vcf \\
                -L {bedFile} \\
                --native-pair-hmm-threads {threads}
            cp {tmpDir}/{sampleID}.htc.vcf {resultsDir}/vcf/{sampleID}.vcf
        """.format(reference=reference, resultsDir=resultsDir, tmpDir=tmpDir, sampleID=sampleID, bedFile=bedFile, threads=threads)
        print(cmd)
        os.system(cmd)

    # VarScan2
    # http://varscan.sourceforge.net/
    def varscan2(self):
        reference = self.reference
        resultsDir = self.output
        sampleID = self.sample
        pairID = self.pair
        bedFile = self.bed
        threads = self.threads
        MAF = self.MAF
        DP = self.filtDP
        minCov = str(int(int(DP) * float(MAF)))

        tmpDir = resultsDir + "/tempFile/varscan2_" + sampleID
        mkdir(tmpDir)
        vs = "/home/bioinfo/ubuntu/software/VarScan.v2.3.9/VarScan.v2.3.9.jar"
        
        if pairID == None:
            print("varscan2仅适用于配对分析，请指定配对样本并重新运行")
            exit()

        cmd = """
            samtools mpileup -B -f {reference} -q 15 -d 10000 \\
                {resultsDir}/bam/{pairID}.bam {resultsDir}/bam/{sampleID}.bam \\
                | java -jar {vs} somatic -mpileup {tmpDir}/{sampleID} \\
                --min-coverage-normal {minCov} --min-coverage-tumor {minCov} \\
                --min-var-freq {MAF} --strand-filter 1 --output-vcf
            bcftools reheader -f {reference}.fai {tmpDir}/{sampleID}.indel.vcf -o {tmpDir}/{sampleID}.indel.fix.vcf
            bcftools reheader -f {reference}.fai {tmpDir}/{sampleID}.snp.vcf -o {tmpDir}/{sampleID}.snp.fix.vcf
        """.format(reference=reference, resultsDir=resultsDir, pairID=pairID, sampleID=sampleID, vs=vs, tmpDir=tmpDir, minCov=minCov, MAF=MAF)
        print(cmd)
        os.system(cmd)

    def varscan_filter(self):
        reference = self.reference
        resultsDir = self.output
        sampleID = self.sample
        pairID = self.pair
        filtDP = self.filtDP
        MAF = self.filtQUAL    

        tmpDir = resultsDir + "/tempFile/varscan2_" + sampleID
        outputFile = tmpDir + "/" + sampleID + ".merge.vcf"
        output = open(outputFile, "w")
        indel = open(tmpDir + "/" + sampleID + ".indel.fix.vcf", "r")
        snp = open(tmpDir + "/" + sampleID + ".snp.fix.vcf", "r")

        for i in indel:
            if i.startswith("#"):
                if "NORMAL\tTUMOR" in i:
                    i = i.replace("NORMAL\tTUMOR", sampleID + "\t" + pairID)
            else:
                ii = i.replace("\n", "").split("\t")
                li = ii[0:9]
                li.append(ii[10])
                li.append(ii[9])
                i = "\t".join(li) + "\n"
            output.write(i)
        indel.close()

        for s in snp:
            if not s.startswith("#"):
                ss = s.replace("\n", "").split("\t")
                si = ss[0:9]
                si.append(ss[10])
                si.append(ss[9])
                s = "\t".join(si) + "\n"
                output.write(s)
        snp.close()
        output.close()

        tmp = tmpDir + "/tmp"
        mkdir(tmp)
        cmd = """
            bcftools sort {outputFile} -O v -o {tmpDir}/{sampleID}.vcf -T {tmp}
            cp {tmpDir}/{sampleID}.vcf {resultsDir}/vcf/{sampleID}_{pairID}.vcf
        """.format(resultsDir=resultsDir, outputFile=outputFile, tmpDir=tmpDir, sampleID=sampleID, pairID=pairID, tmp=tmp)
        print(cmd)
        os.system(cmd)


    # pisces tumor only
    # https://github.com/Illumina/Pisces
    def pisces(self):
        """
        需建立索引
        dotnet CreateGenomeSizeFile.dll \
            -g hg19/ \
            -s "Homo sapiens (UCSC hg19)" \
            -o hg19/
        """
        database = self.databases
        resultsDir = self.output
        sampleID = self.sample
        bedFile = self.bed
        threads = self.threads
        minDP = self.filtDP
        minMAF = self.MAF

        piscesBin = "/home/bioinfo/ubuntu/software/Pisces_5.2.10.49/Pisces.dll"

        tmpDir = resultsDir + "/tempFile/pisces_" + sampleID
        mkdir(tmpDir)

        if bedFile != None:
            cmd = """
                dotnet {piscesBin} -b {resultsDir}/bam/{sampleID}.bam \\
                    -g {database} \\
                    -o {tmpDir} \\
                    -t {threads} \\
                    -i {bedFile} \\
                    --mindp {minDP} \\
                    --minvf {minMAF} \\
                    --minvq 0 --threadbychr true
            """.format(bedFile=bedFile, minDP=minDP, minMAF=minMAF, piscesBin=piscesBin, resultsDir=resultsDir, sampleID=sampleID, database=database, tmpDir=tmpDir, threads=threads)
        else:
            cmd = """
                dotnet {piscesBin} -b {resultsDir}/bam/{sampleID}.bam \\
                    -g {database} \\
                    -o {tmpDir} \\
                    -t {threads} \\
                    --mindp {minDP} \\
                    --minvf {minMAF} \\
                    --minvq 0 --threadbychr true
            """.format(minDP=minDP, minMAF=minMAF, piscesBin=piscesBin, resultsDir=resultsDir, sampleID=sampleID, database=database, tmpDir=tmpDir, threads=threads)
        print(cmd)
        os.system(cmd)

        filt = """
            bcftools view \\
                -e "GT='0/0' | GT='./.' | GT='0/.'" \\
                {tmpDir}/{sampleID}.genome.vcf > {tmpDir}/{sampleID}.muts.vcf
            bcftools view \\
                -e "FILTER='LowDP'" \\
                {tmpDir}/{sampleID}.muts.vcf > {tmpDir}/{sampleID}.pisces.vcf
            cp {tmpDir}/{sampleID}.pisces.vcf {resultsDir}/vcf/{sampleID}.vcf
        """.format(tmpDir=tmpDir, sampleID=sampleID, resultsDir=resultsDir)
        print(filt)
        os.system(filt)

    # bcftools
    # http://samtools.github.io/bcftools/bcftools.html
    # 未测试
    def bcftools(self):
        reference = self.reference
        resultsDir = self.output
        sampleID = self.sample
        bedFile = self.bed
        threads = self.threads

        tmpDir = resultsDir + "/tempFile/bcftools_" + sampleID
        mkdir(tmpDir)

        if bedFile != None:
            cmd = """
                bcftools mpileup -f {reference} \\
                    {resultsDir}/bam/{sampleID}.bam \\
                    | bcftools call -mv -O v \\
                    -o {tmpDir}/{sampleID}.bcftools.vcf \\
                    -t {threads} -R {bedFile}
            """.format(bedFile=bedFile, reference=reference, resultsDir=resultsDir, sampleID=sampleID, tmpDir=tmpDir)
        else:
            cmd = """
                bcftools mpileup -f {reference} \\
                    {resultsDir}/bam/{sampleID}.bam \\
                    | bcftools call -mv -O v \\
                    -o {tmpDir}/{sampleID}.bcftools.vcf \\
                    -t {threads}
            """.format(bedFile=bedFile, reference=reference, resultsDir=resultsDir, sampleID=sampleID, tmpDir=tmpDir)
        print(cmd)
        os.system(cmd)


    # freebayes
    # https://github.com/freebayes/freebayes
    # 为了避免后续annovar注释的bug，此处设定genotyping-max-banddepth
    # 同一位点只输出最多突变条数的突变方向
    def freebayes(self):
        reference = self.reference
        resultsDir = self.output
        sampleID = self.sample
        pairID = self.pair
        bedFile = self.bed
        MAF = self.MAF
        DP = self.filtDP

        checkBQSR = sampleID + ".BQSR.bam"
        if checkBQSR in os.listdir(resultsDir + "/bam"):
            bamFile = checkBQSR
        else:
            bamFile = sampleID + ".bam"

        if pairID != None:
            checkBQSR = pairID + ".BQSR.bam"
            if checkBQSR in os.listdir(resultsDir + "/bam"):
                pairFile = checkBQSR
            else:
                pairFile = pairID + ".bam"

        tmpDir = resultsDir + "/tempFile/freebayes_" + sampleID
        mkdir(tmpDir)

        if bedFile != None:
            if pairID == None:
                cmd = """
                    freebayes -f {reference} \\
                        {resultsDir}/bam/{bamFile} \\
                        -t {bedFile} \\
                        -F {MAF} -C 5 --min-coverage {DP} \\
                        --genotyping-max-banddepth 1 \\
                        > {tmpDir}/{sampleID}.freebayes.vcf
                    sed -i "s/0\\/0/0\\/1/g" {tmpDir}/{sampleID}.freebayes.vcf
                    cp {tmpDir}/{sampleID}.freebayes.vcf {resultsDir}/vcf/{sampleID}.vcf
                """.format(DP=DP, MAF=MAF, tmpDir=tmpDir, bedFile=bedFile, bamFile=bamFile, reference=reference, resultsDir=resultsDir, sampleID=sampleID)
            else:
                cmd = """
                    freebayes -f {reference} \\
                        {resultsDir}/bam/{bamFile} \\
                        {resultsDir}/bam/{pairFile} \\
                        -t {bedFile} \\
                        -F {MAF} -C 5 --min-coverage {DP} \\
                        --genotyping-max-banddepth 1 \\
                        > {tmpDir}/{sampleID}_{pairID}.freebayes.vcf
                    cp {tmpDir}/{sampleID}_{pairID}.freebayes.vcf {resultsDir}/vcf/{sampleID}_{pairID}.vcf
                """.format(DP=DP, MAF=MAF, tmpDir=tmpDir, bedFile=bedFile, bamFile=bamFile, pairFile=pairFile, reference=reference, resultsDir=resultsDir, sampleID=sampleID, pairID=pairID)                

        else:
            if pairID == None:
                cmd = """
                    freebayes -f {reference} \\
                        {resultsDir}/bam/{bamFile} \\
                        -F {MAF} -C 5 --min-coverage {DP} \\
                        --genotyping-max-banddepth 1 \\
                        > {tmpDir}/{sampleID}.freebayes.vcf
                    sed -i "s/0\\/0/0\\/1/g" {tmpDir}/{sampleID}.freebayes.vcf
                    cp {tmpDir}/{sampleID}.freebayes.vcf {resultsDir}/vcf/{sampleID}.vcf
                """.format(DP=DP, MAF=MAF, tmpDir=tmpDir, bedFile=bedFile, bamFile=bamFile, reference=reference, resultsDir=resultsDir, sampleID=sampleID)
            else:
                cmd = """
                    freebayes -f {reference} \\
                        {resultsDir}/bam/{bamFile} \\
                        {resultsDir}/bam/{pairFile} \\
                        -F {MAF} -C 5 --min-coverage {DP} \\
                        --genotyping-max-banddepth 1 \\
                        > {tmpDir}/{sampleID}_{pairID}.freebayes.vcf
                    cp {tmpDir}/{sampleID}_{pairID}.freebayes.vcf {resultsDir}/vcf/{sampleID}_{pairID}.vcf
                """.format(DP=DP, MAF=MAF, tmpDir=tmpDir, bedFile=bedFile, bamFile=bamFile, pairFile=pairFile, reference=reference, resultsDir=resultsDir, sampleID=sampleID, pairID=pairID)
        print(cmd)
        os.system(cmd)

    # GATK4 mutect2 过滤
    def gatk_filter(self):
        small_exac = self.databases + "/" + self.runningInfo["setting"]["Mutation"]["gatk_filter"]["small_exac"]
        bedFile = self.bed
        reference = self.reference
        resultsDir = self.output
        sampleID = self.sample

        checkBQSR = sampleID + ".BQSR.bam"
        if checkBQSR in os.listdir(resultsDir + "/bam"):
            bamFile = checkBQSR
        else:
            bamFile = sampleID + ".bam"       

        tmpDir = resultsDir + "/tempFile/gatk_" + sampleID
        mkdir(tmpDir)

        cmd = """
            gatk GetPileupSummaries \\
                -I {resultsDir}/bam/{bamFile} \\
                -O {tmpDir}/{sampleID}.pileups.table \\
                -V {small_exac} \\
                -L {bedFile} \\
                -R {reference}

            gatk CalculateContamination \\
                -I {tmpDir}/{sampleID}.pileups.table \\
                -O {tmpDir}/{sampleID}.contamination.table

            gatk FilterMutectCalls \\
                -R {reference} \\
                -V {tmpDir}/{sampleID}.m2.vcf \\
                -O {tmpDir}/{sampleID}.m2.contFiltered.vcf \\
                --contamination-table {tmpDir}/{sampleID}.contamination.table

            bcftools view \\
                {tmpDir}/{sampleID}.m2.contFiltered.vcf \\
                -f PASS,clustered_events,slippage \\
                > {tmpDir}/{sampleID}.filter.vcf
            cp {tmpDir}/{sampleID}.filter.vcf {resultsDir}/vcf/{sampleID}.vcf
        """.format(tmpDir=tmpDir, bamFile=bamFile, resultsDir=resultsDir, sampleID=sampleID, small_exac=small_exac, bedFile=bedFile, reference=reference)
        print(cmd)
        os.system(cmd)

    # 过滤
    # bcftools
    # http://samtools.github.io/bcftools/bcftools.html
    def filter(self):
        resultsDir = self.output
        sampleID = self.sample
        pairID = self.pair
        DP = self.filtDP
        QUAL = self.filtQUAL

        if pairID == None:
            cmd = """
                bcftools view \\
                    -i 'MIN(FORMAT/DP)>={DP}' \\
                    {resultsDir}/vcf/{sampleID}.vcf \\
                    > {resultsDir}/vcf/{sampleID}.filter.vcf
            """.format(QUAL=QUAL, DP=DP, resultsDir=resultsDir, sampleID=sampleID)
        else:
            cmd = """
                bcftools view \\
                    -i 'FMT/DP[0]>={DP} && FMT/DP[1]>=100' -s {sampleID} \\
                    {resultsDir}/vcf/{sampleID}_{pairID}.vcf \\
                    > {resultsDir}/vcf/{sampleID}.filter.vcf
            """.format(QUAL=QUAL, DP=DP, resultsDir=resultsDir, sampleID=sampleID, pairID=pairID)
        
        print(cmd)
        os.system(cmd)

# end