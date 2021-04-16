# coding=utf-8
# pzw
# 20210416

# 从Clinvar、LRG、HGNC中获得参考转录本ID
# Clinvar：https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
# LRG：ftp://ftp.ebi.ac.uk/pub/databases/lrgex/list_LRGs_transcripts_xrefs.txt
# HGNC：https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt

import gzip

refTransript = {}
outputFile = open("refTransript.txt", "w")
outputFile.write("#Gene\tTranscript\tSource\n")


HGNC = open("gene_with_protein_product.txt", encoding="utf-8", errors="ignore")
for h in HGNC:
    if not h.startswith("hgnc_id"):
        gene = h.split("\t")[1]
        ts = h.split("\t")[23]
        refTransript[gene] = [ts, "HGNC"]
HGNC.close()

# 次优先
Clinvar = gzip.open("variant_summary.txt.gz", "r")
for c in Clinvar:
    c = c.decode("utf-8")
    if not c.startswith("#"):
        gene_ts = c.split("\t")[2]
        if "NM_" in gene_ts:
            try:
                gene = gene_ts.split("(")[1].split(")")[0]
                ts = gene_ts.split("(")[0]
                refTransript[gene] = [ts, "Clinvar"]
            except:
                continue
Clinvar.close()

# 最优先
LRG = open("list_LRGs_transcripts_xrefs.txt", "r")
for l in LRG:
    if not l.startswith("#"):
        gene = l.split("\t")[1]
        ts = l.split("\t")[4]
        refTransript[gene] = [ts, "LRG"]
LRG.close()

for k in refTransript.keys():
    outputFile.write(k + "\t" + refTransript[k][0] + "\t" + refTransript[k][1] + "\n")
outputFile.close()