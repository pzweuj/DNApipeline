#! /usr/bin/python3
# -*- coding: utf-8 -*-
# bcl2fastq

__Version__ = "0.1"
__Author__ = "pzweuj"
__Date__ = "20210312"


import os
import sys
import argparse

def bcl2fastq(data_dir, output, samplesheet, threads):
	b2f = "/home/bioinfo/ubuntu/software/bcl2fastq2/bin/bcl2fastq"

	cmd = """
		{b2f} -i {data_dir} -o {output} \\
			-r {threads} -p {threads} -w {threads} \\
			--sample-sheet {samplesheet}
	""".format(b2f=b2f, data_dir=data_dir, output=output, threads=threads)
	print(cmd)
	os.system(cmd)

def main(data_dir, output, samplesheet, threads):
	bcl2fastq(data_dir, output, samplesheet, threads)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description="bcl2fastq", prog="bcl2fastq.py",
		usage="python3 bcl2fastq.py [-h] -i <data> -o <output> -s <sample_sheet> -t <threads>",
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser.add_argument("-v", "--version", action="version", version="Version 0.1 20210312")
	parser.add_argument("-i", "--input", type=str, help="bcl data directory")
	parser.add_argument("-o", "--output", type=str, help="fastq data directory")
	parser.add_argument("-s", "--sample", type=str, help="sample sheet csv")
	parser.add_argument("-t", "--threads", type=str, help="running threads")

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(data_dir=args.input, output=args.output, samplesheet=args.sample, threads=args.threads)