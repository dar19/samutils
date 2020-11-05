#!/usr/bin/env python
import pysam
import argparse

parser = argparse.ArgumentParser(description='Count reads from sam/bam file')
parser.add_argument("-i", "--ifile", help="Seq file name with read information as input in sam/bam format ")
parser.add_argument("-f", "--fileformat", default="bam", help="File format of input file - either sam or bam; default is for bam format (rb)")
parser.add_argument("-c", "--ref-col-name", default="reference", help="Header of column wih references, default: reference")
parser.add_argument("-d", "--cnt-col-name", default="count", help="Header of column with read count per reference, default")
parser.add_argument("-g", "--opt-col-name", help="Header of optional column for sample/library name")
parser.add_argument("-l", "--opt-col-val", help="Value of optional column for sample/library name")
parser.add_argument("-s", "--delim", default="\t", help="Delimiter to seperate the columns of the output file, default : TAB")
parser.add_argument("-o", "--ofile", help="Column wise output file with Tab seperated output filename of your choice")
args = parser.parse_args()

ifiletype = "rb" #default for bam format
if args.fileformat == "sam":
    ifiletype = "r"

bamfile = pysam.AlignmentFile(args.ifile, ifiletype)
#bamfile = pysam.AlignmentFile(args.ifile, ifiletype, check_sq=False, check_header=False)
#bamfile = pysam.AlignmentFile(args.ifile, "r", check_sq=False)

reference_counts = {}
for seq in bamfile:
    if seq.is_unmapped:
        continue

    reference = seq.reference_name
    query_length = seq.query_length  #as per bam/sam ifile
    if reference not in reference_counts:
        reference_counts[reference] = 0
    reference_counts[reference] += 1

delim = args.delim
if args.opt_col_name and args.opt_col_val:
    print(delim.join([args.ref_col_name, args.cnt_col_name, args.opt_col_name]) + "\n")
else:
    print(delim.join([args.ref_col_name, args.cnt_col_name]) + "\n")

for reference, count in reference_counts.items():
    if args.opt_col_name and args.opt_col_val:
        print(reference + delim + str(count) + delim + args.opt_col_val + "\n")
    else:
        print(reference + delim + str(count) + "\n")
