#!/usr/bin/env python
import pysam
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("-i","--ifile", help = "Seq file name with read information as input in sam/bam format ")
parser.add_argument("-f","--fileformat", default= "bam", help = "File format of input file - either sam or bam; default is for bam format (rb)")
parser.add_argument("-c","--column1", default= "reference", help = "Header of first output column wih reference names, default: reference")
parser.add_argument("-d","--column2", default= "read_count", help = "Header of second output column wih No of counts per reference, default : count")
parser.add_argument("-g","--column3", default= "sample", help = "Header of fourth output column wih name of sample library, default : sample")
parser.add_argument("-l", "--library_name", default= "sample name", help = "Name of library(sample) that will be used as fouth column name")
parser.add_argument("-s","--col_delimiter", default= "\t", help = "Delimiter to seperate the columns of the output file, default : tab(\t)")
parser.add_argument("-o", "--ofile", help = "Column wise output file with Tab seperated output filename of your choice")
args = parser.parse_args()

ifiletype="rb" #default for bam format
if args.fileformat == "sam":
    ifiletype="r"

bamfile = pysam.AlignmentFile(args.ifile, ifiletype)
#bamfile = pysam.AlignmentFile(args.ifile, ifiletype, check_sq=False, check_header=False)
#bamfile = pysam.AlignmentFile(args.ifile, "r", check_sq=False)

reference_counts =  {}
for seq in bamfile:
    if seq.is_unmapped:
        continue

    reference = seq.reference_name
    query_length= seq.query_length  #as per bam/sam ifile
    if reference not in reference_counts:
        reference_counts[reference] =  [] #FIXME
    reference_counts[reference] +=  [query_length]

with open(args.ofile, 'w') as f:
    delim = args.col_delimiter
    f.write(delim.join([args.column1, args.column2, args.column3]) + "\n")
    for reference,count in reference_counts.items():
        f.write(reference + delim + str(len(count)) + delim + args.library_name + "\n")
