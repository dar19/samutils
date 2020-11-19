#!/usr/bin/env python
import pysam
import argparse
import re

parser = argparse.ArgumentParser(description='Count reads aligned on reference sequence in a SAM/BAM file')
parser.add_argument("-i", "--ifile", help="Input SAM/BAM file")
parser.add_argument("-a", "--all-referenceids", action='store_true', help="Option for all reference-ids incld. zero count")
parser.add_argument("-s", "--sam", action='store_true', help="Input file format only if SAM file; default BAM format")
parser.add_argument("-r", "--ref-col-name", default="reference", help="Name of output column with reference ids, default: reference")
parser.add_argument("-c", "--cnt-col-name", default="count", help="Name of output column with read count, default: count")
parser.add_argument("-n", "--opt-col-name", help="Name of an optional column e.g. sample_name")
parser.add_argument("-v", "--opt-col-val", help="Value for the optional column; same for all rows")
parser.add_argument("-d", "--delim", default="\t", help="Delimiter to separate the columns of the output file, default : TAB")
args = parser.parse_args()

#Check for the file type SAM/BAM
ifiletype = "rb"
if args.sam:
    ifiletype = "r"

#Create the dictionary for the reference ids
reference_counts = {}

#Process the reads and add counts to respective reference-ids
bamfile = pysam.AlignmentFile(args.ifile, ifiletype)
for seq in bamfile:
    if seq.is_unmapped:
        continue
    reference = seq.reference_name

    if reference not in reference_counts:
        reference_counts[reference] = 0
    reference_counts[reference] += 1

#Condition: Get all reference-ids incld. zero count into dictionary
if args.all_referenceids:
    headers = bamfile.header
    for elem in list(str(headers).split("\n")):
        if "@SQ" in elem:
            reference = re.findall('SN:(.+)\s', elem)[0]
            if reference not in reference_counts:
                reference_counts[reference] = 0

#Delimiter for output file
delim = args.delim

#print the header of file
header_ofile = [args.ref_col_name, args.cnt_col_name]
if args.opt_col_name and args.opt_col_val:
    header_ofile += [args.opt_col_name]
print(delim.join(header_ofile) + "\n")

#print the contents of file
for reference, count in reference_counts.items():
    row = [reference, str(count)]
    if args.opt_col_name and args.opt_col_val:
        row += args.opt_col_val
    print(delim.join(row))
