"sam-per-ref-count.py"
Description : Count reads aligned on reference sequence in a SAM/BAM file.
Below options are used as arguments in the script.

usage: sam-per-ref-count.py [-h] [-i IFILE] [-a] [-s] [-r REF_COL_NAME]
                            [-c CNT_COL_NAME] [-n OPT_COL_NAME]
                            [-v OPT_COL_VAL] [-d DELIM]

optional arguments:
  -h, --help            show this help message and exit
  -i IFILE, --ifile IFILE
                        Input SAM/BAM file
  -a, --all-referenceids
                        Option to fetch all reference-ids incld. with zero count
  -s, --sam             Input file format only if SAM file; default BAM format
  -r REF_COL_NAME, --ref-col-name REF_COL_NAME
                        Name of output column with reference ids, default:
                        reference
  -c CNT_COL_NAME, --cnt-col-name CNT_COL_NAME
                        Name of output column with read count, default: count
  -n OPT_COL_NAME, --opt-col-name OPT_COL_NAME
                        Name of an optional column e.g. sample_name
  -v OPT_COL_VAL, --opt-col-val OPT_COL_VAL
                        Value for the optional column; same for all rows
  -d DELIM, --delim DELIM
                        Delimiter to separate the columns of the output file,
                        default : TAB
