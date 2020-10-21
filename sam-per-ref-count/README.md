"sam-per-ref-count.py" counts reads in sab or bam files.
Below options are used as arguments in the script.

-i IFILE, --ifile IFILE
                        Seq file name with read information as input in
                        sam/bam format
  -f FILEFORMAT, --fileformat FILEFORMAT
                        File format of input file - either sam or bam; default
                        is for bam format (rb)
  -c COLUMN1, --column1 COLUMN1
                        Header of first output column wih reference names,
                        default: reference
  -d COLUMN2, --column2 COLUMN2
                        Header of second output column wih No of counts per
                        reference, default : count
  -g COLUMN3, --column3 COLUMN3
                        Header of fourth output column wih name of sample
                        library, default : sample
  -l LIBRARY_NAME, --library_name LIBRARY_NAME
                        Name of library(sample) that will be used as fouth
                        column name
  -s COL_DELIMITER, --col_delimiter COL_DELIMITER
                        Delimiter to seperate the columns of the output file,
                        default : tab( )
  -o OFILE, --ofile OFILE
                        Column wise output file with Tab seperated output
                        filename of your choice

