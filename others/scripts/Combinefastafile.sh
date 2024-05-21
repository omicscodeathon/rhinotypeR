#!/bin/bash 

# Usage
# bash Combinefastafile 'pathtoinputfile' 'pathtofile/output.fasta'
<<com
This script helps analysts to combine the sequence fasta files.
Step 1: List down all files in that directory.
Step 2: Extract the name of the file excluding the extension fasta.
Step 3: Extract the sequence header by using cat and grep to capture just the line starting with >.
Step 4: Combine the basefilename and header to creater new header.
Step 5: Extract the sequence.
Step 6: Concatenate the newheader and sequence separated by new line \n.
Step 7: Create an empty file
Step 8: Write into that new file supplied by the user

Arguments:
$1: input filepath
$2: outputfilename
com



for file in $(ls $1*.fasta)
do
    #echo $file # debug to see the files neing listed within the current directory
    name=$(basename $file | cut -f1 -d '.')
    #echo $name
    # Extract the sequence header
    header=$(cat $file | grep '^>' | cut -f1 -d ' ')
    #echo $header
    #sequence header with accession number and RV type
    seqheader=${header}'_'${name}
    #echo $seqheader
    # Extract the sequence
    seq=$( cat $file | grep -v '^>')
    #echo $seq
    #write them to a file and append the rest of the files
    header_seq="$seqheader\n$seq"
    dir_filename=$2
    touch $dir_filename  # create a file name that is sensible to you and reflects the contents of the file
    echo -e $header_seq >> $dir_filename # echo -e enables interpretation of backslash escapes

done