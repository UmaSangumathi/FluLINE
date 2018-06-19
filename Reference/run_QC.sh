#!/bin/bash
##############################################################
#  Script     : run_QC.sh
#  Author     : uma sangumathi
#  Date       : 24/05/2015
#  Last Edited: 24/05/2015, uma 
#  Description: Quality check of Ion torrent reads and quality trimming
##############################################################
# Purpose: To automate the QC of the orginal fastq files and trim based on quality score of 15 and remove reads below 35bp length 
#     	Fastqc software can be used to do the quality check and it reports a group of statistics 
#	Fastx-Toolkit comprises of a collection of command line tools for short-reads fasta/fastq files preprocessing.http://hannonlab.cshl.edu/fastx_toolkit/commandline.html
#	
# Requirements:
#	1. FastQC package
#	2. fastx-toolkit
#############################################################

echo "Enter the following paramenters:"
echo "./run_QC.sh  'fastq files directory'  'file name to replace' "

if [ $# -ne 2 ]
  then
    echo "Error: No arguments supplied"
    echo "Please enter as : ./run_QC.sh  'fastq files directory'  'Part of file name to change'"
else 
  cd $1
  name_delete=$2
  change_name="-Qc.fastq"
  for x in *.fastq
  do
    echo "Processing.... " 
    echo $x
    fastqc $x
    #fastq_quality_trimmer -v -Q 33  -t 15 -l 35  -i $x -o $(echo $x | sed "s/$name_delete/$change_name/g")  
    trim_galore -q 20 --length 50 --phred33 -fastqc -o $1 $x
    mv $(echo $x | sed "s/.fastq/_trimmed.fq/g")   $(echo $x | sed "s/$name_delete/$change_name/g")
    #fastqc $(echo $x | sed "s/$name_delete/$change_name/g")
  done
fi

