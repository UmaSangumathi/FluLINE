#!/usr/bin/python
##############################################################
#  Script     : createGraphfiles.py
#  Author     : uma sangumathi
#  Date       : 28/04/2015
#  Last Edited: 24/05/2015, uma 
#  Description: Reformat the files to generate graphs 
##############################################################
# Purpose: 
#  Specific for Influenza virus with 8 segments... 
# Requirements:
#  1. Biopython module    
#  2. Samtools 

#############################################################

import os
import sys
import glob
from Bio import SeqIO


# To create coverage graph.... 


if len(sys.argv)==2:
  xd = sys.argv[-1]
  print xd
else:
  print "Incorrect arguments: enter outout directory"
  sys.exit(0)
for x in glob.glob(xd + "/*-sort.bam"):   
  ref_file = x.replace("-sort.bam",".fa") 
  print ref_file

  # Write coverage per position
  sam_cmd = "samtools mpileup -f " + ref_file + " " + x + "  >  " + x.replace(".bam", ".pileup")
  print sam_cmd
  if not os.path.isfile(x.replace(".bam", ".pileup")):os.system(sam_cmd)
  cov_cmd = " awk '{print $1, $2, $2, $4}'  " + x.replace(".bam", ".pileup")  + "  >  " + x.replace(".bam", ".coverage")
  print cov_cmd
  os.system(cov_cmd)

  # find mean coverage
  cmd = "cat " +  x.replace(".bam", ".coverage") + " | awk '{SUM += $4 }  END { print SUM/ NR}' > " + x.replace(".bam", ".covX") 
  os.system(cmd)




