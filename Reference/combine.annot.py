#!/usr/bin/python

import os, sys
import csv
import glob

if len(sys.argv)==3:
  res_dir = sys.argv[1]
  info = sys.argv[2]

else:
  print "Incorrect arguments: enter outout directory"
  sys.exit(0)

# Read the summary info file:
info_list = []
with open(info, 'r') as csvfile:
  reader = csv.reader(csvfile)
  for xi in reader:
    print xi
    info_list = xi 
print info_list
# if one samlple or many samples : fixing the list length issue
if len(info_list[0]) < 4 : subtypes = list[set([c[-1] for c in info_list])]
else: subtypes = [info_list[-1],]

# Merge all Annotation file of the consensus genome
all_annot = []
assembled_cons = [["Sample Id", "Sample Name", "Genome"]]

for sub_type in subtypes:
  for x in glob.glob(res_dir + "/Consensus_genome/*csv"):
    X = x.split("/")
    y = X[-1].replace("-annotation.csv", "")
    with open(x, 'rb') as csvfile:
      r = csv.reader(csvfile)
      genome = "-"

      for a in r:
        if a[0] != "Genome":
	  print X, a
          seg_nam = a[0].split("|")[1]
      #    a.insert(0,y + "." + seg[seg_nam]) 
          all_annot.append(a) 
	  if a[1].split("|")[1] == "Genome": genome = a[-1]      
        else: annot_header = a  
    assembled_cons.append([y, a[1].split("|")[0], genome])      	

  with open(res_dir + '/' + sub_type + '-ConsensusDetail.csv', 'wb') as f:
    writer = csv.writer(f)
    annot_header.insert(0,"Sample Id")
    all_annot.insert(0,annot_header)
    writer.writerows(all_annot) 
         
  with open(res_dir + '/' + sub_type + '-ConsensusSummary.csv', 'wb') as f:
    writer = csv.writer(f)
    writer.writerows(assembled_cons) 


# Merge all SNPs called...
merge_snps = []
for sub_type in subtypes:
  for x in glob.glob(res_dir + "/Snps/" + sub_type + "/*.vcf"):
    X = x.split("/")
    y = X[-1].replace("-genome-snps.vcf", "")
    with open(x, 'rb') as csvfile:
      r = csv.reader(csvfile, delimiter="\t")
      for s in r:
	if not s[0].startswith("#"):
	  print s
          seg_nam = s[0].split("|")[1]
	  merge_snps.append(s)

  with open(res_dir + '/' + sub_type + '-SNPs.csv', 'wb') as f:
    writer = csv.writer(f)
    merge_snps.insert(0, ["Sample Id", "Sample Name", "POS","ID","REF","ALT", "QUAL", "FILTER", "INFO"])
    writer.writerows(merge_snps) 






