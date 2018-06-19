#!/usr/bin/python
import os
import sys
import glob
from Bio import SeqIO


# To create coverage graph.... 

graph_order = [("PB2", "hs15"), ("PB1", "hs13"), ("PA", "hs9"), ("HA", "hs8"), ("NP", "hs7"), ("NA","hs17"), ("MP","hs21"), ("NS", "hs1")]
gene_size = {"PB2": 2300, "PB1" : 2300, "PA": 2230, "NP": 1570, "NA":1460, "HA":1760 , "MP": 1030 , "NS":870 }


if len(sys.argv)==2:
  x = sys.argv[-1]

  print x 
  ref_dict = {}
  ref_file = x.replace("-sort.bam",".fa") 
  print ref_file
  for r in SeqIO.parse(ref_file, "fasta"):
    ref_dict[r.id.split("|")[1]] = [r.id, r.id.split("|")[1] , str(1) , str(len(str(r.seq)))]
  ref_genes = ref_dict.keys() 

  # Write coverage per position
  sam_cmd = "samtools mpileup -f " + ref_file + " " + x + "  >  " + x.replace(".bam", ".pileup")
  print sam_cmd
  if not os.path.isfile(x.replace(".bam", ".pileup")):os.system(sam_cmd)
  cov_cmd = " awk '{print $1, $2, $2, $4}'  " + x.replace(".bam", ".pileup")  + "  >  " + x.replace(".bam", ".coverage")
  print cov_cmd
  os.system(cov_cmd)
  for ll in ref_dict:  # Assumin the gene name is 'B/Sydney/13/2014|NS|'
    lll = ref_dict[ll][0].split("|")
    cmdl = "sed -e 's/" + lll[0].replace("/", "\/") + "//g' " + x.replace(".bam", ".coverage")  + " | sed -e 's/|//g'  > " + x.replace(".bam", "-coverage.txt") 
    print cmdl
    os.system(cmdl)    

  # find mean coverage
  for i in ref_genes:
    cmd = "grep '" + i + "' " +  x.replace(".bam", ".coverage") + " | awk '{SUM += $4 }  END { print SUM/ NR}' > tmp.txt"
    os.system("grep '" + i + "' " +  x.replace(".bam", ".coverage") + " | awk '{SUM += $4 }  END { print SUM, NR, SUM/ NR}'")
    print i
    os.system(cmd)
    if os.path.isfile("tmp.txt"):
      for k in  open("tmp.txt", "r"):
        ref_dict[i].append(k.strip("\n") + "x")
      print "Sum: " , k
      os.system("rm tmp.txt")

  for i in gene_size:
    f = 0
    for j in ref_dict:
      if i == j: f = 1
    if f == 0:  
      ref_dict[i] = [i, i, str(1), str(gene_size[i]), "0x"] 
  print ref_dict 


  # Write ideogram flu coords 
  with open( x.replace("-sort.bam","-flu_coords.txt"), "w") as f:
    for j in graph_order:
      for k in ref_dict:
        if j[0] == k:
          f.write("chr - " + ref_dict[k][1] + "\t" + j[0] + "\t" + ref_dict[k][2] + "\t" + ref_dict[k][3] + "\t" + j[1] + "\r\n") 
          

  # Write mean coverage file...
  with open(x.replace("-sort.bam","-covX.txt"), "w") as f:
    for k in ref_dict:
      print ref_dict[k] 
      f.write( ref_dict[k][1] + "\t" +  ref_dict[k][2] + "\t" +  ref_dict[k][3] + "\t" +  ref_dict[k][4] + "\r\n")





