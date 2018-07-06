#!/usr/bin/python

import os,sys
import glob
import csv

if len(sys.argv)==4:
  out_dir = sys.argv[1]
  fq_dir = sys.argv[2]
  info_file = sys.argv[3]
else:
  print "Incorrect arguments: enter outout directory"
  sys.exit(0)


with open (info_file, "r") as f:
  X = f.readlines()
todir = outdir + "/Summary_result"
if not os.path.exists(to_dir):  os.makedirs(to_dir)



flu_segments = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]

the_key = {}
for x in X:
  xi = x.strip().split(",")
  the_key[xi[0]] = xi

new_key = []
for i in the_key:
  print i, the_key[i]

  for s in flu_segments:
    annot_dir = todir + "/Consensus_genome/" +  the_key[i][-1] 
    if not os.path.exists(annot_dir):  os.makedirs(annot_dir)
    new_dir = todir + "/Consensus_genome/" + the_key[i][-1] + "/" + s
    if not os.path.exists(new_dir):  os.makedirs(new_dir)
    cp_f = "cp " + out_dir + "/*/" + i + "-" + s + "_consensus.fa " + new_dir
    os.system(cp_f)   
    cp_annot = "cp "  + out_dir + "/*/" + i + "-annotation.csv "  + annot_dir
    os.system(cp_annot)
    snp_dir =  todir + "/Snps/" + the_key[i][-1] 
    if not os.path.exists(snp_dir) : os.makedirs(snp_dir)
    cp_snp = "cp "  + out_dir + "/*/" + i + "-genome-snps.vcf "  + snp_dir
    os.system(cp_snp)
    graphs_dir = todir + "/Images/" + the_key[i][-1] # + "/" + i
    if not os.path.exists(graphs_dir) : os.makedirs(graphs_dir)
    #cp_graph = " cp " + hmdir + "/circos/" + i + "-genome_coverage.png "  + graphs_dir
    cp_qc1 = "cp " + fq_dir + "/"  + the_key[i][2] + "-qc-t15-l35_fastqc/Images/per_base_quality.png " +  graphs_dir + "/" + i + "-per_base_quality.png"
    cp_qc2 = "cp " + fq_dir + "/"  + the_key[i][2] + "-qc-t15-l35_fastqc/Images/sequence_length_distribution.png " +  graphs_dir + "/" + i + "-sequence_length_distribution.png"
  #  cp_qc3 = "cp " + hmdir + "/megan_pdf/"  + the_key[i][2] + "-qc-t15-l35_blast_nAll.pdf " +  graphs_dir + "/" + i + "-megan.pdf"
    #os.system(cp_graph)
    os.system(cp_qc1)
    os.system(cp_qc2)
#    os.system(cp_qc3)

# get the mapped summary stats...

  samtool = "samtools flagstat " +  out_dir + "/*/" + i + "-genome-sort.bam  > tmp.txt"
  os.system(samtool)  
  with open ("tmp.txt", "r") as f:
    t=f.readlines()
  tt = []
  for l in t:
    tt.append(l)
  print tt 
  the_key[i].append(tt[0].split(" ")[0])
  the_key[i].append(tt[2].split(" ")[0])
  the_key[i].append((float(the_key[i][-1]) / float(the_key[i][-2])) * 100)

new_key = []
for j in the_key:
  new_key.append(the_key[j])

new_key.insert(0, ["Sample ID", "Designation", "IonXpress barcode", "Subtype" , "Total reads", "Mapped reads", "% of Mapped reads"] )

with open (info_file.replace(".csv", "-res.csv"), "w") as f:
  writer = csv.writer(f)
  writer.writerows(new_key)
    




