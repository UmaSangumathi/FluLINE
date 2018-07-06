#!/usr/bin/python

import os,sys
import glob
import csv
import fileinput


if len(sys.argv)==4:
  info_file = sys.argv[1]
  plot_dir = sys.argv[2]
  script_dir = sys.argv[3]
else:
  print "Incorrect arguments: enter outout directory"
  sys.exit(0)


with open ( info_file, "r") as f:
  X = f.readlines()

os.chdir(plot_dir)

#flu_segments = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]  
the_key = {}
for x in X:
  xi = x.strip().split(",")
  the_key[xi[0]] = xi
 # with open(xi[0] + "-labels.txt", "w") as f:
  #  f.write("NA" + "\t" + "1" + "\t" +	"1407" + "\t" + str(xi))
for i in the_key:
  with open("tmp.conf", "w") as f:
    for line in fileinput.input("hist.conf", inplace=False):
      f.write(line.replace("00000", i).replace("XXXXXX", "10000")) #.replace("dirpath", plot_dir))
      
  circos_cmd = "circos -conf tmp.conf " 
  os.system(circos_cmd)
#  if not os.path.isfile(hmdir + "/graphs/circos_png/" + i + "-genome_coverage.png"):
#    print i
  
