#!/usr/bin/python
##############################################################
#  Script     : N_consensus.py
#  Author     : Uma Sangumathi
#  Date       : 28/04/2015
#  Last Edited: 24/05/2015, Uma 
#  Description: 5' and 3' region N's are replaced with the recovered sequences from the softclipped regions
##############################################################
# Purpose: 
#  Specific for the output consensus created from pipeline
# Requirements:
#  1. Biopython module    
#  2. Samtools 
#  3. extractSoftclipped tool
#  4. gunzip
#
# Note: For lenient parameters used in the consensus generation, there is no grt improvement.  
#############################################################


# To edit the consensus genome based on the reads mapped
import os, sys
import glob
import csv
from Bio.Seq import Seq
from Bio import SeqIO
from operator import itemgetter
from itertools import groupby


def write_csv(data, file_nam):
    import csv, itertools
    with open(file_nam, 'wb') as f:
        writer = csv.writer(f, delimiter=",")
        print "Witing to :", file_nam
        writer.writerows(list(itertools.chain(t)) for t in data)

# Get the  all the softclipped sequences 
def get_softclip(bam_file):
  softclipped_seqs = {}
  ext_softclip = "extractSoftclipped"  
  cmd = ext_softclip + " " + bam_file  + " " + " > " + bam_file.replace(".bam", "-softclipped.fq.gz")
  os.system(cmd)
  os.system("gunzip " + bam_file.replace(".bam", "-softclipped.fq.gz"))
  X = SeqIO.parse(bam_file.replace(".bam", "-softclipped.fq"), "fastq")
  for x in X:
    if x.id not in softclipped_seqs: softclipped_seqs[x.id.replace("|:", "|*")] = [str(x.seq)] 
    else: softclipped_seqs[x.id.replace("|:", "|*")].append(str(x.seq)) 
  return (softclipped_seqs)

# From the VIPR pipeline find the N's in the consensus genome
def get_Ncoordinates(consensus_file, Numextension):
  consen = SeqIO.parse(consensus_file, "fasta")
  cons_dict = {}
  for x in consen:
    cons_dict[x.id] = "N" * int(Numextension)  + str(x.seq) + "N" * int(Numextension)
    #print cons_dict[x.id], x.id
  N_pos = {}
  def findOccurences(s, ch):
    return [i for i, letter in enumerate(s) if letter == ch]
  for gene_name in cons_dict:
    N_pos[gene_name] = findOccurences(cons_dict[gene_name], "N")   
  return N_pos, cons_dict

# Groups N positions 
def groupNs(data):
    ranges = []
    for k, g in groupby(enumerate(data), lambda (i,x):i-x):
      group = map(itemgetter(1), g)
      ranges.append((group[0], group[-1]))
    return ranges

# get the consensus from the softclipped sequences
def recover_consensus(bam_file, Npos, cons_dict):
  softclip_seq = get_softclip(bam_file)   # Get the softclipped sequences [starts from number 1]
  #print softclip_seq.keys()
  softclip_seq_Nt = {}
  for g in Npos:
    subs = []   
    Nrange = groupNs(Npos[g]) 
    #print "Nrange", Nrange
    for r in Nrange:
      ntAfterN = [cons_dict[g][r[0]-1:r[0]-6], cons_dict[g][r[1]+1:r[1]+6]]
      end_g = len(cons_dict[g]) -1 
      s = r[0] 
      e = r[1] 
   #   if s < 0: s = 0
   #   if e > end_g: e = end_g  
      if s ==0 : 
 	N_location = "5prime"
	ck = "reverse"
        pt = e + 1
      elif e == end_g : 
	N_location = "3prime"
	ck = "forward"
        pt = s - 1
      elif s != 0 and e != end_g : N_location = "middle"; ck = "middle"
      # process the N's at the 3 amd 5 prime ends ; this code doesnt correct the N's in the middle of the genome
      len_r = r[1] - r[0] +2 
      tmp_lk = {}
      if ck == "forward" : 
        if g + "*" + str(pt) in softclip_seq: 
          #print g + "*" + str(pt), softclip_seq[g + "*" + str(pt)]
          for j in range(0,len(softclip_seq[g + "*" + str(pt)])):
            for k in range(0,len_r):
              if len(softclip_seq[g + "*" + str(pt)][j]) <= k  : 
		if k not in tmp_lk: tmp_lk[k] =["-"]
		else: tmp_lk[k].append("-")
              else:
 		if k not in tmp_lk: tmp_lk[k] = [softclip_seq[g + "*" + str(pt)][j][k]]
                else: tmp_lk[k].append(softclip_seq[g + "*" + str(pt)][j][k]) 
          new_tmp = ["".join(tmp_lk[s]) for s in range(0,len_r)] 
          if g+  "*" + str(pt) not in softclip_seq_Nt: softclip_seq_Nt[g + "*" + str(pt) + "*" + N_location] = new_tmp
          else: softclip_seq_Nt[g + "*" + str(pt) + "*" + N_location].extend(new_tmp)
        else:
          print "missing" , g + "*" + str(pt)
      elif ck == "reverse":
        if g + "*" + str(pt) in softclip_seq:
          #print g + "*" + str(pt), softclip_seq[g + "*" + str(pt)]
	  for j in range(0,len(softclip_seq[g + "*" + str(pt)])):
            for k in range(-1,(len_r+1)*-1,-1):
              if len(softclip_seq[g + "*" + str(pt)][j]) <= k*-1: 
                if k not in tmp_lk: tmp_lk[k] =["-"]
                else: tmp_lk[k].append("-")
              else:
                if k not in tmp_lk: tmp_lk[k] = [softclip_seq[g + "*" + str(pt)][j][k]]
                else: tmp_lk[k].append(softclip_seq[g + "*" + str(pt)][j][k])
         # print tmp_lk, k , tmp_lk[k] 
          new_tmp = ["".join(tmp_lk[s]) for s in range(-1,len_r*-1, -1)] 
          if g+  "*" + str(pt) not in softclip_seq_Nt: softclip_seq_Nt[g + "*" + str(pt) + "*" + N_location] = new_tmp
          else: softclip_seq_Nt[g + "*" + str(pt) + "*" + N_location].extend(new_tmp)
        else:
          print "missing" , g + "*" + str(pt)
      print N_location, Nrange , r, g , pt 
  print softclip_seq_Nt  
  print [(x, len(softclip_seq_Nt[x])) for x in softclip_seq_Nt]     
  return softclip_seq_Nt

       

# Get the high occuring nt and also check the orientation
def get_Ntcomp_new(softclip_seq_Nt):
  recovered_nt = {}
  for x in softclip_seq_Nt:
    nt_seq = []
    for i in range(0,len(softclip_seq_Nt[x])):
      A = softclip_seq_Nt[x][i].count("A")
      G = softclip_seq_Nt[x][i].count("G")
      C = softclip_seq_Nt[x][i].count("C")
      T = softclip_seq_Nt[x][i].count("T")
      tot = A + G + C + T
      ck = [A, G, C, T]
      cons_ind =  ck.index(max(ck))
      if cons_ind == 0: cons = ["a" , str(A)]
      if cons_ind == 1: cons = ["g" , str(G)]
      if cons_ind == 2: cons = ["c" , str(C)]
      if cons_ind == 3: cons = ["t" , str(T)]
      if A == T == G == C == 0: cons = ["N" , str(softclip_seq_Nt[x][i].count("-"))]
      nt_seq.extend(cons[0]) 
    if "5prime" in x: recovered_nt[x] = "".join(nt_seq)[::-1] 	# change the orientation of the 5prime end
    elif "3prime" in x: recovered_nt[x] = "".join(nt_seq)
  print recovered_nt
  return recovered_nt

# replace the N's with recovered sequence
def combine_genome(Npos, cons_dict, softclip_seq_Nt, recovered_nt):
  New_gen = {}
  for s in cons_dict:
    New_gen[s] = list(cons_dict[s])
  for x in recovered_nt:
    [g, pos, loc] = x.split("*")  
    if loc == "5prime": 
      for i in range(0, int(pos)):
        New_gen[g][i] = recovered_nt[x][i]
    elif loc == "3prime":
      for i in range(int(pos), len(New_gen[g])):
        New_gen[g][i] = recovered_nt[x][i - int(pos)]
    with open(bam_file.replace("genome-sort.bam", g.split("|") [-2] + "-new_cons.fa"), "w") as f :
      f.write(">" + g + "\r\n" + "".join(New_gen[g]) + "\r\n")
  return New_gen        

def create_summary_table(New_gen):
  summary = []
  for g in New_gen:
    cons = "F"
    start_pos = 0
    end_pos = len(New_gen[g])
    for i in range(0,len(New_gen[g])):
      if New_gen[g][i] != "N" and cons == "F": 
	start_pos = i
	cons = "T"
        break
    cons = "F"
    for i in range(-1, len(New_gen[g])*-1, -1):
      #print "inloop :", i   
      if New_gen[g][i] != "N" and cons == "F": 
        end_pos = len(New_gen[g]) + i + 1
        cons = "T" 
        break
    indices_N = [i for i, x in enumerate(New_gen[g]) if x == "N"]
    rec_nt_p = [j for j in range(0,len(New_gen[g])) if New_gen[g][j].islower() ]
    rec_nt_r = groupNs(rec_nt_p)   
    rec_nt_s = [ str(rec_nt_r[l][0] + 1) + "-" + str(rec_nt_r[l][1] + 1) for l in range(0,len(rec_nt_r))]
    assem_nt_p = [j for j in range(0,len(New_gen[g])) if New_gen[g][j].isupper() ]
    a_b = list(set(assem_nt_p) - set(indices_N))
    assem_nt_r = groupNs(a_b)
    assem_nt_s = [ str(assem_nt_r[l][0] + 1) + "-" + str(assem_nt_r[l][1] + 1) for l in range(0,len(assem_nt_r))]
    #print rec_nt_r, start_pos, end_pos , New_gen[g][start_pos : end_pos],  assem_nt_r  
    with open(bam_file.replace("genome-sort.bam", g.split("|") [-2]  + "_consensus.fa"), "w") as f :
      f.write(">" + g + "\r\n" + "".join(New_gen[g][start_pos : end_pos]) + "\r\n")
    if start_pos == 0 : clip_5 = "NA"
    else: clip_5 = "1-" + str(start_pos)
    if end_pos == len(New_gen[g]) : clip_3 = "NA"
    else: clip_3 = str(end_pos + 1) + "-" + str(len(New_gen[g]))
    summary.append([g, clip_5 + ";", clip_3 + ";" , "; ".join(rec_nt_s) + ";", "; ".join(assem_nt_s) + ";"])
  summary.insert(0, ["Genome", "5' clipped", "3' clipped", "Recovered portion", "Assembled consensus"])
  writer = csv.writer(open(bam_file.replace("genome-sort.bam", "annotation.csv"), 'w'))
  for m in summary:
    writer.writerow(m)

# The main script
if len(sys.argv)==3:
  ext_num = sys.argv[1]
  xd = sys.argv[-1]
 # print xd
 # print glob.glob(xd + "/*-genome.fa")
  for fil in glob.glob(xd + "/*-genome.fa"):
#for fil in glob.glob("./Result/*/*-genome.fa"):
 # filen = fil.split("-") 
    consensus_file = fil
    bam_file = fil.replace("-genome.fa",  "-genome-sort.bam")
   # print "get_Ncoordinates(consensus_file)"
    [Npos, cons_dict] = get_Ncoordinates(consensus_file, ext_num)
   # print "Npos", Npos, cons_dict.keys()
    softclip_seq_Nt = recover_consensus(bam_file, Npos, cons_dict)
    recovered_nt = get_Ntcomp_new(softclip_seq_Nt)
    New_gen = combine_genome(Npos, cons_dict, softclip_seq_Nt, recovered_nt)
    create_summary_table(New_gen)
else:
  print "Incorrect arguments: enter outout directory"
  sys.exit(0)


