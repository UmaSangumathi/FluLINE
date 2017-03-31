#!/usr/bin/python
##############################################################
#  Script     : GenerateConsensusGenome.py
#  Author     : uma sangumathi
#  Date       : 28/04/2015
#  Last Edited: 24/05/2015, uma 
#  Description: Consensus genome and snps 
##############################################################
# Purpose: 
#  Create guided consensus genome from the reads and nearest reference genome     
#  Map reads using bowtie2    
#  Snps detected using Lofreq2     
# Requirements:
#  1. Biopython module    
#  2. Bwa
#  3. Bowtie2
#  4. components of Vipr pipeline  
#  5. Samtools 
#  6. NCBI blast commandline tool and database
#############################################################

import os
import sys
import glob
import csv
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import Entrez
from optparse import OptionParser


def cmdline_parser():
    """
    creates an OptionParser instance
    """

    # http://docs.python.org/library/optparse.html
    usage = "usage: %prog [options]\nType %prog -h to get more help"
    parser = OptionParser(usage=usage)
    parser.add_option("-r", "--refdir",
                      dest="ref_dir",
                      help="Directory of Reference fasta file")
    parser.add_option("-q", "--fastq_dir",
                      dest="fq_dir",
                      help="Directory of Fastq files")
    parser.add_option("-f", "--fastq_file",
                      dest="f1",
                      help="Fastq file name eg) IonExpress001-Qc.fastq")
    parser.add_option("-b", "--blastn",
                      dest="blast_n",
                      help="Path of Blastn NCBI commanline tool")
    parser.add_option("-d", "--database",
                      dest="db",
                      help="NCBI nt database")
    parser.add_option("-o", "--outdir",
                      dest="out_dir",
                      help="Output directory")
    parser.add_option("-i", "--info",
                      dest="csv_file",
                      help="csv file containing the sample and reference information")
    parser.add_option("-e", "--email",
                      dest="Entrez_email",
                      help="Email ID is required to access the NCBI website online")
    return parser

def main():
  """The main function
  """
  try:
    parser = cmdline_parser()
    (opts, args) = parser.parse_args()
    print opts,args
    if len(args):
      parser.error("Unrecognized arguments found: %s." %( 
		' '.join(args)))
      sys.exit(1)
      
  except:
    parser.print_help()
    sys.exit(0)
  Entrez.email = opts.Entrez_email
  threads = str(3)
  tab_fmt = "'6 qseqid sseqid qstart qend sstart send pident qcovs evalue bitscore stitle'"
  print opts.f1, args
# Read the summary file:
  info_dict = {}
  with open(opts.csv_file, 'r') as csvfile:
    reader = csv.reader(csvfile)
    for xi in reader:
      print xi
      info_dict[xi[-2]] = xi 
  print info_dict.keys()

  os.chdir(opts.fq_dir)
  ref = opts.ref_dir + "/" + info_dict[opts.f1.replace("-Qc.fastq", "")][-1] + ".fa" 
  if not os.path.exists(opts.out_dir): os.mkdir(opts.out_dir)
  res_dir = opts.out_dir + "/" + opts.f1.replace("-Qc.fastq", "")
  if not os.path.exists(res_dir): os.mkdir(res_dir)
  flag = 0
  for x in SeqIO.parse(open(ref, "r"), "fasta"):
    split_ref = res_dir + "/" + x.id.replace("_", "").replace("|","").replace("/","") + ".ref" 
    with open(split_ref, "w") as ff:
      ff.write(">" + x.id + "\r\n" + str(x.seq) + "\r\n") 
    
    # Map reads to a reference genome...
    cons_ref = res_dir + "/" + opts.f1.replace("-Qc.fastq", x.id.replace("_", "").replace("|","").replace("/","") +".fa")    
    cons_cmd = "bam2cons_iter-lenient.sh  -f " + opts.f1 + " -r " + split_ref + " -t " + threads + " --force -o " + cons_ref    
    print cons_cmd
    if not os.path.isfile(cons_ref): os.system(cons_cmd) 
     
    cons = SeqIO.parse(cons_ref, "fasta")
    nee = [xl for xl in cons]
    if len(nee) !=0: 
      X = []
      for xl in nee:
        nnn = str(xl.seq).replace("N","")
        if len(nnn) > 30 : X.append(nnn)
    else: 
      X=[]     
    print "NEE", nee, X
    if len(nee) != 0 and len(X) != 0:
      # Blast the consensus formed.... 
      b_out = res_dir + "/" + opts.f1.replace(".fastq", ".tmp")
      blast_run_uni = NcbiblastnCommandline(cmd=opts.blast_n, task="megablast", db=opts.db, max_target_seqs=1, query=cons_ref , outfmt=tab_fmt, out=b_out)
      print "Runnin blast...", blast_run_uni
      nearest_ref1 = cons_ref.replace(".fa", "-nref1.fa")
      if not os.path.isfile(nearest_ref1): 
        stdout, stderr = blast_run_uni()
        c = 0
        for line in open(b_out, "r"):
          c = c+1
          if c ==1 : u_id=line.split("|")
        print u_id, u_id[3] 
        handle = Entrez.efetch(db="nuccore", id=u_id[3], rettype="gb", retmode="text",  idtype="acc")  # Use the PrimaryID instead of GI
 	#handle = Entrez.efetch(db="nucleotide", id="AY851612", rettype="gb", retmode="xml")
	#handle = "elink -db nuccore -query " + u_id[3] + " |efetch -format gb "
 	#print handle
	record = SeqIO.read(handle, 'genbank')
        #record = Entrez.read(handle) #, validate=False)
	#print record
        with open( nearest_ref1, "w") as f :
          #f.write(">" + record[0]["GBSeq_primary-accession"] + "\r\n" + record[0]["GBSeq_sequence"] + "\r\n")
   	  f.write(">" + record.id + "\r\n" + str(record.seq) + "\r\n")	
    # 2-iteration Map reads to nearest reference...
    # IonXpress_012-Qc.fastq
      ddd = opts.f1.replace("-Qc.fastq","")
      cons_ref2 = res_dir + "/" + info_dict[ddd][0] + "-"  + x.id + "cons2.fa" 
      cons_cmd2 = "bam2cons_iter-lenient.sh -f " + opts.f1 + " -r " + nearest_ref1 + " -t " + threads + " --force -o " + cons_ref2    
      print cons_cmd2
      if not os.path.isfile(cons_ref2):os.system(cons_cmd2)


    # combine all seqments into one ref...
      flag = 1
      Cons_ref_ful = res_dir + "/" + info_dict[ddd][0] + "-genome.fa"
      for i in SeqIO.parse(open(cons_ref2, "r"), "fasta"):
        with open(Cons_ref_ful, "a") as f :   
	  print i
          f.write(">" + info_dict[ddd][1].replace(" ","") + "|" +  i.id.split("-")[1].replace("cons2","") + "|" + "\r\n" + str(i.seq) + "\r\n")
  else:
    print "No Coverage :( "

  # Map reads to consensus genome.... Bowtie2
  if flag == 1:
    cons_ref2_indx = "bowtie2-build " + Cons_ref_ful  + " tmp_idx" 
    os.system(cons_ref2_indx)
    out_sam = Cons_ref_ful.replace(".fa", ".sam")
    bowtie2_cmd = "bowtie2 --local --fast-local  -x  tmp_idx  -U " + opts.f1 + " -S  " + out_sam + " -p 8"
    os.system(bowtie2_cmd)

  # run lofreq2....
    bam_cmd = "samtools view -b -S " + out_sam  + " > " + out_sam.replace(".sam", ".bam")
    bam_sort = "samtools sort " + out_sam.replace(".sam", ".bam") + " " + out_sam.replace(".sam", "-sort")
    cmd = "lofreq   call -C 100 -f  " + Cons_ref_ful  + " -o " + out_sam.replace(".sam", "-snps.vcf") + " " + out_sam.replace(".sam", "-sort.bam")
    # 100 is min depth required to call SNPs
    print bam_cmd
    os.system(bam_cmd)
    print bam_sort
    os.system(bam_sort)
    print cmd
    os.system(cmd)
    print "\n\n"
  else: print "Not this strain..." 

if __name__ == "__main__":
    main()



