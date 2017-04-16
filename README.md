# Influenza-Analysis-pipeline
FluLine.py is a wrapper script  

The main steps in the pipeline are 
- i) Filtering of the sequencing reads by cutadapt and FastQC 
-- Quality filter with quality 20 and minimum length 50bp.
-- code = /bin/run_QC.sh

- ii) Find the nearest sequence in NCBI database for each read
 -- Download the NCBI database locally (ftp://ftp.ncbi.nlm.nih.gov/blast/db/). If the nearest reference genome is unknown, then use "/bin/FindSpeciesInSample.py" to generate a XML file with the Blast of all sequence reads against the NCBI database
 -- code = bin//bin/FindSpeciesInSample.py [This is not included in the pipeline, run seperately]
 
- iii) Cluster and identify the viral species, 
 -- MEGAN5 (http://ab.inf.uni-tuebingen.de/software/megan5/) can be used to view the XML file generated by the blast of each read against the NCBI database  

- iv) Generate consensus genomic sequence iteratively
 -- VIPR pipeline is used to iteratively get the consensus. The iterative mapping is done by BWA and the consensus is based on the maximum occurance of the nucleotide at a given position
 -- code = bin/GenerateConsensusGenome_withBlast.py
 
- v) Map the reads to the final consensus genome
  -- Uses Bowtie2 to map reads to the consensus genome
  
- vi) Identify SNVs and 
  -- Uses Lofreq2 to identify the SNVs (genome positions should atleadt have 100 reads mapped)
  
- vii) Visualize the coverage of the genome 
  -- Circos plot is used to visualize the different segments of Influenza
  -- code = bin/createGraphfiles-full.py

