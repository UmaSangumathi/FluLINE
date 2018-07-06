#!/usr/bin/R
library(csv)
library(ggplot2)
library(Gviz)
library(GenomicRanges)
library(GenomicFeatures)

args <- commandArgs(trailingOnly = TRUE)
coord_file = args[1] #"/home/administrator/Desktop/Package_FluAnalysis/RSV/RSV_coords.txt"
depth = args[2] #"./Desktop/Package_FluAnalysis/Run80testanalysis/IonXpress_027/18S0922-genome-sort.coverage"

print(args)
print(coord_file)
print(depth)

coord = read.table(coord_file, sep=",", header=F)
head(coord)
X = read.table(depth, header = F, sep= " ")
head(X)
gr <- GRanges(seqnames = coord$V1, 
              ranges = IRanges(start = coord$V3, width = coord$V4 - coord$V3), id=coord$V2)
dat <- GRanges(seqnames = "chr-1", 
               ranges = IRanges(start = X$V2, width = 1), depth=X$V4)

gtrack <- GenomeAxisTrack(frame=TRUE, col="black", fontcolor="black")
grtrack <- AnnotationTrack(#GeneRegionTrack(
  gr,
#  chromosome = gr, start = st, end = en,
  showId = TRUE, stacking = "dense", # gene=gr@elementMetadata$gene, showFeatureId = TRUE,
  name = "Gene Annotation", fontcolor="black", col.title="black", showOverplotting=TRUE, 
  fill="lightblue",showFeatureId = TRUE,  col.id="black", cex=0.75,frame=TRUE
)


dtrack <- DataTrack(
  range = dat,
  type = "a",
  genome = 'Genome',frame=TRUE,
  name = "Depth of coverage", fontcolor="black", col.title="black", col.axis="black"
)


png(gsub(".coverage", ".png",depth), height=300,width=800)
plotTracks(
  list(dtrack, gtrack, grtrack)
)
dev.off()
