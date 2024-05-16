# GeneInfoReader

A simple program to parse GFF files and produce some output. 

Usage: 
```
inputFile=... specify gff3 input file (required)
out=/tmp/     specify an output directory (optional)
-noqc         provide this flag to suppress QC output
-geneTrack    provide this flag to generate a serialized GeneTrack for Genvisis
-genesXln     provide this flag to generate an xln file of genes
-bedExons     provide this flag to generate a bed file of exons
-bedIntrons   provide this flag to generate a bed file of introns
-bedAll       provide this flag to generate three bed files. One of exons, one of introns, and one containing both.
```