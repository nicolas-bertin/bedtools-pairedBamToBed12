=================================================
   (Paired-End Bam to BED12 modified) BEDTools         
=================================================

Addition of a pairedBamToBed12 utility by Nicolas Bertin (OSC RIKEN Yokohama) to BEDTools Created by Aaron Quinlan Spring 2009.

Copyright 2009,2010,2011 Aaron Quinlan. All rights reserved.
Released under GNU public license version 2 (GPL v2).

Repository:      https://github.com/nicolas-bertin/bedtools-pairedBamToBed12


Summary
-------
BEDTools is a collection of utilities for comparing, summarizing, and 
intersecting genomic features in BED, GTF/GFF, VCF and BAM formats. 
In this branch is added *pairedBamToBed12* which converts a name sorted paired-end sequences bam file into a 2-block bed12 formatted file.


Manual
------
See the extensive PDF manual included at: http://code.google.com/p/bedtools/downloads/detail?name=BEDTools-User-Manual.v4.pdf.

This manual covers many common usage examples.  There are also examples available at:
http://code.google.com/p/bedtools/wiki/Usage
http://code.google.com/p/bedtools/wiki/UsageAdvanced
Note that the manual does not reflect the pairedBamToBed12 addition.

Installation
------------
  #. Unpack the source downloaded tarball.
  #. cd into the expanded folder.
  #. Type "make clean" and hit enter.
  #. Type "make all" and hit enter.
  #. If you encountered no errors, then all of the BED Tools should now be in bin/
     If not, try to troubleshoot then email me: aaronquinlan at gmail dot com
  #. Copy the files in bin/ to ~/bin or if you have the privileges, to /usr/local/bin.
  #. Use the tools.


List of available tools
-----------------------

=========================  =======================================================================================================
Utility                    Description
=========================  =======================================================================================================
*intersectBed (BAM)*       Returns overlaps between two BED/GFF/VCF files. 
*pairToBed (BAM)*          Returns overlaps between a paired-end BED file and a regular BED/VCF/GFF file.
*bamToBed (BAM)*           Converts BAM alignments to BED6, BED12, or BEDPE format.
*bedToBam (BAM)*           Converts BED/GFF/VCF features to BAM format.
*bed12ToBed6*              Converts "blocked" BED12 features to discrete BED6 features.
*bedToIgv*                 Creates IGV batch scripts for taking multiple snapshots from BED/GFF/VCF features.
*coverageBed (BAM)*        Summarizes the depth and breadth of coverage of features in one BED versus features (e.g, "windows", exons, etc.) defined in another BED/GFF/VCF file. 
*genomeCoverageBed (BAM)*  Creates either a histogram, BEDGRAPH, or a "per base" report of genome coverage. 
*unionBedGraphs*           Combines multiple BedGraph files into a single file, allowing coverage/other comparisons between them. 
*annotateBed*              Annotates one BED/VCF/GFF file with overlaps from many others. 
*groupBy*                  Summarizes data in a file/stream based on common columns.
*overlap*                  Returns the number of bases pairs of overlap b/w two features on the same line.
*pairToPair *              Returns overlaps between two paired-end BED files. 
*closestBed*               Returns the closest feature to each entry in a BED/GFF/VCF file. 
*subtractBed*              Removes the portion of an interval that is overlapped by another feature. 
*windowBed (BAM)*          Returns overlaps between two BED/VCF/GFF files based on a user-defined window. 
*mergeBed*                 Merges overlapping features into a single feature. 
*complementBed*            Returns all intervals _not_ spanned by the features in a BED/GFF/VCF file. 
*fastaFromBed*             Creates FASTA sequences based on intervals in a BED/GFF/VCF file. 
*maskFastaFromBed*         Masks a FASTA file based on BED coordinates. 
*shuffleBed*               Randomly permutes the locations of a BED file among a genome. 
*slopBed*                  Adjusts each BED entry by a requested number of base pairs. 
*sortBed*                  Sorts a BED file by chrom, then start position. Other ways as well. 
*linksBed*                 Creates an HTML file of links to the UCSC or a custom browser. 
*pairedBamToBed12*         Convert a name sorted paired-end sequences bam file into a 2-block bed12 formatted file.
=========================  =======================================================================================================
