/*****************************************************************************
  pairedBamToBed12.cpp

  2011 - Nicolas Bertin
  OSC RIKEN Yokohama
  nbertin@gsc.riken.jp

  directly inspired from 
     bamToBed.cpp
     (c) 2009 - Aaron Quinlan
     Hall Laboratory
     Department of Biochemistry and Molecular Genetics
     University of Virginia
     aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "version.h"
#include "BamReader.h"
#include "BamWriter.h"
#include "BamAncillary.h"
#include "BamAux.h"
#include "bedFile.h"
using namespace BamTools;

#include <vector>
#include <map>
#include <algorithm>    // for swap()
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;


// define our program name
#define PROGRAM_NAME "pairedBamToBed12"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)


// function declarations
void ShowHelp(void);
void ConvertPairedBamToBed12(const string &bamFile, bool showSummary, bool trackUnprocessed, const string &unprocessedBamFile,int minMapQuality, bool delAsBlock, const string &color);
void PrintPairedBed12(const BamAlignment &bam1, const BamAlignment &bam2, const RefVector &refs, bool delAsBlock, string color = "255,0,0");
void ParseCigarBed12(const vector<CigarOp> &cigar, bool delAsBlock, int &currStart, vector<int> &blockStarts, vector<int> &blockEnds, int &alignmentEnd);


int main(int argc, char* argv[]) {

  // our configuration variables
  bool showHelp = false;
  bool showSummary = true;
  // input files
  string bamFile            = "stdin";
  string color              = "255,0,0";
  bool haveBam              = true;
  bool delAsBlock           = false;
  bool trackUnprocessed     = false;
  string unprocessedBamFile = "unprocessedPair.bam";
  int minMapQuality         = 0;


  // do some parsing (all of these parameters require 2 strings)
  for(int i = 1; i < argc; i++) {

    int parameterLength = (int)strlen(argv[i]);
    if((PARAMETER_CHECK("-h", 2, parameterLength)) || (PARAMETER_CHECK("--help", 5, parameterLength))) {
      showHelp = true;
    }
    else if(PARAMETER_CHECK("-quiet", 6, parameterLength)) {
      showSummary = false; 
    }
    else if(PARAMETER_CHECK("-i", 2, parameterLength)) {
      if ((i+1) < argc) {
	bamFile = argv[i + 1];
	i++;
      }
    }
    else if(PARAMETER_CHECK("-color", 6, parameterLength)) {
      if ((i+1) < argc) {
	color = argv[i + 1];
	i++;
      }
    }
    else if(PARAMETER_CHECK("-qual", 5, parameterLength)) {
      if ((i+1) < argc) {
	minMapQuality = atoi(argv[i + 1]);
	i++;
      }
    }
    else if(PARAMETER_CHECK("-dblock", 7, parameterLength)) {
      delAsBlock = true;
    }
    else if(PARAMETER_CHECK("-x", 2, parameterLength)) {
      trackUnprocessed = true;
      if ((i+1) < argc) {
	unprocessedBamFile = argv[i + 1];
	i++;
      }
    }
    else {
      cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
      showHelp = true;
    }
  }
  
  // make sure we have an input files
  if (haveBam == false) {
    cerr << endl << "*****" << endl << "*****ERROR: Need -i (BAM) file. " << endl << "*****" << endl;
    showHelp = true;
  }
  
  // if there are no problems, let's convert BAM to BED12
  if (!showHelp) {
    ConvertPairedBamToBed12(bamFile, showSummary, trackUnprocessed, unprocessedBamFile, minMapQuality, delAsBlock, color);
  }
  else {
    ShowHelp();
  }
}

void ShowHelp(void) {

    cerr << endl << "Program: " << PROGRAM_NAME << " (based on BedTools v" << VERSION << ")" << endl;
    cerr << "Author:  Nicolas Bertin (directly inspired from Aaron Quinlan original bamToBed)" << endl;
    cerr << "Summary: Converts 'properly paired' BAM alignments to BED12 format." << endl;
    cerr << "         Typically producing a 2 blocks BED12 entry for each 'properly paired' BAM pair" << endl;
    cerr << "         Additional blocks are produced when an aligment contains long deletion (CIGAR N-op)" << endl;
    cerr << "         The BAM input file must be grouped/sorted by query name (not alignment position)" << endl << endl;
    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <bam> " << endl << endl;
    cerr << "Options: " << endl;
    cerr << "\t-help\t"    << "Show this help." << endl;
    cerr << "\t-quiet\t"   << "Do not print to sdterr the summary of the number of BAM processed." << endl;
    cerr << "\t-dblock\t"  << "Triggers the creation of a new block when an aligment contains short deletion from reference (CIGAR D-op)" << endl;
    cerr << "\t-color\t"   << "An R,G,B string for the color used with BED12 format." << endl;
    cerr                   << "\t\tDefault is (255,0,0)." << endl;
    cerr << "\t-qual\t"    << "The minimum (inclusive) mapQ sum for reporting the paired BAM into a BED12." << endl;
    cerr                   << "\t\tDefault is (0)." << endl;
    cerr << "\t-x\t"       << "Optional filename where unprocessed mapped pairs can be stored." << endl << endl;

    // end the program here
    exit(1);
}

 
/*
  Assumptions:
     1.  The BAM file is grouped/sorted by query name,
         not alignment position
     2.  Yet  after filtering for quality or requesting a range
      prior to parsing into BED12 there can be some "unpaired-in-name" entries
         Up to 1000 consecutive (aka not interupted by a "paired-in-name" InProperPairBam)
      such entry are tolerated before coming to an halt
*/
void ConvertPairedBamToBed12(const string &bamFile, bool showSummary, bool trackUnprocessed, const string &unprocessedBamFile, int minMapQuality, bool delAsBlock, const string &color) {
  // counter for the number of processed aligments
  int notInProperPairBam = 0;
  int notPairedInNameBam = 0;
  int notHasMinMapQualityBED12 = 0;
  int ProcessedBED12 = 0;
  unsigned int bamBufferMaxSize = 1000;
  // open the BAM file
  BamReader reader;
  reader.Open(bamFile);

  // get header & reference information
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();

  // BamWriter *writer = new BamWriter();
  BamWriter writer;
  if (trackUnprocessed == true){
    writer.Open(unprocessedBamFile, header, refs);
  }

  // bamBuffer is a map <bam.Name, bam>
  map<std::string,BamAlignment> bamBuffer;
  map<std::string,BamAlignment>::iterator bamItr;
  // bamRet allow for tracking "paired-in-name" entries
  pair<map<std::string,BamAlignment>::iterator,bool> bamRet;
  BamAlignment bam;
  while (reader.GetNextAlignment(bam)) {
    if (bam.IsProperPair() == true) {     
      bamRet = bamBuffer.insert(pair<std::string,BamAlignment>(bam.Name,bam));
      if (bamRet.second == false){
	BamAlignment bam1, bam2;
	// processing is ordered 1st (in position) bam then second
	if (bam.Position < bamRet.first->second.Position){ bam1 = bam;                   bam2 = bamRet.first->second; }
	else                                             { bam1 = bamRet.first->second;  bam2 = bam;    }
	// skip pair below the combined MapQ threshold
	if (bam1.MapQuality + bam2.MapQuality >= minMapQuality ){
	  ProcessedBED12++;
	  PrintPairedBed12(bam1, bam2, refs, delAsBlock, color); 
	}
	else{
	  //  pair not satisfying the MapQ threshold can be saved into the unprocessedBamFile
	  notHasMinMapQualityBED12++;
	  if (trackUnprocessed == true){
	    writer.SaveAlignment(bam1);
	    writer.SaveAlignment(bam2);
	  }
	}
	//  Once a "paired-in-name" BAM pair has been seen the content of the buffer should considered irrelevant
	// its content can be saved into the unprocessedBamFile
	if (trackUnprocessed == true){
	  bamBuffer.erase(bam.Name);
	  for (bamItr=bamBuffer.begin(); bamItr != bamBuffer.end(); bamItr++){
	    writer.SaveAlignment( bamItr->second );
	    notPairedInNameBam++;
	  }
	}
	// clear bamBuffer
	bamBuffer.clear();
      }
      if (bamBuffer.size() > bamBufferMaxSize){  
	cerr << PROGRAM_NAME << " *****ERROR: Requires BAM to be sorted/grouped by query name."<< endl;
	cerr << PROGRAM_NAME << " *****ERROR: More than " << bamBufferMaxSize << " consecutive ProperPair BAM with different name were seen before premature halt at " <<  bam.Name.c_str() << endl;
	cerr << PROGRAM_NAME << " *****ERROR: Are you sure the input BAM is properly sorted?" << endl;
	exit(1);
      }
    }
    else{
      notInProperPairBam++;
      // notInProperPairBam can be saved into the unprocessedBamFile
      if (trackUnprocessed == true) writer.SaveAlignment(bam);
    }
  }
  reader.Close();

  if (trackUnprocessed == true) writer.Close();
  if (showSummary){
    int bamTot =  notInProperPairBam +  notPairedInNameBam + (2*notHasMinMapQualityBED12) + (2*ProcessedBED12) ;
    int Bed12Tot = notHasMinMapQualityBED12+ProcessedBED12;
    int notPairedInNamePercent =  int(100 * notPairedInNameBam / bamTot);
    int notInProperPairPercent =  int(100 * notInProperPairBam / bamTot);
    int notHasMinMapQualityPercent = int(100 * notHasMinMapQualityBED12 / Bed12Tot);
    cerr << PROGRAM_NAME << " processed " << bamTot << " BAM alignements producing "<< ProcessedBED12 << " BED12 entries" << endl;
    cerr << "\t" << notInProperPairBam << "(" << notInProperPairPercent << "%) were not considered 'properly paired' BAM alignemts" << endl;
    cerr << "\t" << notPairedInNameBam  << "(" << notPairedInNamePercent << "%) were not 'paired in mame' and were thus skipped" << endl;
    cerr << "\t" << Bed12Tot << " BED12 were produced out of which " << notHasMinMapQualityBED12 << "(" << notHasMinMapQualityPercent << "%) were excluded because their combined MapQ was below " <<  minMapQuality << endl;
  }
}
	




void ParseCigarBed12(const vector<CigarOp> &cigar, bool delAsBlock, unsigned int &currStart, vector<int> &blockStarts, vector<int> &blockLengths, unsigned int &alignmentEnd) {

  int blockLength  = 0;
  //  Rip through the CIGAR ops and figure out if there is more
  //  than one block for this alignment
  vector<CigarOp>::const_iterator cigItr = cigar.begin();
  vector<CigarOp>::const_iterator cigEnd = cigar.end();
  for (; cigItr != cigEnd; ++cigItr) {
    switch (cigItr->Type) {
    case ('M') :
      blockLength  += cigItr->Length;
      currStart += cigItr->Length;
    case ('I') : break;
    case ('S') : break;
    case ('D') : 
      if (delAsBlock == true){
	blockStarts.push_back(currStart + cigItr->Length);
	blockLengths.push_back(blockLength);
	currStart += cigItr->Length;
	blockLength = 0;
      }
      else{
      blockLength  += cigItr->Length;
      currStart += cigItr->Length;
      }
      break;
    case ('P') : break;
    case ('N') :
      blockStarts.push_back(currStart + cigItr->Length);
      blockLengths.push_back(blockLength);
      currStart += cigItr->Length;
      blockLength = 0;
    case ('H') : break;                           // for 'H' - do nothing, move to next op
    default    :
      printf("%s *****ERROR: Invalid Cigar op type\n", PROGRAM_NAME);   // shouldn't get here
      exit(1);
    }
  }
  // add the last block and set the alignment end 
  blockLengths.push_back(blockLength);
  alignmentEnd = currStart;
}
 


void PrintPairedBed12(const BamAlignment &bam1, const BamAlignment &bam2, const RefVector &refs, bool delAsBlock, string color) {

  // set the chrom
  string chrom = refs.at(bam1.RefID).RefName;

  // set the strand
  string strand = "+";
  if (bam1.IsFirstMate()){ if (bam1.IsReverseStrand()) strand = "-"; }
  else                   { if (bam2.IsReverseStrand()) strand = "-"; }

  // set the name of the BED12
  string name = bam1.Name;

  // parse the CIGAR string and figure out the alignment blocks
  unsigned int bam1_alignmentEnd;
  unsigned int bam2_alignmentEnd;
  vector<int> blockLengths;
  vector<int> blockStarts;
  // very first block
  unsigned int currPosition1 = 0;
  blockStarts.push_back(currPosition1);
  // extract the additional block starts and lengths from the CIGAR string of bam1
  ParseCigarBed12(bam1.CigarData, delAsBlock, currPosition1, blockStarts, blockLengths, bam1_alignmentEnd);
  // new block with bam2.position start
  unsigned int currPosition2 =  bam2.Position-bam1.Position;
  blockStarts.push_back(currPosition2);
  // extract the additional block starts and lengths from the CIGAR string of bam2
  ParseCigarBed12(bam2.CigarData, delAsBlock, currPosition2 , blockStarts, blockLengths, bam2_alignmentEnd);

  // set the start and and position of the BED12
  unsigned int alignmentStart;
  unsigned int alignmentEnd;
  alignmentStart = bam1.Position;
  alignmentEnd = alignmentStart + bam2_alignmentEnd;

  // write BED6 portion
  // the score is the sum of the MapQ
  printf("%s\t%d\t%d\t\%s\t%d\t%s\t",
	 chrom.c_str(),
	 alignmentStart,
	 alignmentEnd,
	 name.c_str(),
	 bam1.MapQuality + bam2.MapQuality,
	 strand.c_str()
	 );
  // write the remaining BED12 fields
  // write the txStart / txEnd mark the extend of the 5'read block(s)
  if (bam1.IsFirstMate()){ printf("%d\t%d\t",  bam1.Position, bam1.Position+bam1_alignmentEnd ); }
  else                   { printf("%d\t%d\t",  bam2.Position,  alignmentEnd);  }
  printf("%s\t%d\t",
	 color.c_str(),
	 (int) blockStarts.size()
	 );

    // write the comma delimited blockSizes
    unsigned int b;
    for (b = 0; b < blockLengths.size() - 1; ++b) {
        printf("%d,", blockLengths[b]);
    }
    printf("%d\t", blockLengths[b]);

    // write the comma delimited blockStarts
    for (b = 0; b < blockStarts.size() - 1; ++b) {
        printf("%d,", blockStarts[b]);
    }
    printf("%d\n", blockStarts[b]);
}


