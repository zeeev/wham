//
//  main.cpp
//  wham
//
//  Created by Zev Kronenberg on 12/31/12.
//  Copyright (c) 2012 Zev Kronenberg. All rights reserved.
//

#include <iostream>
#include <string>

#include  "api/api_global.h"
#include  "api/BamReader.h"

#include  "read_pileup.h"
#include  "randomregion.h"
#include  "flag.h"

#include "boost/lexical_cast.hpp"
#include "boost/math/distributions/binomial.hpp"

using namespace std;
using namespace BamTools;


//------------------------------------------------------------
//------------------------------------------------------------

void printFlagCounts(map<string, int> flags){

  int total = 0;

  for(map<string , int>::iterator it = flags.begin() ; it != flags.end(); it++){
    total += it->second;
  }

  for(map<string , int>::iterator it = flags.begin() ; it != flags.end(); it++){
    
    float frq = static_cast<float>(it->second) / static_cast<float>(total);

    std::cerr << it->first << "\t"  << frq << "\n";
  }
}

//------------------------------------------------------------
//------------------------------------------------------------

void prepFileIterator(const std::string filename, const std::string fileindex, const std::string seqid, BamTools::BamReader & reader, int rstart, int rend){

  if ( ! reader.Open(filename)){
    std::cerr << "FATAL: cannot open bam index\n";
    exit(EXIT_FAILURE);
  }
  std::cerr << "INFO: " << filename << " opened for reading\n";
  if( ! reader.OpenIndex(fileindex)){
    std::cerr << "FATAL: cannot open bam index\n";
    exit(EXIT_FAILURE);
  }
  std::cerr << "INFO: " << fileindex << " opened for reading\n";

  int i =0 ;
  BamTools::SamHeader           header = reader.GetHeader();
  BamTools::SamSequenceDictionary seqs = header.Sequences;

  if(! seqs.Contains(seqid)){
    std::cerr << "FATAL: cannot find seqid in BAM header" << "\n";
    exit(EXIT_FAILURE);
  }

  BamTools::SamSequence seq;

  BamTools::SamSequenceConstIterator seqIter = seqs.ConstBegin();
  BamTools::SamSequenceConstIterator seqEnd  = seqs.ConstEnd();
  for ( ; seqIter != seqEnd; ++seqIter ) {
    seq = (*seqIter);
    std::string sname = seq.Name;
    if(sname.compare(seqid) == 0){
      break;
    }
    i++;
  }
  
  if(! reader.SetRegion(i, rstart, i, rend)){
    std::cerr << "FATAL: cannot set region\n";
    exit(EXIT_FAILURE);
  }

  std::cerr << "INFO: " << "scaffold: " << seq.Name << " ID: " << i << " Seq len: " << seq.Length << "\n";

}

//------------------------------------------------------------
//------------------------------------------------------------
void collectFlags(BamTools::BamReader & reader, map <string, int> & flag_count){

  read_pileup zk;
  BamTools::BamAlignment al;




  while(reader.GetNextAlignment(al)){

  
    string flag = boost::lexical_cast<string>(al.AlignmentFlag);


   flag_count[flag]++;

    
    
  }

  std::cerr << "flagstat done\n";
}

//------------------------------------------------------------
//------------------------------------------------------------

int main(int argc, const char * argv[])
{
  std::string filename  = argv[1];
  std::string fileindex = argv[2];
  std::string seqid     = argv[3];
  std::string start     = argv[4];
  std::string end       = argv[5];


  int istart = boost::lexical_cast<int>(start);
  int iend   = boost::lexical_cast<int>(end);

  BamTools::BamReader reader_a;

  prepFileIterator(filename, fileindex, seqid, reader_a, istart, iend);

  map< string, int> flag_count;

  collectFlags(reader_a, flag_count);

  printFlagCounts(flag_count);

  
}

//------------------------------------------------------------
//------------------------------------------------------------
