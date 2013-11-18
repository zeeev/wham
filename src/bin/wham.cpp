//
//  main.cpp
//  wham
//
//  Created by Zev Kronenberg on 12/31/12.
//  Copyright (c) 2012 Zev Kronenberg. All rights reserved.
//

#include  "api/api_global.h"
#include  "api/BamReader.h"

#include  "read_pileup.h"
#include  "randomregion.h"
#include  "flag.h"

#include "boost/lexical_cast.hpp"
#include "boost/math/distributions/binomial.hpp"

using namespace std;
using namespace BamTools;
using namespace boost;

//------------------------------------------------------------
//------------------------------------------------------------

struct flagPileupInfo{
  int good;
  int bad;
  int total_reads;
};

//------------------------------------------------------------
//------------------------------------------------------------

void printFlagCounts(std::map<int, int> flags){    
  for(map<int, int>::iterator it = flags.begin() ; it != flags.end(); it++){
    std::cerr << "flag:\t" << it->first << "\t" << "count:\t" << it->second << "\n";
  }
}

//------------------------------------------------------------
//------------------------------------------------------------

flagPileupInfo processFlagCounts(std::map<int, int> flags){
    
  flagPileupInfo c = {0};
    
  for(map<int, int>::iterator it = flags.begin() ; it != flags.end(); it++){
        
    c.total_reads += it->second;
        
    int passFail = 0;     
    flag currentflag;

    currentflag.addFlag(it->first);
    if(currentflag.returnFlag() == 0){
      passFail = 1;
    }
    if(! currentflag.isPaired()){
      passFail = 1;
    }
    if(! currentflag.isPairMapped()){
      passFail = 1;
    }
    if(! currentflag.sameStrand()){
      passFail = 1;
    }
    if(passFail == 0){
      c.good += it->second;
    }
    else{
      c.bad += it->second;
    }
  }
  return c;
}

//------------------------------------------------------------
//------------------------------------------------------------

void prepFileIterator(const std::string filename, const std::string fileindex, const std::string seqid, BamTools::BamReader & reader){

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

  int i ;
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

  if(! reader.SetRegion(i, 0, i, boost::lexical_cast<int>(seq.Length))){
    std::cerr << "FATAL: cannot set region\n";
    exit(EXIT_FAILURE);
  }

  std::cerr << "INFO: " << "scaffold: " << seq.Name << " ID: " << i << " Seq len: " << seq.Length << "\n";

}

//------------------------------------------------------------
//------------------------------------------------------------
void collectFlags(BamTools::BamReader & reader, map <int, int> & flag_count){

  read_pileup zk;
  BamTools::BamAlignment al;

  while(reader.GetNextAlignment(al)){
    flag_count[al.AlignmentFlag]++;
  }

  std::cerr << "flagstat done\n";
}

//------------------------------------------------------------
//------------------------------------------------------------

void runPileup(BamTools::BamReader & reader, flagPileupInfo & global, const std::string seqid){
                                                                                                                
  read_pileup zk;
  BamTools::BamAlignment al;
  
  while(reader.GetNextAlignment(al)){

      zk.proccess_alignment(al);
   
    if(zk.currentStart() > zk.currentPos()){
      int tmp_pos = zk.currentPos();
      std::list<BamAlignment> dat =  zk.pileup();
            
      std::map<int, int> countingflags;
      int treads = 0;
            
      for(std::list<BamAlignment>::iterator it = dat.begin(); it != dat.end(); it++){
  countingflags[it->AlignmentFlag]++;
      }
      flagPileupInfo local =  processFlagCounts(countingflags);
      float b = boost::lexical_cast<float>(local.bad);
      float t = boost::lexical_cast<float>(local.total_reads);
      float p = b / t;
            
      float gp = boost::lexical_cast<float>(global.bad) / boost::lexical_cast<float>(global.total_reads);
                 
      float sign = 1;
                    
      float pp = boost::math::cdf(complement(boost::math::binomial_distribution<float>(t, gp), b));
               
      std::cout << seqid << "\t" << al.Position << "\t" << b << "\t" << t << "\t" << p << "\t" << pp << "\t" << gp << "\n";
    }
  }
  std::cerr << "wham wham done\n";
}

//------------------------------------------------------------
//------------------------------------------------------------

int main(int argc, const char * argv[])
{
  std::string filename  = argv[1];
  std::string fileindex = argv[2];
  std::string seqid     = argv[3];

  BamTools::BamReader reader_a;
  BamTools::BamReader reader_b;

  prepFileIterator(filename, fileindex, seqid, reader_a);
  prepFileIterator(filename, fileindex, seqid, reader_b);

  map<int, int> flag_count;

  collectFlags(reader_a, flag_count);
  printFlagCounts(flag_count);

  flagPileupInfo global = processFlagCounts(flag_count);

  runPileup(reader_b, global, seqid);

}

//------------------------------------------------------------
//------------------------------------------------------------
