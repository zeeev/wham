//
//  main.cpp
//  wham
//
//  Created by Zev Kronenberg on 12/31/12.
//  Copyright (c) 2012 Zev Kronenberg. All rights reserved.
//

// starting development

//g++ -o read_oddity -I ../bamtools-1.0.2/include/ -L ../bamtools-1.0.2/lib/ -I ../lib/ -I boost -lbamtools wham.cpp ../lib/*cpp

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


void printFlagCounts(std::map<int, int> flags){
    
  for(map<int, int>::iterator it = flags.begin() ; it != flags.end(); it++){
    std::cerr << "flag:\t" << it->first << "\t" << "count:\t" << it->second << "\n";
  }
}

struct flagPileupInfo{
  int good;
  int bad;
  int total_reads;
};


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


int main(int argc, const char * argv[])
{
  std::string filename  = argv[1];
  std::string fileindex = argv[1]; 
  std::string seqid     = argv[2];  
  fileindex += ".bai";

  BamTools::BamReader reader;
  
  if ( ! reader.Open(filename)){
    std::cerr << "FATAL: cannot open bam index\n";     
    return 1;
  }
  std::cerr << "INFO: " << filename << " opened for reading\n";
  if( ! reader.OpenIndex(fileindex)){
    std::cerr << "FATAL: cannot open bam index\n";
    return 1;
  }
  std::cerr << "INFO: " << fileindex << " opened for reading\n";
  BamTools::SamHeader header;
  
  header = reader.GetHeader();

  BamTools::SamSequenceDictionary seqs;
  
  seqs = header.Sequences;
  
  if(! seqs.Contains(seqid)){
    std::cerr << "FATAL: cannot find seqid in SAM header" << "\n";
    return 1; 
  }

   SamSequence seq;
   int i = 0;
   BamTools::SamSequenceConstIterator seqIter = seqs.ConstBegin();
   BamTools::SamSequenceConstIterator seqEnd  = seqs.ConstEnd();
   for ( ; seqIter != seqEnd; ++seqIter ) {
     seq = (*seqIter);
     //     cout << "Name: " << seq.Name << "\tLength: " << seq.Length << endl;

     std::string sname = seq.Name;
     if(sname.compare(seqid) == 0){
       break;
     }
     i++;
   }

   std::cerr << "INFO: " << "scaffold: " << seq.Name << " ID: " << i << " Seq len: " << seq.Length << "\n"; 

   if(! reader.SetRegion(i, 0, i, boost::lexical_cast<int>(seq.Length))){
     std::cerr << "FATAL: cannot set region\n";    
     return 1;
   }

  map<int, int> flag_count;
    
  read_pileup zk;
  BamTools::BamAlignment al;
    
  while(reader.GetNextAlignment(al)){
    //    if(!al.IsPaired()){
    //  continue;
    //}
    flag_count[al.AlignmentFlag]++;
  }
  std::cerr << "flagstat done\n";
    
  printFlagCounts(flag_count);
  flagPileupInfo global = processFlagCounts(flag_count);
    
  BamTools::BamReader readerb;                                                                                
                                                                                                             
  if ( ! readerb.Open(filename)){                                                                             
    std::cerr << "FATAL: cannot open bam index\n";                                                           
    return 1;                                                                                                
  }                                                                                                          
  std::cerr << "INFO: " << filename << " opened for reading\n";                                              
  if( ! readerb.OpenIndex(fileindex)){                                                                        
    std::cerr << "FATAL: cannot open bam index\n";                                                           
    return 1;                                                                                                
  } 

    
  if(! readerb.SetRegion(i, 0, i, boost::lexical_cast<int>(seq.Length))){                                    
    std::cerr << "FATAL: cannot set region\n";                                                              
    return 1;                                                                                               
  }  


  while(readerb.GetNextAlignment(al)){
    // if(!al.IsPaired()){
    //  continue;
    // }

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
                 
      //      boost::math::binomial_distribution<float> pdf =  boost::math::binomial_distribution<float>(t, gp);

      //      float pp = boost::math::pdf(boost::math::binomial_distribution<float>(t, gp), b);
            
      float sign = 1;
            
      //if(p < gp){
      //	sign = -1;
      //}     
        
      float pp = boost::math::cdf(complement(boost::math::binomial_distribution<float>(t, gp), b));

                
      std::cout << seq.Name << "\t" << al.Position << "\t" << b << "\t" << t << "\t" << p << "\t" << pp << "\t" << gp << "\n";
    }
  }
  std::cerr << "pileup prop done\n";
}
