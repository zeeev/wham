//
//  main.cpp
//  wham
//
//  Created by Zev Kronenberg on 12/31/12.
//  Copyright (c) 2012 Zev Kronenberg. All rights reserved.
//

#include <stdio.h>
#include "api/BamMultiReader.h"
#include  "read_pileup.h"
#include  "randomregion.h"
#include  "flag.h"
#include "math.h"

#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/program_options.hpp"

using namespace std;
using namespace BamTools;
using namespace boost;
namespace 
{ 
  const size_t ERROR_IN_COMMAND_LINE     = 1; 
  const size_t SUCCESS                   = 0; 
  const size_t ERROR_UNHANDLED_EXCEPTION = 2;  
} 

struct posInfo
{
  int nreads;
  int mateunmapped;
  int samestrand;
  int otherscaffold;
  vector<int> mapq;
  vector<int> flags;
  vector<int> fragl;
};

//------------------------------------------------------------
//------------------------------------------------------------

/// split a string and return a vector

vector<string> parse_groups(string  files){

  vector<string> split_files;
  split(split_files, files, is_any_of(","));

  return split_files;
}

//------------------------------------------------------------
//------------------------------------------------------------

/// print a vector and a base string for the tag

void printvec(string base, vector<string> & data){

  for(vector<string>::iterator it = data.begin(); it != data.end(); it++){
    cerr << base  << ":  " << *it << endl;
  }
}
//------------------------------------------------------------
//------------------------------------------------------------

/// test if bam indecies exists and if they doesn't it creates them

void check_index(vector<string> & group){

  BamMultiReader mreaderz;

  if(! mreaderz.Open( group ) ){
    cerr << "couldn't open file\n";
  }
  
  if(! mreaderz.LocateIndexes()){
    cerr << "INFO: wham didn't find index, creating one.\n";
    mreaderz.CreateIndexes();
  }
   

  //cerr << "INFO: done checking index\n";
}

//------------------------------------------------------------
//------------------------------------------------------------

/// takes a vector and calculates the mean

float mean(vector<int> & dat){

  int n   = 0;
  int sum = 0;

  for(vector<int>::iterator datum = dat.begin(); datum != dat.end(); datum++){

    n++;
    sum += *datum;
      
  }

  float mean = static_cast<float>(sum) / static_cast<float>(n);
  
  return mean;

}

//------------------------------------------------------------
//------------------------------------------------------------

/// load up the pileup information

void load_info_struct(posInfo  *info, BamAlignment & read){

  (*info).nreads++;
  (*info).mapq.push_back(read.MapQuality);
  (*info).flags.push_back(read.AlignmentFlag);

  flag alflag;

  (*info).fragl.push_back(abs(read.InsertSize));

  alflag.addFlag(read.AlignmentFlag);
  //  cout << "readflag: " << read.AlignmentFlag << endl;

  if(!(read.RefID == read.MateRefID)){
    (*info).otherscaffold++;
  }

  if(! alflag.isPairMapped()){
    (*info).mateunmapped++;
    //  cout << "\tread pair is not mapped" << endl;
  }
  if(alflag.sameStrand()){
    (*info).samestrand++;
    // cout << "\tread is on the same strand" << endl;
  }
}

//------------------------------------------------------------
//------------------------------------------------------------

/// Initialize the pileup data struct

void initPosInfo (posInfo * info){

  (*info).nreads        = 0;
  (*info).mateunmapped  = 0;
  (*info).samestrand    = 0;
  (*info).otherscaffold = 0;

}
//------------------------------------------------------------
//------------------------------------------------------------

/// proces the pileup information

void process_pileup(list<BamAlignment> & data, map<string, int> & target_info, int pos, string & seqid){

  posInfo target, background;

  initPosInfo (&target);
  initPosInfo (&background);


  //  cout << "nreads: " << "\t" << data.size() << endl;


  for(list<BamAlignment>::iterator read = data.begin(); read != data.end(); read++){
    
    int amItarget = target_info.count(read->Filename);
    if(amItarget == 0){
      load_info_struct( &background, *read);
    }
    else{
      load_info_struct( &target, *read);
    }    
  }

  double t, b, tss, bss, tn, bn, tos, bos ;

  tn = static_cast<float>(target.nreads);
  bn = static_cast<float>(background.nreads);

  t = static_cast<float>(target.mateunmapped)     ;
  b = static_cast<float>(background.mateunmapped) ;
  
  tss = static_cast<float>(target.samestrand)    ;
  bss = static_cast<float>(background.samestrand);

  tos = static_cast<float>(target.otherscaffold);
  bos = static_cast<float>(background.otherscaffold);
    

  double mti = mean(target.fragl);
  double mbi = mean(background.fragl);

  double mtm = mean(target.mapq);
  double mbm = mean(background.mapq);

  double uclid = 0;

  double mdif = (mtm - mbm) ; // ((mbm+mtm)/2);
  double idif = (mti - mbi) ; // ((mbi + mti)/2);
  double ddif = (tn - bn  ) ;
  double ndif = (t - b    ) ;
  double sdif = (tss - bss) ;
  double odif = (tos - bos) ;

  if(! isnan(ddif)){
    // uclid += pow(ddif, 2);
  }
  if(! isnan(mdif)){
    // uclid += pow(idif, 2);
  }
  if(! isnan(mdif)){
    uclid += pow(mdif, 2);
  }
  if(! isnan(odif) ){
    uclid += pow(odif, 2);
  }
  if(! isnan(ndif)){
    uclid += pow(ndif, 2);
  }
  if(! isnan(sdif)){
    uclid += pow(sdif, 2);
  }
  double u = sqrt(uclid);

  cout << seqid << "\t" << pos << "\t" << t << "\t" << b << "\t" << tss << "\t" 
       << bss << "\t" << tos << "\t" << bos << "\t" << target.nreads << "\t" 
       << background.nreads << "\t" << mti << "\t" << mbi 
       << "\t" << mtm << "\t" << mbm << "\t" << u << endl;

}

//------------------------------------------------------------
//------------------------------------------------------------

/// run the pileup

void pileup(BamMultiReader & mreader, map<string, int> & target_info){

  BamAlignment al;
  read_pileup PileUp;

  BamTools::RefVector seqids = mreader.GetReferenceData();


  while(mreader.GetNextAlignment(al)){
    if(al.IsDuplicate()){
    	continue;
    }
    PileUp.proccess_alignment(al);
    string seqid = seqids[al.RefID].RefName;
    // cout << al.Position << "\t" << al.Filename << endl;
    if(PileUp.currentStart() > PileUp.currentPos()){
      list<BamAlignment> dat =  PileUp.pileup();
      process_pileup(dat, target_info, PileUp.currentPos(), seqid);
    } 

  }
}

//------------------------------------------------------------
//------------------------------------------------------------

/// run the regions

void run_regions(vector<string> & target, vector <string> & background, string & seqid, BamRegion & seqr){

  vector <string> total  = target;

  map <string, int> target_info;

  total.insert( total.end(), background.begin(), background.end() );

  for(vector<string>::iterator it = target.begin(); it != target.end(); it++){
    target_info[*it]++;
  }

   printvec("total grouping", total);

  // cerr << "INFO: entering check_index\n";

   check_index(total);

  // cerr << "INFO: left check_index\n"; 

  BamMultiReader mreader;

  // cerr << "INFO: opening BAMs and indices\n" << endl;

  if(! mreader.Open(total)){
    cerr << "cannot open bams.\n";
    exit(EXIT_FAILURE);
  }

  if(! mreader.LocateIndexes()){
    cerr << "cannot create or locate indicies" << endl;
    exit(EXIT_FAILURE);  
  }

  RefVector seqids = mreader.GetReferenceData();

  int nseqs = mreader.GetReferenceCount() -1;

  if(seqid.compare("NA") == 0 ){
    pileup(mreader, target_info);
    cerr << "Finished Whole Genome!" << endl;
  }
  else{
    if(! mreader.SetRegion(seqr)){
      cerr << "failed to jump to seqid!" << endl;
      exit(EXIT_FAILURE);
    }
    pileup(mreader, target_info);
  }
  mreader.Close();

}

//------------------------------------------------------------
//------------------------------------------------------------

/// given a seqid and a file/set of files find the index for the seqid

BamRegion getseqidn (string & seqid, vector<string> & target){

  BamRegion region;

  cerr << "WTF: " << seqid << endl;

  BamMultiReader mreader;

  if(! mreader.Open(target)){
    cerr << "cannot open bams.\n";
  }
  
  int i = 0 ;

  BamTools::SamHeader           header = mreader.GetHeader();
  BamTools::SamSequenceDictionary seqs = header.Sequences;

  BamTools::SamSequence seq;
  BamTools::SamSequenceConstIterator seqIter = seqs.ConstBegin();
  BamTools::SamSequenceConstIterator seqEnd  = seqs.ConstEnd();

  for ( ; seqIter != seqEnd; ++seqIter ) {

    seq = (*seqIter);
    std::string sname = seq.Name;
    if(sname.compare(seqid) == 0){
      cerr << "seqid: " << sname << "\t" << "is index: " << "\t" << i << endl;

      region.LeftRefID     = i;
      region.RightRefID    = i;
      region.LeftPosition  = 0;
      region.RightPosition = boost::lexical_cast<int>(seq.Length);

      break;
    }
    i++;
  }
  
  return region;

}

//------------------------------------------------------------
//------------------------------------------------------------

int main(int argc,  char * argv[]){

  string seqid = "NA";
  BamRegion seqr;
  
  try{
  
    namespace po = boost::program_options;
    po::options_description desc ("Allowed options");
       
    desc.add_options()
      ("help,h",       "produce help info")
      ("target,t",     po::value<string>() ,   "The target bam files, comma sep list")
      ("background,b", po::value<string>() ,   "The background bam files, comma sep list")
      ("seqid,s",      po::value<string>() ,   "Confine the analysis to a single seqid" );
      
    po::variables_map vm;
    
    try{
    
      po::store(po::parse_command_line(argc, argv, desc), vm);
      
      if(vm.count("help")){
	cout << "Usage: whammy -t a.bam,b.bam,c.bam -b d.bam,e.bam,f.bam" << endl;
	return SUCCESS;
      }
      if(! vm.count("target")){
	cout << "failure to specify target correctly" << endl;
	cout << "Usage: whammy -t a.bam,b.bam,c.bam -b d.bam,e.bam,f.bam" << endl;
	return ERROR_IN_COMMAND_LINE;
      }
      if(! vm.count("background")){
      	cout << "failure to specify background correctly" << endl;
	cout << "Usage: whammy -t a.bam,b.bam,c.bam -b d.bam,e.bam,f.bam" << endl;
	return ERROR_IN_COMMAND_LINE;
      }
      po::notify(vm);
    }
    catch(po::error & e){
      cerr << "ERROR: " << e.what() << endl << endl;
      cerr << desc << endl;
      return ERROR_IN_COMMAND_LINE;	  
    }
 
    vector<string> target     = parse_groups(vm["target"].as<string>());
    vector<string> background = parse_groups(vm["background"].as<string>());

    if(vm.count("seqid")){
      seqid = vm["seqid"].as<string>();
      seqr  = getseqidn(seqid, target);
    }

    //    cerr << "if region is set: " << "\t" << seqn << "\t" << seqid << endl;

    printvec("INFO: target bam",     target    );
    printvec("INFO: background bam", background);
    cerr << "INFO: starting to run regions\n";
    
    run_regions(target, background, seqid, seqr);
    cerr << "INFO: Whamy has finished normally\n";
    return 0;    
  }
  catch(std::exception& e){
    std::cerr << "Unhandled Exception reached the top of main: " 
	      << e.what() << ", application will now exit" << std::endl; 
    return ERROR_UNHANDLED_EXCEPTION; 
  }

  
  

}

//------------------------------------------------------------
//------------------------------------------------------------
