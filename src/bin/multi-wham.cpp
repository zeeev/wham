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
  vector<int> mapq;
  vector<int> flags;
};

//------------------------------------------------------------
//------------------------------------------------------------

vector<string> parse_groups(string  files){

  vector<string> split_files;
  split(split_files, files, is_any_of(","));

  return split_files;
}

//------------------------------------------------------------
//------------------------------------------------------------

void printvec(string base, vector<string> & data){

  for(vector<string>::iterator it = data.begin(); it != data.end(); it++){
    cerr << base  << ":  " << *it << endl;
  }
}
//------------------------------------------------------------
//------------------------------------------------------------

void check_index(vector<string> & group){

  BamMultiReader mreaderz;

  if(! mreaderz.Open( group ) ){
    cerr << "couldn't open file\n";
  }


  
  if(! mreaderz.HasIndexes()){
    mreaderz.CreateIndexes();
  }
   

  //cerr << "INFO: done checking index\n";
}

//------------------------------------------------------------
//------------------------------------------------------------

void load_info_struct(posInfo  *info, BamAlignment & read){

  (*info).nreads++;
  (*info).mapq.push_back(read.MapQuality);
  (*info).flags.push_back(read.AlignmentFlag);

  flag alflag;

  alflag.addFlag(read.AlignmentFlag);
  //  cout << "readflag: " << read.AlignmentFlag << endl;
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
void initPosInfo (posInfo * info){

  (*info).nreads       = 0;
  (*info).mateunmapped = 0;
  (*info).samestrand   = 0;

}
//------------------------------------------------------------
//------------------------------------------------------------
void process_pileup(list<BamAlignment> & data, map<string, int> & target_info){

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

  float t, b;

  t = (float)(target.mateunmapped) / (float)(target.nreads);
  b = (float)(background.mateunmapped) / (float)(background.nreads);

  cout  << t << "\t" << b << endl;

}

//------------------------------------------------------------
//------------------------------------------------------------

void pileup(BamMultiReader & mreader, map<string, int> & target_info){

  BamAlignment al;
  read_pileup PileUp;

  while(mreader.GetNextAlignment(al)){
    PileUp.proccess_alignment(al);
    // cout << al.Position << "\t" << al.Filename << endl;
    if(PileUp.currentStart() > PileUp.currentPos()){
      list<BamAlignment> dat =  PileUp.pileup();
      process_pileup(dat, target_info);
    } 

  }
}

//------------------------------------------------------------
//------------------------------------------------------------

void run_regions(vector<string> & target, vector <string> & background){

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
  }

  mreader.LocateIndexes();

  int nseqs = mreader.GetReferenceCount() -1;

  for(int i; i <= nseqs; i++){
    mreader.Jump(i, 0);
    pileup(mreader, target_info);
    cerr << "Finished seqid: " << i << endl;
  }


  mreader.Close();

}

//------------------------------------------------------------
//------------------------------------------------------------

int main(int argc,  char * argv[]){
    
  try{
  
    namespace po = boost::program_options;
    po::options_description desc ("Allowed options");
       
    desc.add_options()
      ("help",       "produce help info")
      ("target",     po::value<std::string>() ,  "The target bam files, comma sep list")
      ("background", po::value<std::string>() ,  "The background bam files, comma sep list");
    
    po::variables_map vm;
    
    try{
    
      po::store(po::parse_command_line(argc, argv, desc), vm);
      
      if(vm.count("help")){
	cout << "Usage: whammy -t a.bam,b.bam,c.bam -b d.bam,e.bam,f.bam" << endl;
	return SUCCESS;
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

    printvec("INFO: target bam",     target    );
    printvec("INFO: background bam", background);


    cerr << "INFO: starting to run regions\n";
    
    run_regions(target, background);

  }
  catch(std::exception& e){
    std::cerr << "Unhandled Exception reached the top of main: " 
	      << e.what() << ", application will now exit" << std::endl; 
    return ERROR_UNHANDLED_EXCEPTION; 
  }
  return SUCCESS;
}

//------------------------------------------------------------
//------------------------------------------------------------
