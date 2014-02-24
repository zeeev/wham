//
//  main.cpp
//  wham
//
//  Created by Zev Kronenberg on 12/31/12.
//  Copyright (c) 2012 Zev Kronenberg. All rights reserved.
//

#include <stdio.h>
#include "api/BamMultiReader.h"
#include "read_pileup.h"
#include "randomregion.h"
#include "flag.h"
#include "math.h"

#include "boost/math/distributions/normal.hpp"
#include "boost/math/distributions/negative_binomial.hpp"
#include "boost/math/distributions/poisson.hpp"
#include "boost/math/distributions/binomial.hpp"
#include "boost/math/distributions/chi_squared.hpp"
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
  int             nreads;
  int       mateunmapped;
  int         samestrand;
  int      otherscaffold;
  vector<int>       mapq;
  vector<int>      flags;
  vector<double>   fragl;
  map<string,int>  depth;
  vector<int>     depths;
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

double mean(vector<double> & dat){

  double n   = 0;
  double sum = 0;

  for(vector<double>::iterator datum = dat.begin(); datum != dat.end(); datum++){
    n += 1.0 ;  
    sum += *datum;
  
}

  double mean = sum / n;
  return mean;

}

//------------------------------------------------------------
double sd(vector<double> & dat, double mean){
  
  double ss = 0;
  double n  = 0;
  
  for(vector<double>::iterator datum = dat.begin(); datum != dat.end(); datum++){    
    ss += pow((*datum - mean), 2);
    n  += 1;
  }

  double sd = sqrt(ss/(n-1));
  return sd + 1;

}

//------------------------------------------------------------
double lldnorm(vector<double> & dat, double mu, double sdev){
  double llsum = 0;

  double ss = pow(sdev, 2);
  
  double p = mu / ss; 
  double r = pow(mu, 2) / (ss - mu);

  cerr << "INFO " << mu << "\t" << sdev << endl;
  //  cerr << "INFO " << mu << "\t" << sdev << "\t" << r << "\t" << p << "\t" << endl;

  for(vector<double>::iterator datum = dat.begin(); datum != dat.end(); datum++){

    //    if(r > 0){
    //  llsum += log ( boost::math::pdf(boost::math::negative_binomial_distribution<double>( r, p ), *datum) );
    // }
    // else{
      llsum += log ( boost::math::pdf(boost::math::normal_distribution<double>( mu, sdev ), *datum) );
      //}
  }

  return -1 * llsum;

}

//------------------------------------------------------------
void combine_info_struct(posInfo *t, posInfo *b, posInfo *a){
  (*a).nreads         =  (*t).nreads         +    (*b).nreads       ;
  (*a).otherscaffold  =  (*t).otherscaffold  +    (*b).otherscaffold; 
  (*a).samestrand     =  (*t).samestrand     +    (*b).samestrand   ;
  (*a).mateunmapped   =  (*t).mateunmapped   +    (*b).mateunmapped ;
  (*a).fragl.insert((*a).fragl.end(), (*t).fragl.begin(),(*t).fragl.end());
  (*a).fragl.insert((*a).fragl.end(), (*b).fragl.begin(),(*b).fragl.end());

}

//------------------------------------------------------------

/// load up the pileup information



void load_info_struct(posInfo  *info, BamAlignment & read, double rolling_mean){

  (*info).nreads++;
  (*info).mapq.push_back(read.MapQuality);
  (*info).flags.push_back(read.AlignmentFlag);

  (*info).depth[read.Filename]++;



  double ins =  abs(boost::lexical_cast<double>(read.InsertSize)) ;

    ins = abs(ins / rolling_mean);


  (*info).fragl.push_back(ins);

  flag alflag;

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
void loadDepths(posInfo *info){

  map<string, int>::iterator it;

  for(it = (*info).depth.begin(); it != (*info).depth.end(); ++it){    
    (*info).depths.push_back(it->second);
  }

}

//------------------------------------------------------------
double llbinom(posInfo *info , double mp, double sp, double op, int flag){
  
  double m = boost::lexical_cast<double> ((*info).mateunmapped);
  double o = boost::lexical_cast<double> ((*info).otherscaffold);
  double s = boost::lexical_cast<double> ((*info).samestrand);
  double n = boost::lexical_cast<double>((*info).nreads);
  
  if(flag == 1){
     mp = m / n;
     sp = s / n;
     op = o / n;
  }
  

  double psum = 0.0;

  psum += log(boost::math::pdf(boost::math::binomial_distribution<double>(n, mp), m ));
  psum += log(boost::math::pdf(boost::math::binomial_distribution<double>(n, op), o ));
  psum += log(boost::math::pdf(boost::math::binomial_distribution<double>(n, sp), s ));

  
  return( -1 * psum);

}


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

void process_pileup(list<BamAlignment> & data, map<string, int> & target_info, int pos, string & seqid, double rolling_meant, double rolling_meanb, double rolling_meana){

  //  cerr << "INFO: processing pileup\n"; 

  posInfo target, background, all;

  initPosInfo (&target);
  initPosInfo (&background);
  initPosInfo (&all);

  for(list<BamAlignment>::iterator read = data.begin(); read != data.end(); read++){
    int amItarget = target_info.count(read->Filename);
    if(amItarget == 0){
      load_info_struct( &background, *read, rolling_meanb);
    }
    else{
      load_info_struct( &target, *read, rolling_meant);
    }    
  }
  combine_info_struct(&target, &background, &all);

  if( background.nreads < 10 ){
    return;
  }
  if( target.nreads    < 10  ){
    return;
  }

  double nreads  = boost::lexical_cast<double>(all.nreads);
  double treads  = boost::lexical_cast<double>(target.nreads);
  double breads  = boost::lexical_cast<double>(background.nreads);
  double mp      = boost::lexical_cast<double>(all.mateunmapped) / nreads;
  double sp      = boost::lexical_cast<double>(all.samestrand) / nreads;
  double op      = boost::lexical_cast<double>(all.otherscaffold) / nreads;
  
  double fraglmu = mean(all.fragl);
  double fraglsd = sd(all.fragl,  fraglmu);

  if(fraglsd > 20){
    return;
  }

  double fglt = mean(target.fragl);
  double fglb = mean(background.fragl);
  double fgst = sd(target.fragl, fglt);
  double fgsb = sd(background.fragl, fglb);

  //  double FiST  =  (((treads/nreads) * fgst  ) + ((breads/nreads) * fgsb )) / fraglsd;

  double btn    = llbinom(&target,     mp, sp, op, 0);
  double bbn    = llbinom(&background, mp, sp, op, 0);
  double bta    = llbinom(&target,     mp, sp, op, 1);
  double bba    = llbinom(&background, mp, sp, op, 1);

  double fta    = lldnorm(target.fragl, fglt, fgst           );  
  double fba    = lldnorm(background.fragl, fglb, fgsb       );
  double ftn    = lldnorm(target.fragl, fraglmu, fraglsd     );  
  double fbn    = lldnorm(background.fragl, fraglmu, fraglsd );  


  //  double lrt =  2 * ((btn + bbn) - (bta + bba));
  // double lrtb = 2 * ((ftn + fbn) - (fta + fba)); 

  double lrt = 2 * ((btn + bbn + ftn + fbn) - (bta + bba + fta + fba));

  boost::math::chi_squared_distribution<double> chisq(5);
  //   boost::math::chi_squared_distribution<double> chisqb(3.0);

  if(lrt <= 0){
    lrt = 0.0001;
  }
  

  double lp  = 1 - boost::math::cdf(chisq, lrt);


  //  cout << pos << "\t" << "IM" << "\t" << fglt << "\t" << fglb << endl;
  // cout << pos << "\t" << "IS" << "\t" << fgst << "\t" << fgsb << endl;

  cout << seqid << "\t" << pos << "\t" << lrt << "\t" << lp << endl ;

}

//------------------------------------------------------------
//------------------------------------------------------------

/// run the pileup

void pileup(BamMultiReader & mreader, map<string, int> & target_info){

  BamAlignment al;
  read_pileup PileUp;

  BamTools::RefVector seqids = mreader.GetReferenceData();

  cerr << "INFO: starting pileup\n";


  double  Nt = 1;
  double  St = 1;
  double  Nb = 1;
  double  Sb = 1;
  double  Na = 1;
  double  Sa = 1;

  while(mreader.GetNextAlignment(al)){

    //    cerr << al.Position << "\n";

    //cerr << "DIST: " << al.InsertSize << "\n";

    if(al.IsDuplicate()){
      continue;
    }

    double ins = abs(boost::lexical_cast<double>(al.InsertSize));

    int amItarget = target_info.count(al.Filename);
    if(amItarget == 0){

      Nb += 1 ;
      Sb += ins;
    }
    else{
      Nt += 1 ;
      St += ins;
    }
    Na += 1;
    Sa += ins;

    double rolling_meant = St / Nt; 
    double rolling_meanb = Sb / Nb; 
    double rolling_meana = Sa / Na; 

    //    cerr << "ROLLING: " <<  rolling_meant << "\t" << rolling_meanb << "\t" <<  rolling_meana << endl;

    PileUp.proccess_alignment(al);
    string seqid = seqids[al.RefID].RefName;
    // cout << al.Position << "\t" << al.Filename << endl;                                                                     
    if(PileUp.currentStart() > PileUp.currentPos()){
      list<BamAlignment> dat =  PileUp.pileup();
      process_pileup(dat, target_info, PileUp.currentPos(), seqid, rolling_meant, rolling_meanb, rolling_meana);
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

   cerr << "INFO: entering check_index\n";

   check_index(total);

   cerr << "INFO: left check_index\n"; 

  BamMultiReader mreader;

   cerr << "INFO: opening BAMs and indices\n" << endl;

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

  cerr << "INFO region set to: " << seqid << endl;

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
