//
//  main.cpp
//  wham
//
//  Created by Zev Kronenberg on 12/31/12.
//  Copyright (c) 2012 Zev Kronenberg. All rights reserved.
//


#include "stdio.h"
#include "api/BamMultiReader.h"
#include "read_pileup.h"
#include "randomregion.h"
#include "entropy.h"
#include "flag.h"
#include "math.h"
#include <stdio.h>      
#include <stdlib.h>     
#include <time.h>       


#include "boost/asio/io_service.hpp"
#include "boost/bind.hpp"
#include "boost/thread/thread.hpp"

#include "boost/math/distributions/normal.hpp"
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
  double           fragm;
  vector<double>   fragl;
};


boost::mutex print_guard;

double running_mb = 0;


//------------------------------------------------------------

/// split a string and return a vector

vector<string> parse_groups(string  files){

  vector<string> split_files;
  split(split_files, files, is_any_of(","));

  return split_files;
}

//------------------------------------------------------------

/// print a vector

void printansvec(string base, vector<string> & data){
  for(vector<string>::iterator it = data.begin(); it != data.end(); it++){
    cout <<  *it ;
  }
}

//------------------------------------------------------------

/// print a vector and a base string for the tag

void printvec(string base, vector<string> & data){
  for(vector<string>::iterator it = data.begin(); it != data.end(); it++){
    cerr << base  << ":  " << *it << endl;
  }
}

//------------------------------------------------------------

/// test if bam indecies exists and if they doesn't it creates them

void check_index(vector<string> & group){

  BamMultiReader mreaderz;

  if(! mreaderz.Open( group ) ){
    cerr << "couldn't open bams" << endl;
  }
  
  if(! mreaderz.LocateIndexes()){
    cerr << "INFO: wham didn't find indices, creating them." << endl;
    mreaderz.CreateIndexes();
  }
   
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

/// take a vector and a mean and calculate the standard deviation

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

/// loop over a vector of data and return the total log likelihood

double lldnorm(vector<double> & dat, double mu, double sdev){
  double llsum = 0;

  for(vector<double>::iterator datum = dat.begin(); datum != dat.end(); datum++){

      llsum += log ( boost::math::pdf(boost::math::normal_distribution<double>( mu, sdev ), *datum) );
  }
  return -1 * llsum;
}

//------------------------------------------------------------

/// take target and background information from the pileup and combine into total

void combine_info_struct(posInfo *t, posInfo *b, posInfo *a){
  (*a).nreads         =  (*t).nreads         +    (*b).nreads       ;
  (*a).otherscaffold  =  (*t).otherscaffold  +    (*b).otherscaffold; 
  (*a).samestrand     =  (*t).samestrand     +    (*b).samestrand   ;
  (*a).mateunmapped   =  (*t).mateunmapped   +    (*b).mateunmapped ;
  (*a).fragm          =  (*t).fragm          +    (*b).fragm;
  (*a).fragl.insert((*a).fragl.end(), (*t).fragl.begin(),(*t).fragl.end());
  (*a).fragl.insert((*a).fragl.end(), (*b).fragl.begin(),(*b).fragl.end());

}

//------------------------------------------------------------

/// load up the pileup information

void load_info_struct(posInfo  *info, BamAlignment & read){

  (*info).nreads++;

  double ins =  abs(double(read.InsertSize)) ;
  
  (*info).fragm += ins;

  //  ins = abs(ins / rolling_mean);

  (*info).fragl.push_back(ins);
  
  flag alflag;

  alflag.addFlag(read.AlignmentFlag);
  
  if(!(read.RefID == read.MateRefID)){
    (*info).otherscaffold++;
  }
  if(! alflag.isPairMapped()){
    (*info).mateunmapped++;
  }
  if(alflag.sameStrand()){
    (*info).samestrand++;
  }
}

//------------------------------------------------------------

/// give the paramters of the binomial mp,sp,op calculate total log likelihood

double llbinom(posInfo *info , double mp, double sp, double op, int flag){
  
  double m = double (((*info).mateunmapped));
  double o = double (((*info).otherscaffold));
  double s = double (((*info).samestrand));
  double n = double (((*info).nreads));
  
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
  (*info).fragm         = 0;

}
//------------------------------------------------------------

/// proces the pileup information

double score(vector<BamAlignment> & dat, map<string, int> & target_info ){

  posInfo target, background, all;

  initPosInfo (&target);
  initPosInfo (&background);
  initPosInfo (&all);

  double nreads = double (dat.size());

  map<string, int>::iterator amItarget;

  for(int d = 0; d < nreads; d++){
    amItarget = target_info.find(dat[d].Filename);
    if(amItarget == target_info.end()){
      load_info_struct( &background, dat[d]);
    }
    else{
      load_info_struct( &target, dat[d]);
    }    
  }
  combine_info_struct(&target, &background, &all);

  double treads  = double (target.nreads);
  double breads  = double (background.nreads);


  if(treads < 10){
    return 1;
  }
  if(breads < 10){
    return 1;
  }

  double mp      = double (all.mateunmapped) / nreads;
  double sp      = double (all.samestrand) / nreads;
  double op      = double (all.otherscaffold) / nreads;
  

  double fraglmu = all.fragm / nreads;
  double fraglsd = sd(all.fragl,  fraglmu);

  double fglt = target.fragm / treads ;  
  double fglb = background.fragm / breads; 
  double fgst = sd(target.fragl, fglt);
  double fgsb = sd(background.fragl, fglb);


  double btn    = llbinom(&target,     mp, sp, op, 0);
  double bbn    = llbinom(&background, mp, sp, op, 0);
  double bta    = llbinom(&target,     mp, sp, op, 1);
  double bba    = llbinom(&background, mp, sp, op, 1);

  double fta    = lldnorm(target.fragl, fglt, fgst           );  
  double fba    = lldnorm(background.fragl, fglb, fgsb       );
  double ftn    = lldnorm(target.fragl, fraglmu, fraglsd     );  
  double fbn    = lldnorm(background.fragl, fraglmu, fraglsd );  

  double lrt = 2 * ((btn + bbn + ftn + fbn) - (bta + bba + fta + fba));

  if(isinf(lrt)){
    return 1;
  }
  if(isnan(lrt)){
    return 1;
  }
  return lrt;
}

//------------------------------------------------------------                                                                

/// permute the score function by suffling the reads' filenames
       
double permute(double lrt, vector<BamAlignment> & dat, map<string, int> & target_info){

  double success = 0;
  double ntrials = 0;

  while(ntrials < 100000 && success < 20){
    ntrials += 1;

    for(int i = dat.size() -1 ;  i > 0 ; --i){
      
      if(target_info.count(dat[i].Filename)){
	swap(dat[i].Filename, dat[ rand() % dat.size() ].Filename);
      }
    }

    double newlrt = score(dat, target_info);
    if(newlrt > lrt){
      success += 1;
    }
  }
  return success / ntrials;
}

//------------------------------------------------------------ 

/// load random LRTS 

void pileupLRT(int s, int j, int e,  map <string, int> target_info, vector<string> total, vector<double>& LRT){

  BamMultiReader mreader_thread;
  BamRegion      region_thread;

  if(! mreader_thread.Open(total)){
    cerr << "ERROR: cannot open bams" << endl;
    exit(EXIT_FAILURE);
  }

  if(! mreader_thread.LocateIndexes()){
    cerr << "ERROR: cannot create or locate indicies" << endl;
    exit(EXIT_FAILURE);
  }

  region_thread.LeftRefID     = s;
  region_thread.RightRefID    = s;
  region_thread.LeftPosition  = j;
  region_thread.RightPosition = e;
  
  if(! mreader_thread.SetRegion(region_thread)){
    cerr << "ERROR: a thread was unable to set region" << endl;
    exit(EXIT_FAILURE);
  }

  BamAlignment al;
  read_pileup PileUp;

  BamTools::RefVector seqids = mreader_thread.GetReferenceData();

  std::vector <string> buffer;

  double nreads = 0;

  map <string , int> chucker;
  
  while(mreader_thread.GetNextAlignment(al)){
    nreads += 1;

    string rname = al.Name;

    map<string, int>::iterator check = chucker.find(rname);

    if(check == chucker.end()){
      size_t found = al.QueryBases.find("N");
      if(found != std::string::npos){
        chucker[rname] = 1;
        continue;
      }
      fastQ read_dat;
      read_dat.setDNA(al.QueryBases);
      double en = read_dat.entropy(4);
      if(en < 1){
        chucker[rname] = 1;
        continue;
      }
    }
    else{
      chucker.erase(rname);
      continue;
    }
    
    if(! al.IsMapped()){
      continue;
    }
    if(! al.IsPrimaryAlignment()){
      continue;
    }
 
    PileUp.proccess_alignment(al);
    string seqid = seqids[al.RefID].RefName;
    if(PileUp.currentStart() > PileUp.currentPos()){
      list<BamAlignment> dat =  PileUp.pileup();
      
      if(dat.size() < 20){
	continue;
      }

      vector <BamAlignment> data(dat.begin(), dat.end());
      
      double results = score(data, target_info);
      
      boost::math::chi_squared_distribution<double> chisq(5);
      
       if(results < 2){
 	continue;
       }
 
        double pv = 1 - boost::math::cdf(chisq, results);
        double pr = 1;
       
        if(pv > 0.01){
       	continue;
        }
      

      
      LRT.push_back(results);
    }
  }
}

//------------------------------------------------------------

/// run the pileup

void pileup(int s, int j, int e,  map <string, int> target_info, vector<string> total, double cut){

  BamMultiReader mreader_thread;
  BamRegion      region_thread;

  if(! mreader_thread.Open(total)){
    cerr << "ERROR: cannot open bams" << endl;
    exit(EXIT_FAILURE);
  }

  if(! mreader_thread.LocateIndexes()){
    cerr << "ERROR: cannot create or locate indicies" << endl;
    exit(EXIT_FAILURE);
  }

  region_thread.LeftRefID     = s;
  region_thread.RightRefID    = s;
  region_thread.LeftPosition  = j;
  region_thread.RightPosition = e;

 if(! mreader_thread.SetRegion(region_thread)){
   cerr << "ERROR: a thread was unable to set region" << endl;
   exit(EXIT_FAILURE);
 }

  BamAlignment al;
  read_pileup PileUp;

  BamTools::RefVector seqids = mreader_thread.GetReferenceData();

  std::vector <string> buffer;

  double nreads = 0;

  map <string , int> chucker;

  while(mreader_thread.GetNextAlignment(al)){
    nreads += 1;
    
    string rname = al.Name;

    map<string, int>::iterator check = chucker.find(rname);
    
    if(check == chucker.end()){
     
      size_t found = al.QueryBases.find("N");
      if(found != std::string::npos){
	chucker[rname] = 1;
	continue;	
      }

      fastQ read_dat;
      read_dat.setDNA(al.QueryBases);
      double en = read_dat.entropy(4);     

      if(en < 1){
        chucker[rname] = 1;
        continue;
      }

    }
    else{
      chucker.erase(rname);
      continue;
    }

    if(! al.IsMapped()){
      continue;
    }
    if(! al.IsPrimaryAlignment()){
      continue; 
    }

    PileUp.proccess_alignment(al);
    string seqid = seqids[al.RefID].RefName;
    if(PileUp.currentStart() > PileUp.currentPos()){
      list<BamAlignment> dat =  PileUp.pileup();

      if(dat.size() < 20){
	continue;
      }

      vector <BamAlignment> data(dat.begin(), dat.end());
      
      double results = score(data, target_info);

      boost::math::chi_squared_distribution<double> chisq(5);
      
      if(results < 2){
	continue;
      }

      double pv = 1 - boost::math::cdf(chisq, results);
      string pr = "NA";
      
      if(pv > 0.01){
      	continue;
      }
      if(results > cut){
      	double pvp = permute(results, data, target_info);
	pr = lexical_cast<string>(pvp);
      }

      string ans;
      ans.append(seqid);
      ans.append("\t");
      ans.append(lexical_cast<string>(al.Position));
      ans.append("\t");
      ans.append(lexical_cast<string>(data.size()));
      ans.append("\t");
      ans.append(lexical_cast<string>(results));
      ans.append("\t");
      ans.append(lexical_cast<string>(pv));
      ans.append("\t");
      ans.append(pr);
      ans.append("\n");

      buffer.push_back(ans);
    
    }
  }

  mreader_thread.Close();
  
  print_guard.lock();
  running_mb += (nreads / 1000000);
  cerr << "INFO:" << running_mb << " Million reads finished" << endl; 
  printansvec("", buffer);
  print_guard.unlock();
  

 
}

//------------------------------------------------------------

/// set the mean mapping distance between matepairs

void set_mu_i(map<string, int> & target_info, vector<string> & all, double *mut, double *mub){

  BamMultiReader mreader;

  if(! mreader.Open(all)){
    cerr << "cannot open bams." << endl;
  }

  double bs = 0;
  double bn = 0;

  double ts = 0;
  double tn = 0;

  BamAlignment al;

  while(mreader.GetNextAlignment(al)){
    
    if(! al.IsMapped()){
      continue;
    }

    if(! al.IsPrimaryAlignment()){
      continue;
    }

    double ins = abs(boost::lexical_cast<double>(al.InsertSize));

    map<string, int>::iterator amItarget = target_info.find( al.Filename );
    if(amItarget == target_info.end() ){
      bs += ins;
      bn += 1;
    }
    else{
      ts += ins;
      tn += 1;
    }
    if(bn > 1000000 && tn > 1000000){
      break;
    }
  }

  (*mut) = ts / tn;
  (*mub) = bs / bn;

  mreader.Close();

}
//------------------------------------------------------------

/// run the regions

void run_regions(vector<string> & target, vector <string> & background, string & seqid, int seqr, int cpu){

  vector <string> total  = target;

  map<string, int> target_info;

  total.insert( total.end(), background.begin(), background.end() );

  for(vector<string>::iterator it = target.begin(); it != target.end(); it++){
    target_info[*it] = 1;
  }

  cerr << "INFO: locating and testing bam indices" << endl;

   check_index(total);

   cerr << "INFO: bam indices checked out" << endl; 

   BamMultiReader mreader;

   cerr << "INFO: opening bam files for reading" << endl;

  if(! mreader.Open(total)){
    cerr << "ERROR: cannot open bams." << endl;
    exit(EXIT_FAILURE);
  }

  double mu_i_target    ;
  double mu_i_background;

  // set_mu_i(target_info, total, & mu_i_target, & mu_i_background);

  //   cerr << "INFO: target mate pair mapping distance average over 1million reads: "     << mu_i_target     << endl;
  //   cerr << "INFO: background mate pair mapping distance average over 1million reads: " << mu_i_background << endl;

  if(! mreader.LocateIndexes()){
    cerr << "ERROR: cannot create or locate indicies" << endl;
    exit(EXIT_FAILURE);  
  }

  RefVector seqids = mreader.GetReferenceData();

  int nseqs = mreader.GetReferenceCount() -1;
  
  SamHeader           header = mreader.GetHeader();
  SamSequenceDictionary seqs = header.Sequences;
  
  vector<double> randLRT;

  SamSequence seq;

  while(randLRT.size() < 100000){
   
    int rseqid = rand() % nseqs;
    int ni     = 0;

    SamSequenceConstIterator seqIter = seqs.ConstBegin();
    SamSequenceConstIterator seqEnd  = seqs.ConstEnd();
       
    for ( ; seqIter != seqEnd; ++seqIter ) {
      if(ni != rseqid){
	ni += 1;
	continue;
      }
      else{
	seq = (*seqIter);
	break;
      }
    }

    int slen = lexical_cast<int>(seq.Length);
    int rstart = rand() % (slen - 50000);
    if(rstart < 0){
      	continue;
    }
    
    int rend   = rstart + 50000;
   
    pileupLRT(rseqid, rstart, rend, target_info, total, randLRT );
    double per = double(randLRT.size());

    cerr << "INFO: " << "assayed " << (100 * (per / 50000)) << " % of LRT genomic baseline" << endl;
  }

  BamTools::SamSequenceConstIterator seqIter = seqs.ConstBegin();
  BamTools::SamSequenceConstIterator seqEnd  = seqs.ConstEnd();


  double rmu = mean(randLRT);
  double rsd = sd(randLRT, rmu);
  
  double cut = rmu + (3*rsd);
    
  asio::io_service io_service;
  auto_ptr<asio::io_service::work> work(new asio::io_service::work(io_service));

  boost::thread_group threads;
  
  for(int i = 0; i < cpu; i++){
    threads.create_thread(boost::bind(&asio::io_service::run, &io_service));
  }
  
  cerr << "INFO: launched " << cpu << " threads" << endl;
  
  int s      = 0;
  double sum = 0;
    
  for ( ; seqIter != seqEnd; ++seqIter ) {
    if(seqr != -5){
      if( seqr != s ){
	s += 1;
	  continue;
      }
    }
    
    seq = (*seqIter);
    std::string sname = seq.Name;
    
    cerr << "INFO: loading " << sname << " into thread pool " << endl;
    int j = 0;
    
    for (; j + 2000000 <=  boost::lexical_cast<int>(seq.Length); j += 2000000){
      io_service.post(boost::bind(pileup, s, j, j+2000000, target_info, total, cut)); 
    }
    io_service.post(boost::bind(pileup, s, j, boost::lexical_cast<int>(seq.Length), target_info, total, cut)); 
    s += 1.0;
  }
  
  work.reset();

  threads.join_all();
  
  mreader.Close();

}

//------------------------------------------------------------


/// given a seqid and a file/set of files find the index for the seqid

int getseqidn (string & seqid, vector<string> & target){

  BamRegion region;

  cerr << "INFO region set to: " << seqid << endl;

  BamMultiReader mreader;

  if(! mreader.Open(target)){
    cerr << "cannot open bams." << endl;
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
      break;
    }
    i++;
  }
  return i;
}

//------------------------------------------------------------

/// main's main of course.  Using boost program options flag errors.

int main(int argc,  char * argv[]){

  srand(time(NULL));

  string seqid = "NA";
  int seqr     = -5;
  int cpu      = 2;

  try{
  
    namespace po = boost::program_options;
    po::options_description desc ("Allowed options");
       
    desc.add_options()
      ("help,h",       "produce help info")
      ("target,t",     po::value<string>() ,   "The target bam files, comma sep list")
      ("background,b", po::value<string>() ,   "The background bam files, comma sep list")
      ("seqid,s",      po::value<string>() ,   "Confine the analysis to a single seqid" )
      ("cpu,c",        po::value<int>()    ,      "The number of threads to use");
    po::variables_map vm;
    
    try{
    
      po::store(po::parse_command_line(argc, argv, desc), vm);
      
      if(vm.count("help")){
	cout << "Usage: raw -t a.bam,b.bam,c.bam -b d.bam,e.bam,f.bam" << endl;
	return SUCCESS;
      }
      if(! vm.count("target")){
	cout << "failure to specify target correctly" << endl;
	cout << "raw -c 5 -t a.bam,b.bam,c.bam -b d.bam,e.bam,f.bam" << endl;
	return ERROR_IN_COMMAND_LINE;
      }
      if(! vm.count("background")){
      	cout << "failure to specify background correctly" << endl;
	cout << "raw -c 5 -t a.bam,b.bam,c.bam -b d.bam,e.bam,f.bam" << endl;
	return ERROR_IN_COMMAND_LINE;
      }
      if(! vm.count("cpu")){
        cout << "ERROR: failure to specify cpus correctly" << endl;
        cout << "Usage: raw -c 5 -t a.bam,b.bam,c.bam -b d.bam,e.bam,f.bam" << endl;
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

    cpu = vm["cpu"].as<int>();

    if(vm.count("seqid")){
      seqid = vm["seqid"].as<string>();
      seqr  = getseqidn(seqid, target);
    }

    //    cerr << "if region is set: " << "\t" << seqn << "\t" << seqid << endl;

    printvec("INFO: target bam",     target    );
    printvec("INFO: background bam", background);
    cerr <<  "INFO: starting to run regions" << endl;
    
    run_regions(target, background, seqid, seqr, cpu);
    cerr << "INFO: raw has finished" << endl;
    return 0;    
  }
  catch(std::exception& e){
    std::cerr << "Unhandled Exception reached the top of main: " 
	      << e.what() << ", application will now exit" << endl; 
    return ERROR_UNHANDLED_EXCEPTION; 
  }
}
