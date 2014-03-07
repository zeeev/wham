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
#include "flag.h"
#include "math.h"

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
  vector<double>   fragl;
};


boost::mutex print_guard;

double running_mb = 0;


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

void printansvec(string base, vector<string> & data){
  for(vector<string>::iterator it = data.begin(); it != data.end(); it++){
    cout <<  *it ;
  }
}

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
    cerr << "couldn't open file" << endl;
  }
  
  if(! mreaderz.LocateIndexes()){
    cerr << "INFO: wham didn't find index, creating one." << endl;
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

  for(vector<double>::iterator datum = dat.begin(); datum != dat.end(); datum++){

      llsum += log ( boost::math::pdf(boost::math::normal_distribution<double>( mu, sdev ), *datum) );
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

  double ins =  abs(boost::lexical_cast<double>(read.InsertSize)) ;
  
     ins = abs(ins / rolling_mean);

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
double bound(double v){
  if(v <= 0){
    return 0.000001;
  }
  if(v >= 1){
    return 0.999999;
  }
  return v;
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
  //  psum += log(boost::math::pdf(boost::math::binomial_distribution<double>(n, op), o ));
  //  psum += log(boost::math::pdf(boost::math::binomial_distribution<double>(n, sp), s ));

  
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
/// proces the pileup information

double score(vector<BamAlignment> & dat, map<string, int> & target_info , double rolling_meant, double rolling_meanb){

  posInfo target, background, all;

  initPosInfo (&target);
  initPosInfo (&background);
  initPosInfo (&all);

  for(int d = 0; d < dat.size(); d++){
    int amItarget = target_info.count(dat[d].Filename);
    if(amItarget == 0){
      load_info_struct( &background, dat[d], rolling_meanb);
    }
    else{
      load_info_struct( &target, dat[d], rolling_meant);
    }    
  }
  combine_info_struct(&target, &background, &all);

  double nreads  = boost::lexical_cast<double>(all.nreads);
  double treads  = boost::lexical_cast<double>(target.nreads);
  double breads  = boost::lexical_cast<double>(background.nreads);
  double mp      = boost::lexical_cast<double>(all.mateunmapped) / nreads;
  double sp      = boost::lexical_cast<double>(all.samestrand) / nreads;
  double op      = boost::lexical_cast<double>(all.otherscaffold) / nreads;
  
  if(treads < 10 ){
    return 1;
  }
  if(breads < 10 ){
    return 1;
  }

  double fraglmu = mean(all.fragl);
  double fraglsd = sd(all.fragl,  fraglmu);

  double fglt = mean(target.fragl);
  double fglb = mean(background.fragl);
  double fgst = sd(target.fragl, fglt);
  double fgsb = sd(background.fragl, fglb);

  if(abs(fglt - fglb) < 0.05 && (mp + sp + op) < 0.05){
    return 1;
  }

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
//------------------------------------------------------------                                                                       
double permute(double lrt, vector<BamAlignment> & dat, map<string, int> & target_info, double rolling_meant, double rolling_meanb){

  double success = 0;

  for(int rep = 0; rep < 1000; rep++){

    for(int i = dat.size() -1 ;  i > 0 ; --i){
      
      if(target_info.count(dat[i].Filename)){
	swap(dat[i].Filename, dat[ rand() % dat.size() ].Filename);
      }
    }

    double newlrt = score(dat, target_info, rolling_meant, rolling_meanb);
    //    cerr << "newlrt" << "\t" << newlrt << endl;
    if(newlrt > lrt){
      success += 1;
    }
  }
  return success / 1000;
}

//------------------------------------------------------------ 
//------------------------------------------------------------

/// run the pileup

void pileup(int s, int j, int e,  map <string, int> target_info, vector<string> total, double rolling_meant, double rolling_meanb){

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

  std::map<string, int> readname;

  while(mreader_thread.GetNextAlignment(al)){
  
    if(! al.IsMapped()){
      continue;
    }
    if(! al.IsPrimaryAlignment()){
      continue; 
    }
    if(al.MapQuality < 20){
      continue;
    }
    if(! al.IsFirstMate()){
      continue;
    }
    if(al.IsDuplicate()){
      continue;
    }
   
    PileUp.proccess_alignment(al);
    string seqid = seqids[al.RefID].RefName;
    if(PileUp.currentStart() > PileUp.currentPos()){
      list<BamAlignment> dat =  PileUp.pileup();

      vector <BamAlignment> data(dat.begin(), dat.end());
      
      //      cerr << "about to score\n";

      double results = score(data, target_info, rolling_meant, rolling_meanb);

      //      cerr << "about to score\n";

      boost::math::chi_squared_distribution<double> chisq(5);
      
      if(results < 0){
	results = 1;
      }

      double pv = 1 - boost::math::cdf(chisq, results);
      double pr = 1;
      

      if(pv < 0.01){
	pr = permute(results, data, target_info, rolling_meant, rolling_meanb);
      }
      
      if(pv > 0.05){
	continue;
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
      ans.append(lexical_cast<string>(pr));
      ans.append("\n");

      buffer.push_back(ans);
    
    }
  }

  print_guard.lock();
  running_mb += (lexical_cast<double>(buffer.size()) / 1000000);
  cerr << "INFO:" << running_mb << "Mb finished" << endl; 
  if(buffer.size() > 0){
    printansvec("", buffer);
    print_guard.unlock();
    buffer.clear();
  }
}

//------------------------------------------------------------
//------------------------------------------------------------
void set_mu_i(map<string, int> & target_info, vector<string> & target, double *mut, double *mub){

  BamMultiReader mreader;

  if(! mreader.Open(target)){
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

    if(al.MapQuality < 20){
      continue;
    }

    if(! al.IsFirstMate()){
      continue;
    }

    if(al.IsDuplicate()){
      continue;
    }

    double ins = abs(boost::lexical_cast<double>(al.InsertSize));

    int amItarget = target_info.count(al.Filename);
    if(amItarget == 0){
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

  map <string, int> target_info;

  total.insert( total.end(), background.begin(), background.end() );

  for(vector<string>::iterator it = target.begin(); it != target.end(); it++){
    target_info[*it]++;
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

  set_mu_i(target_info,total, & mu_i_target, & mu_i_background);

   cerr << "INFO: target mate pair mapping distance average over 1million reads: "     << mu_i_target     << endl;
   cerr << "INFO: background mate pair mapping distance average over 1million reads: " << mu_i_background << endl;

  if(! mreader.LocateIndexes()){
    cerr << "ERROR: cannot create or locate indicies" << endl;
    exit(EXIT_FAILURE);  
  }

  RefVector seqids = mreader.GetReferenceData();

  int nseqs = mreader.GetReferenceCount() -1;
  
  BamTools::SamHeader           header = mreader.GetHeader();
  BamTools::SamSequenceDictionary seqs = header.Sequences;
  
  BamTools::SamSequence seq;
  BamTools::SamSequenceConstIterator seqIter = seqs.ConstBegin();
  BamTools::SamSequenceConstIterator seqEnd  = seqs.ConstEnd();
     
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
      io_service.post(boost::bind(pileup, s, j, j+2000000, target_info, total, mu_i_target, mu_i_background)); 
    }
    io_service.post(boost::bind(pileup, s, j, boost::lexical_cast<int>(seq.Length), target_info, total, mu_i_target, mu_i_background)); 
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
//------------------------------------------------------------

int main(int argc,  char * argv[]){

  string seqid = "NA";
  int seqr = -5;
  int cpu = 2;

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
	cout << "Usage: whammy -t a.bam,b.bam,c.bam -b d.bam,e.bam,f.bam" << endl;
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

//------------------------------------------------------------
//------------------------------------------------------------
