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
  //  vector<int>       mapq;
  // vector<int>      flags;
  vector<double>   fragl;
  // map<string,int>  depth;
  // vector<int>     depths;
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
  //  (*info).mapq.push_back(read.MapQuality);
  //  (*info).flags.push_back(read.AlignmentFlag);

  //  (*info).depth[read.Filename]++;

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
// void loadDepths(posInfo *info){
// 
//   map<string, int>::iterator it;
// 
//   for(it = (*info).depth.begin(); it != (*info).depth.end(); ++it){    
//     (*info).depths.push_back(it->second);
//   }
// 
// }

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

std::string process_pileup(list<BamAlignment> & data, map<string, int> & target_info, int pos, string & seqid, double rolling_meant, double rolling_meanb, double rolling_meana){

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

  double nreads  = boost::lexical_cast<double>(all.nreads);
  double treads  = boost::lexical_cast<double>(target.nreads);
  double breads  = boost::lexical_cast<double>(background.nreads);
  double mp      = boost::lexical_cast<double>(all.mateunmapped) / nreads;
  double sp      = boost::lexical_cast<double>(all.samestrand) / nreads;
  double op      = boost::lexical_cast<double>(all.otherscaffold) / nreads;
  
  double fraglmu = mean(all.fragl);
  double fraglsd = sd(all.fragl,  fraglmu);

  if(fraglsd <= 0){
    fraglsd = 0.0001;
  }

  if(treads < 10){
    return "";
  }
  if(breads < 10){
    return "";
  }
  if(fraglsd > 20){
    return "";
  }

  double fglt = mean(target.fragl);
  double fglb = mean(background.fragl);
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

  boost::math::chi_squared_distribution<double> chisq(5);

  if(isinf(lrt)){
    return "";
  }
  if(isnan(lrt)){
    return "";
  }

  if(lrt <= 0){
    return "";
  }
  
  double lp  = 1 - boost::math::cdf(chisq, lrt);

  string printdat =  seqid ;

  printdat.append("\t") ;
  printdat.append(boost::lexical_cast<string>(pos)) ;
  printdat.append("\t") ;
  printdat.append(boost::lexical_cast<string>(lrt)) ;
  printdat.append("\t") ;
  printdat.append(boost::lexical_cast<string>(lp)) ;
  printdat.append("\n");
  
  if(lp > 0.05){
    printdat = "";
  }

  return printdat; 
  
}

//------------------------------------------------------------
//------------------------------------------------------------

/// run the pileup

void pileup(int s, int j, int e,  map <string, int> target_info, vector<string> total){

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

  double  Nt = 1;
  double  St = 1;
  double  Nb = 1;
  double  Sb = 1;
  double  Na = 1;
  double  Sa = 1;

  std::vector <string> buffer;
  double onedump = 100000;
  onedump = (onedump / 1000000);

  while(mreader_thread.GetNextAlignment(al)){
    if(! al.IsFirstMate()){
      continue;
    }

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
   

    PileUp.proccess_alignment(al);
    string seqid = seqids[al.RefID].RefName;
    if(PileUp.currentStart() > PileUp.currentPos()){
      list<BamAlignment> dat =  PileUp.pileup();
      string ans = process_pileup(dat, target_info, PileUp.currentPos(), seqid, rolling_meant, rolling_meanb, rolling_meana);
      buffer.push_back(ans);
      if(buffer.size() > 100000){
	print_guard.lock();
	running_mb += onedump;
	cerr << "INFO:" << running_mb << "Mb finished" << endl; 
	printansvec("", buffer);	
	print_guard.unlock();
	buffer.clear();
      }
    }
  }
  print_guard.lock();
  running_mb += (lexical_cast<double>(buffer.size()) / 1000000);
  cerr << "INFO:" << running_mb << "Mb finished" << endl; 
  printansvec("", buffer);
  print_guard.unlock();
  buffer.clear();
}

//------------------------------------------------------------
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
  
  //boost::asio::io_service io_service;
  //auto_ptr<boost::asio::io_service::work> work(new boost::asio::io_service::work(io_service));

  
   
  asio::io_service io_service;
  auto_ptr<asio::io_service::work> work(new asio::io_service::work(io_service));

    //asio::io_service::work work (io_service);
  
  
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
      io_service.post(boost::bind(pileup, s, j, j+2000000, target_info, total)); 
    }
    io_service.post(boost::bind(pileup, s, j, boost::lexical_cast<int>(seq.Length), target_info, total)); 
    s += 1.0;
  }
  
  work.reset();

  threads.join_all();
  
  mreader.Close();

}

//------------------------------------------------------------
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
