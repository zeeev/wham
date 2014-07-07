#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <cmath>
#include <time.h>
#include "split.h"

#include "api/BamMultiReader.h"
#include "readPileUp.h"

using namespace std;
using namespace BamTools;

struct indvDat{
  int nReads;
  int notMapped;
  int mateMissing;
  int mateCrossChromosome;
};

struct global_opts {
  vector<string> targetBams;
  vector<string> backgroundBams;
  vector<string> region   ; 
  vector<string> regionDat;
} globalOpts;

static const char *optString ="ht:b:r:";

void printIndv(indvDat s, int t){
  cerr << "INDV: " << t << endl;
  cerr << " nReads      :" << s.nReads      << endl;
  cerr << " notMapped   :" << s.notMapped   << endl;
  cerr << " mateMissing :" << s.mateMissing << endl;
  cerr << " crossChrom  :" << s.mateCrossChromosome << endl;

  cerr << endl;
}

string join(vector<string> strings){

  string joined = "";

  for(vector<string>::iterator sit = strings.begin(); sit != strings.end(); sit++){
    joined = joined + " " + (*sit) + "\n";
  }

  return joined;

}

void printVersion(void){
  cerr << "Version 0.0.1 ; Zev Kronenberg; zev.kronenberg@gmail.com " << endl;
  cerr << endl;
}

void printHelp(void){
  cerr << "usage: wham -t <STRING> -b <STRING>" << endl << endl;
  cerr << "option: t <STRING> -- comma separated list of target bam files"           << endl ;
  cerr << "option: b <STRING> -- comma separated list of background bam files"       << endl ;
  cerr << "option: r <STRING> -- a genomic region in the format \"seqid:start-end\"" << endl ;
  cerr << endl;
  printVersion();
}

void parseOpts(int argc, char** argv){
  int opt = 0;

  opt = getopt(argc, argv, optString);

  while(opt != -1){
    switch(opt){
    case 't':
      {
	globalOpts.targetBams     = split(optarg, ",");
	cerr << "INFO: target bams:\n" << join(globalOpts.targetBams) ;
	break;
      }
    case 'b':
      {
	globalOpts.backgroundBams = split(optarg, ",");
	cerr << "INFO: background bams:\n" << join(globalOpts.backgroundBams) ;
	break;
      }
    case 'h':
      {
	printHelp();
	exit(1);
      }
    case 'r':
      {
	vector<string> tmp_region = split(optarg, ":");
	vector<string> start_end = split(tmp_region[1], "-");

	globalOpts.region.push_back(tmp_region[0]);
	globalOpts.region.push_back(start_end[0]);
	globalOpts.region.push_back(start_end[1]);
		
	cerr << "INFO: region set to: " <<   globalOpts.region[0] << ":" <<   globalOpts.region[1] << "-" <<  globalOpts.region[2] << endl;
	
	if(globalOpts.region.size() != 3){
	  cerr << "FATAL: incorrectly formatted region." << endl;
	  cerr << "FATAL: wham is now exiting."          << endl;
	  exit(1);
	}
	break;
      }
    case '?':
      {
	printHelp();
	exit(1);
      }  
  default:
    {
      cerr << "FATAL: Failure to parse command line options." << endl;
      cerr << "FATAL: Now exiting wham." << endl;
      printHelp();
      exit(1);
    }
    }
    opt = getopt( argc, argv, optString );
  }
  if( globalOpts.targetBams.empty() || globalOpts.backgroundBams.empty() ){
    cerr << "FATAL: Failure to specify target or background bam files." << endl;
    cerr << "FATAL: Now exiting wham." << endl;
    printHelp();
    exit(1);
  }
}

void prepBams(BamMultiReader & bamMreader, string group){

  int errorFlag    = 3;
  string errorMessage ;

  vector<string> files;

  if(group == "target"){
    files = globalOpts.targetBams;
  }
  if(group == "background"){
    files = globalOpts.backgroundBams;
  }
  if(files.empty()){
    cerr << "FATAL: no files ?" << endl;
    exit(1);
  }

  if(! bamMreader.Open(files)){
    errorMessage = bamMreader.GetErrorString();
    cerr << "FATAL: issue opening bams" << endl;
    exit(1);
  }
  if(! bamMreader.LocateIndexes()){
    errorMessage = bamMreader.GetErrorString();
    cerr << "FATAL: locating bam indicies"<< endl;
    exit(1); 
  }
}

double bound(double v){
  if(v <= 0.00001){
    return  0.00001;
  }
  if(v >= 0.99999){
    return 0.99999;
  }
  return v;
}


double mean(vector<double> & data){
  
  double sum = 0;
  double n   = 0;
  
  for(vector<double>::iterator it = data.begin(); it != data.end(); it++){
    sum += (*it);
    n += 1;
  }
  return sum / n;
}

double var(vector<double> & data, double mu){
  double variance = 0;

  for(vector<double>::iterator it = data.begin(); it != data.end(); it++){
    variance += pow((*it) - mu,2);
  }

  return variance / (data.size() - 1);
}


double logLbinomial(double x, double n, double p){

  double ans = lgamma(n+1)-lgamma(x+1)-lgamma(n-x+1) + x * log(p) + (n-x) * log(1-p);
  return ans;

}


double logDbeta(double alpha, double beta, double x){
  
  double ans = 0;

  ans += lgamma(alpha + beta) - ( lgamma(alpha) + lgamma(beta) );
  ans += log( pow(x, (alpha -1) ) ) + log( pow( (1-x) , (beta-1)));
  
  cerr << "alpha: " << alpha << "\t" << "beta: " << beta << "\t" << "x: " << x << "\t" << ans;
  cerr << endl;

  return ans;

}

double totalLL(vector<double> & data, double alpha, double beta){

  double total = 0;

  for(vector<double>::iterator it = data.begin(); it != data.end(); it++){
    total += logDbeta(alpha, beta, (*it));
  }
  return total;
}

double methodOfMoments(double mu, double var, double * aHat, double * bHat){

  double mui = 1 - mu;
  double right = ( (mu*mui / var) - 1 );

  (*aHat)  = mu *  right;
  (*bHat)  = mui * right; 
  
}


bool printGeno(map<string, indvDat> & dat, vector<string> & keys, string * toprint){
  
  bool nonRef = false;

  stringstream all;

  for(int i = 0; i < keys.size(); i++){
    
    stringstream gl;
    stringstream mn;
    stringstream genotypeField;

    string genotype = "0/0";
    
    if(dat[keys[i]].nReads < 2){
      genotype = "./.";
      gl << ".,.,.";
    }
    
    else{
      double max;

      double aa =  logLbinomial(dat[keys[i]].mateMissing, dat[keys[i]].nReads, 0.01);
      max = aa;
      double ab = logLbinomial(double(dat[keys[i]].mateMissing), double(dat[keys[i]].nReads), 0.50);
      if(ab > max){
	max = ab;
	genotype = "0/1";
	nonRef = true;
      }
      double bb = logLbinomial(double(dat[keys[i]].mateMissing), double(dat[keys[i]].nReads), 0.99);
      if(bb > max){
	max = bb;
	genotype = "1/1";
	nonRef = true;
      }

      gl << aa << "," << ab << "," << bb ;

    }
    mn << dat[keys[i]].mateMissing << "," << dat[keys[i]].nReads; 
    genotypeField << genotype << ":" << mn.str() << ":" << gl.str() ;
    
    
    if(i != (keys.size() - 1)){
      genotypeField << "\t";
    }
    all << genotypeField.str();
  }

  (*toprint) = all.str();

  return nonRef;
}

bool score(string seqid, long int pos, readPileUp & targDat, readPileUp & backDat, double * s, long int cp){
  
  map < string, indvDat > ti, bi;
  
  for(int t = 0; t < globalOpts.targetBams.size(); t++){
    indvDat i = {0};
    ti[globalOpts.targetBams[t]] = i;
  }
  for(int b = 0; b < globalOpts.backgroundBams.size(); b++){
    indvDat i = {0};
    bi[globalOpts.backgroundBams[b]] = i;
  }

  double totalTcount, totalBcount = 0;

  for(list<BamAlignment>::iterator tit = targDat.currentData.begin(); tit != targDat.currentData.end(); tit++){
    string fname = (*tit).Filename;
    totalTcount += 1;
    if((*tit).IsMapped()){
      ti[fname].nReads++;
    }
    ti[fname].notMapped   += (!(*tit).IsMapped());
    ti[fname].mateMissing += (!(*tit).IsMateMapped());
    ti[fname].mateCrossChromosome += ((*tit).RefID == (*tit).MateRefID);
  }

  for(list<BamAlignment>::iterator bit = backDat.currentData.begin(); bit != backDat.currentData.end(); bit++){
    string fname = (*bit).Filename;
    totalBcount += 1;
    if((*bit).IsMapped()){
      bi[fname].nReads++;
    }
    bi[fname].notMapped   += (*bit).IsMapped();
    bi[fname].mateMissing += (!(*bit).IsMateMapped());
    bi[fname].mateCrossChromosome += ((*bit).RefID == (*bit).MateRefID);
  }
  
  if(totalTcount < 10 || totalBcount < 10){
    return true;
  }

  vector<double> tiFrq, biFrq, totFrq;


  for(map<string, indvDat>::iterator tiIt = ti.begin(); tiIt != ti.end(); tiIt++){
    if((tiIt->second).nReads < 2){
      continue;
    }
    else{
      double frq = double(tiIt->second.mateMissing) / double(tiIt->second.nReads);
      tiFrq.push_back(frq);
      totFrq.push_back(frq);
    }
  }

  for(map<string, indvDat>::iterator biIt = bi.begin(); biIt != bi.end(); biIt++){
    if((biIt->second).nReads < 2){
      continue;
    }
    else{
      double frq = double(biIt->second.mateMissing) / double(biIt->second.nReads);
      biFrq.push_back(frq);
      totFrq.push_back(frq);
    }
  }

  if(totFrq.size() < 2 || tiFrq.size() < 2 || biFrq.size() < 2 ){
    return true;
  }

  string genosT, genosB;
  
  bool altT = printGeno(ti, globalOpts.targetBams,     &genosT);
  bool altB = printGeno(bi, globalOpts.backgroundBams, &genosB);
  
  if(altT && altB){
    cout << seqid << "\t" << pos << "\t" << genosT << "\t" << genosB << endl;
  }

  double tiMean  =  mean(tiFrq);
  double biMean  =  mean(biFrq);
  double totMean =  mean(totFrq);

  (*s) = tiMean - biMean;

  if((tiMean + biMean) < 0.2){
    return true;
  }

  
//  double tiVar  = var(tiFrq, tiMean);
//  double biVar  = var(biFrq, biMean);
//  double totVar = var(totFrq, totMean);
//
//  double tiAhat, tiBhat, biAhat, biBhat, totAhat, totBhat;
//  
//  methodOfMoments(tiMean, tiVar, &tiAhat, &tiBhat);
//  methodOfMoments(biMean, biVar, &biAhat, &biBhat);
//  methodOfMoments(totMean, totVar, &totAhat, &totBhat);
//
//  cerr << endl;
//  cerr << tiAhat << "\t" << tiBhat << endl;
//
//  double alt  = totalLL(tiFrq, tiAhat, tiBhat)   + totalLL(biFrq, biAhat, biBhat);
//  double null = totalLL(tiFrq, totAhat, totBhat) + totalLL(biFrq, totAhat, totBhat); 

  //  cout << cp << "\t" <<  tiMean << "\t" << biMean << "\t" << totMean << "\t" << alt/null << endl;

  return true;
}

bool runRegion(RefData region, int seqidIndex, vector< RefData > seqNames){

  BamMultiReader targetReader, backgroundReader;
  
  prepBams(targetReader, "target");
  prepBams(backgroundReader, "background");
  
  int setRegionErrorFlag = 0;

  setRegionErrorFlag += targetReader.SetRegion(seqidIndex, 0, seqidIndex, region.RefLength);
  setRegionErrorFlag += backgroundReader.SetRegion(seqidIndex, 0, seqidIndex, region.RefLength);
  
  if(setRegionErrorFlag != 2){
    return false;
  }
    
  BamAlignment alt, alb;
  
  readPileUp targetPileUp, backgroundPileUp;

  if(! targetReader.GetNextAlignment(alt)){
    targetReader.Close();
    backgroundReader.Close();
    return false;
  }
  if(! backgroundReader.GetNextAlignment(alb)){
    targetReader.Close();
    backgroundReader.Close();
    return false;
  }

  long int currentPos = -1;

  targetPileUp.processAlignment(alt,     currentPos);
  backgroundPileUp.processAlignment(alb, currentPos);

  bool stillReads = true;
  
  while(stillReads){
    while(currentPos > targetPileUp.currentStart() || currentPos > backgroundPileUp.currentStart()){
      
      bool getTarget, getBackground;

      if(currentPos > targetPileUp.currentStart()){
	getTarget = targetReader.GetNextAlignment(alt);
	targetPileUp.processAlignment(alt, currentPos);
      }
      if(currentPos > backgroundPileUp.currentStart()){
	getBackground = backgroundReader.GetNextAlignment(alb);
	backgroundPileUp.processAlignment(alb, currentPos);
      }
      if(getTarget == false && getBackground == false){
	stillReads = false;
	break;
      }
    }
   
    targetPileUp.purgePast();
    backgroundPileUp.purgePast();


    double s = 0;

    //    cout << seqNames[seqidIndex].RefName << "\t" << currentPos << "\t" << s << "\t" << "GT:GL" << "\t";

    if(! score(seqNames[seqidIndex].RefName, currentPos, targetPileUp, backgroundPileUp, &s, currentPos )){
      cerr << "FATAL: problem during scoring" << endl;
      cerr << "FATAL: wham exiting"           << endl;
      exit(1);
    }

    currentPos += 50;


    
    



  }


  targetReader.Close();
  backgroundReader.Close();
  return true;
}


int main(int argc, char** argv) {

  srand((unsigned)time(NULL));

  parseOpts(argc, argv);

  BamMultiReader targetReader, backgroundReader;
  
  prepBams(targetReader, "target");
  prepBams(backgroundReader, "background");
  
  RefVector sequences = targetReader.GetReferenceData();

  targetReader.Close();
  backgroundReader.Close();

  int seqidIndex = 0;
    
  for(vector< RefData >::iterator sit = sequences.begin(); sit != sequences.end(); sit++){
    //    cerr << (*sit).RefName << endl;
    
    if(globalOpts.region.size() == 3){
      if((*sit).RefName == globalOpts.region[0]){
	if(! runRegion((*sit), seqidIndex, sequences)){
	  cerr << "FATAL: region failed to run properly: " << (*sit).RefName << endl;
	  cerr << "FATAL: Wham exiting" << endl;
	}
      }
    }
    else{
      if(! runRegion((*sit), seqidIndex, sequences)){
	cerr << "FATAL: region failed to run properly: " << (*sit).RefName << endl;
	cerr << "FATAL: Wham exiting" << endl;
      }
    }
    seqidIndex += 1;
  }

}
