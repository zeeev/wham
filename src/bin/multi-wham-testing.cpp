#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <cmath>
#include <time.h>
#include "split.h"

#include "api/BamMultiReader.h"
#include "readPileUp.h"

#define MATHLIB_STANDALONE

using namespace std;
using namespace BamTools;

struct indvDat{
  int nReads;
  int nAboveAvg;
  int notMapped;
  int mateMissing;
  int mateCrossChromosome;
  vector<BamAlignment> data;
};

struct insertDat{
  map<string, double> mus; // mean of insert length for each indvdual across 1e6 reads
  map<string, double> sds;  // standard deviation
  map<string, double> lq ;  // 25% of data
  map<string, double> up ;  // 75% of the data
} insertDists;

struct global_opts {
  vector<string> targetBams    ;
  vector<string> backgroundBams;
  vector<string> all           ;
  string         seqid         ;
  vector<int>    region        ; 
} globalOpts;

static const char *optString ="ht:b:r:";

void initIndv(indvDat * s){
    s->nReads              = 0;
    s->nAboveAvg           = 0;
    s->notMapped           = 0;
    s->mateMissing         = 0;
    s->mateCrossChromosome = 0;
    s->data.clear();
}

void printIndv(indvDat * s, int t){
  cout << "#INDV: " << t << endl;
  cout << "# nReads      :" << s->nReads      << endl;
  cout << "# notMapped   :" << s->notMapped   << endl;
  cout << "# mateMissing :" << s->mateMissing << endl;
  cout << "# crossChrom  :" << s->mateCrossChromosome << endl;
  cout << "# n 1sd above :" << s->nAboveAvg           << endl;

  double sum = 0;

  string readNames = "# read names: ";

  for(int i = 0; i <  s->data.size(); i++){
    sum += double (s->data[i].MapQuality);
    readNames.append(s->data[i].Name);
    readNames.append(" ");
  }

  cout << "# average mapping q : " << sum/double(s->nReads) << endl;
  cout << readNames << endl;
  cout << endl;
}

void printHeader(void){
  cout << "##fileformat=VCFv4.1"                                                                                                                  << endl;
  cout << "#INFO=<LRT,Number=1,type=Float,Description=\"Likelihood Ratio Test Statistic\">"                                                       << endl;
  cout << "#INFO=<EAF,Number=3,type=Float,Description=\"Allele frequency approximation based on mapping quality of: target,background,combined\">" << endl;
  cout << "#INFO=<AF,Number=3,type=Float,Description=\"Allele frequency of: target,background,combined\">" << endl;
  cout << "#INFO=<NALT,Number=2,type=Int,Description=\"Number of alternative pseuod alleles for target and background \">" << endl;
  cout << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Pseudo genotype\">"                                                                 << endl;
  cout << "##FORMAT=<ID=MM,Number=1,Type=Int,Description=\"Number of Missing Mates\">"                                                            << endl;
  cout << "##FORMAT=<ID=DP,Number=1,Type=Int,Description=\"Number of reads with mapping quality greater than 0\">"                                << endl;
  cout << "##FORMAT=<GL,Number=A,type=Float,Desciption=\"Genotype likelihood \">"                                                                 << endl;
  cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" << "\t";

  for(int b = 0; b < globalOpts.all.size(); b++){
    cout << globalOpts.all[b] ;
    if(b < globalOpts.all.size() - 1){
      cout << "\t";
    }
  }
  cout << endl;  
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
  cerr << "usage: WHAM-BAM -t <STRING> -b <STRING> -r <STRING>" << endl << endl;
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

	globalOpts.seqid = tmp_region[0];
	globalOpts.region.push_back(atoi(start_end[0].c_str()));
	globalOpts.region.push_back(atoi(start_end[1].c_str()));
		
	cerr << "INFO: region set to: " <<   globalOpts.seqid << ":" <<   globalOpts.region[0] << "-" <<  globalOpts.region[1] << endl;
	
	if(globalOpts.region.size() != 2){
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

bool grabInsertLengths(string file){

  BamReader bamR;
  BamAlignment al;

  vector<double> alIns;

  bamR.Open(file);

  int i = 1;

  while(i < 100000 && bamR.GetNextAlignment(al)){
    if(al.IsFirstMate() && al.IsMapped() && al.IsMateMapped() && abs(double(al.InsertSize)) < 10000){
      i += 1;
      alIns.push_back(abs(double(al.InsertSize)));

    }
  }

  bamR.Close();

  double mu       = mean(alIns        );
  double variance = var(alIns, mu     );
  double sd       = sqrt(variance     );

  insertDists.mus[file] = mu; 
  insertDists.sds[file] = sd; 

  cerr << "INFO: mean insert length, number of reads, file : " 
       << insertDists.mus[file] << ", " 
       << insertDists.sds[file] << ", "
       << i  << ", " 
       << file << endl; 

  return true;
  
}



bool getInsertDists(void){

  bool flag = true;

  for(int i = 0; i < globalOpts.all.size(); i++){
    flag = grabInsertLengths(globalOpts.all[i]);
  }
  
  return true;
  
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

double unphred(double p){
  return pow(10, (-p/10));
}


bool processGenotype(indvDat * idat, vector<double> & gls, string * geno, double * nr, double * na, double * ga, double * gb, double * aa, double * ab, double * gc){
  
  bool alt = false;
  
  double aal = 0;
  double abl = 0;
  double bbl = 0;
  (*geno) = "./.";
  
  double nreads = idat->data.size();

  if(nreads < 3){
    gls.push_back(-255.0);
    gls.push_back(-255.0);
    gls.push_back(-255.0);

    return false;
    
  }

  double nref = 0;
  double nalt = 0;

  stringstream m;


  bool odd = false;
  if( idat->mateMissing > 2 || idat->nAboveAvg > 2 ){
    odd = true;
  }

  for(vector<BamAlignment>::iterator it = idat->data.begin(); it != idat->data.end(); it++ ){
    
    double mappingP = unphred((*it).MapQuality);

    bool sameStrand = false;
    
    if((*it).IsMateMapped() && ( (*it).IsReverseStrand() && (*it).IsMateReverseStrand() ) || ( !(*it).IsReverseStrand() && !(*it).IsMateReverseStrand() ) || sameStrand ){
      sameStrand = true;
    }
    
    m << "\t" << mappingP ;
    
    double iDiff = abs ( abs ( double ( (*it).InsertSize ) - insertDists.mus[(*it).Filename] ));
    
    
    if( ( odd &&  (iDiff > 3 * insertDists.sds[(*it).Filename] ))  || ! (*it).IsMateMapped() || sameStrand ){ 
      nalt += 1;
      aal += log((2-2) * (1-mappingP) + (2*mappingP)) ;
      abl += log((2-1) * (1-mappingP) + (1*mappingP)) ;
      bbl += log((2-0) * (1-mappingP) + (0*mappingP)) ;
    }
    else{
      nref += 1;
      aal += log((2 - 2)*mappingP + (2*(1-mappingP)));
      abl += log((2 - 1)*mappingP + (1*(1-mappingP)));
      bbl += log((2 - 0)*mappingP + (0*(1-mappingP)));
    } 
  }
  
  aal = aal - log(pow(2,nreads));
  abl = abl - log(pow(2,nreads));
  bbl = bbl - log(pow(2,nreads));

  if(nref == 0){
    aal = -255.0;
    abl = -255.0;
  }
  if(nalt == 0){
    abl = -255.0;
    bbl = -255.0;
  }

  double max = aal;
  (*geno) = "0/0";

  if(abl > max){
    (*geno) = "0/1";
    max = abl;
    alt = true;
    *nr+=1;
    *na+=1;
  }
  if(bbl > max){
    (*geno) = "1/1";
    *na += 2 ;
    alt = true;
  }
  if((*geno) == "0/0"){
    *nr += 2;
  }
  
  *ga += 2* exp(aal);
  *aa += 2* exp(aal);
  *ga += exp(abl);
  *aa += exp(abl);
  
  *gb += 2* exp(bbl);
  *ab += 2* exp(bbl);
  *gb += exp(abl);
  *ab += exp(abl);
  
  gls.push_back(aal);
  gls.push_back(abl);
  gls.push_back(bbl);

  *gc += 1;

  //  cerr << (*geno) << "\t" << aal << "\t" << abl << "\t" << bbl << "\t" << nref << "\t" << nalt << m.str() << endl;

  return alt;

}

void pseudoCounts(vector< vector <double> > &dat, double *alpha, double * beta){
  
  for(int i = 0; i < dat.size(); i++){

    (*alpha) += 2 * exp(dat[i][0]);
    (*alpha) +=     exp(dat[i][1]);
    (*beta ) +=     exp(dat[i][1]);
    (*beta ) += 2 * exp(dat[i][2]);

  }
}

bool score(string seqid, long int pos, readPileUp & targDat, readPileUp & backDat, double * s, long int cp){

  map < string, indvDat*> ti, bi;
  
  for(int t = 0; t < globalOpts.targetBams.size(); t++){
    indvDat * i;
    i = new indvDat;
    initIndv(i);
    ti[globalOpts.targetBams[t]] = i;
  }
  for(int b = 0; b < globalOpts.backgroundBams.size(); b++){
    indvDat * i;
    i = new indvDat;
    initIndv(i);
    bi[globalOpts.backgroundBams[b]] = i;
  }

  int td = 0;
  int bd = 0;
  
  for(list<BamAlignment>::iterator tit = targDat.currentData.begin(); tit != targDat.currentData.end(); tit++){
    string fname = (*tit).Filename;
    if((*tit).IsMapped() && (pos > ( (*tit).Position) ) && ( (*tit).MapQuality > 0 )){
            
      double insdiff = double ( (*tit).InsertSize) - insertDists.mus[fname];

      if(insdiff > 3.5 * insertDists.sds[fname]){
	ti[fname]->nAboveAvg++;
      }

      ti[fname]->nReads++; 
      ti[fname]->mateMissing += (!(*tit).IsMateMapped());
      ti[fname]->data.push_back(*tit);
      
      td += 1;
    }
  }

  for(list<BamAlignment>::iterator bit = backDat.currentData.begin(); bit != backDat.currentData.end(); bit++){
    string fname = (*bit).Filename;

    double insdiff = double ( (*bit).InsertSize) - insertDists.mus[fname];

    if(insdiff > 3.5 * insertDists.sds[fname]){
      bi[fname]->nAboveAvg++;
    }
   
    if((*bit).IsMapped() && (pos > ( (*bit).Position) ) && ( (*bit).MapQuality > 0 )){
      bi[fname]->nReads++;
      bi[fname]->mateMissing += (!(*bit).IsMateMapped());
      bi[fname]->data.push_back(*bit);
   
      bd += 1;
      
    }
  }


  if(td < 5 || bd < 5){
    for(vector<string>::iterator all = globalOpts.all.begin(); all != globalOpts.all.end(); all++ ){
      delete ti[*all];
      delete bi[*all];
    }
    return true;
  
  }

  vector< vector < double > > targetGls, backgroundGls, totalGls;
  vector<string> genotypes;
  vector<int> depth;
  vector<int> nMissingMate;

  int nAlt = 0;

  double targetRef = 0;
  double targetAlt = 0;
  
  double backgroundRef = 0;
  double backgroundAlt = 0;

  double ta  = 0.0001;
  double tb  = 0.0001;
  double ba  = 0.0001;
  double bb  = 0.0001;
  double aa  = 0.0001;
  double ab  = 0.0001;
  double tgc = 0;
  double bgc = 0;
  
  for(int t = 0; t < globalOpts.targetBams.size(); t++){
    string genotype; 
    vector< double > Targ;
    nMissingMate.push_back(ti[globalOpts.targetBams[t]]->mateMissing);
    depth.push_back(ti[globalOpts.targetBams[t]]->nReads);
    nAlt += processGenotype(ti[globalOpts.targetBams[t]], Targ, &genotype, &targetRef, &targetAlt, &ta, &tb, &aa, &ab, &tgc);
    genotypes.push_back(genotype);
    
    targetGls.push_back(Targ);
    totalGls.push_back(Targ);
    
  }
  for(int b = 0; b < globalOpts.backgroundBams.size(); b++){
    string genotype; 
    vector< double > Back;
    nMissingMate.push_back( bi[globalOpts.backgroundBams[b]]->mateMissing);
    depth.push_back(bi[globalOpts.backgroundBams[b]]->nReads);
    nAlt += processGenotype(bi[globalOpts.backgroundBams[b]], Back, &genotype, &backgroundRef, &backgroundAlt, &ba, &bb, &aa, &ab, &bgc);
    genotypes.push_back(genotype);

    backgroundGls.push_back(Back);
    totalGls.push_back(Back);

  }
  
  double trueTaf = targetAlt / (targetAlt + targetRef) ;
  double trueBaf = backgroundAlt / (backgroundAlt + backgroundRef);

  if(nAlt < 2 ){
    for(vector<string>::iterator all = globalOpts.all.begin(); all != globalOpts.all.end(); all++ ){
      delete ti[*all];
      delete bi[*all];
    }
    return true;
  }



//  pseudoCounts(targetGls, &ta, &tb);
//  pseudoCounts(backgroundGls, &ba, &bb);
//  pseudoCounts(totalGls, &aa, &ab);
  
  double taf = bound(tb / (ta + tb));
  double baf = bound(bb / (ba + bb));
  double aaf = bound((tb + bb) / (ta + tb + ba + bb));


  double tafG = bound(targetAlt / (targetAlt + targetRef));
  double bafG = bound(backgroundAlt / (backgroundAlt + backgroundRef));
  double aafG = bound((targetAlt + backgroundAlt) / (targetAlt + targetRef + backgroundAlt + backgroundRef));

  double alt  = logLbinomial(tb, (tb+ta), taf) + logLbinomial(bb, (bb+ba), baf);
  double null = logLbinomial(tb, (tb+ta), aaf) + logLbinomial(bb, (bb+ba), aaf);

  double lrt = 2 * (alt - null);

  double altG  = logLbinomial(targetAlt, (targetAlt+targetRef), tafG) + logLbinomial(backgroundAlt, (backgroundAlt+backgroundRef), bafG);
  double nullG = logLbinomial(targetAlt, (targetAlt+targetRef), aafG) + logLbinomial(backgroundAlt, (backgroundAlt+backgroundRef), aafG);

  double lrtG = 2 * (altG - nullG);

  
  if(lrtG < lrt){
    lrt = lrtG;
  }

  if(lrt <= 0){
    for(vector<string>::iterator all = globalOpts.all.begin(); all != globalOpts.all.end(); all++ ){
      delete ti[*all];
      delete bi[*all];
    }
    return true;
  }

  cout << seqid   << "\t" ;       // CHROM
  cout << pos     << "\t" ;       // POS
  cout << "."     << "\t" ;       // ID
  cout << "NA"    << "\t" ;       // REF
  cout << "SV"    << "\t" ;       // ALT
  cout << "."     << "\t" ;       // QUAL
  cout << "."     << "\t" ;       // FILTER
  cout << "LRT="  << lrt  ;        // INFO
  cout << ";EAF="  << taf ;
  cout << ","     << baf  ;
  cout << ","     << aaf  ;
  cout << ";AF="  << trueTaf << "," << trueBaf ;
  cout << ";NALT=" << targetAlt << "," << backgroundAlt;
  cout << ";GC="   << tgc       << "," << tgc << "\t";
  cout << "GT:GL:MM:DP" << "\t" ;
      
  int index = 0;

  for(vector<string>::iterator it = genotypes.begin(); it != genotypes.end(); it++){
    
    stringstream ss ; 
    
    cout << (*it);

    ss  << ":"  <<  nMissingMate[index] << ":" << depth[index] << ":" << totalGls[index][0] << "," << totalGls[index][1] << "," << totalGls[index][2];
    cout << ss.str();

    if(it + 1 != genotypes.end()){
      cout << "\t";
    }
    index += 1;
  }
  cout << endl;
  
//  if(lrt > 3){
//    printIndv(ti[globalOpts.targetBams[0]],0 );
//    printIndv(ti[globalOpts.targetBams[1]],1 );
//    printIndv(ti[globalOpts.targetBams[2]],2 );
//    printIndv(ti[globalOpts.targetBams[3]],3 );
//    printIndv(ti[globalOpts.targetBams[4]],4 );
//    printIndv(ti[globalOpts.targetBams[5]],5 );
//    printIndv(ti[globalOpts.targetBams[6]],6 );
//    printIndv(ti[globalOpts.targetBams[7]],7 );
//    printIndv(ti[globalOpts.targetBams[8]],8 );
//    printIndv(ti[globalOpts.targetBams[9]],9 );
//  }

  for(vector<string>::iterator all = globalOpts.all.begin(); all != globalOpts.all.end(); all++ ){
    delete ti[*all];
    delete bi[*all];
  }
  
  return true;
}

  



bool runRegion(int seqidIndex, int start, int end, vector< RefData > seqNames){

  BamMultiReader targetReader, backgroundReader;
  
  prepBams(targetReader, "target");
  prepBams(backgroundReader, "background");
  
  int setRegionErrorFlag = 0;

  setRegionErrorFlag += targetReader.SetRegion(seqidIndex, start, seqidIndex, end);
  setRegionErrorFlag += backgroundReader.SetRegion(seqidIndex, start, seqidIndex, end);
  
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


  bool getTarget     = true;
  bool getBackground = true;  

  while(1){
    
    while(currentPos > targetPileUp.CurrentStart && getTarget){
      getTarget = targetReader.GetNextAlignment(alt);
      targetPileUp.processAlignment(alt, currentPos);
    }
    while(currentPos > backgroundPileUp.CurrentStart && getBackground){
      getBackground = backgroundReader.GetNextAlignment(alb);
      backgroundPileUp.processAlignment(alb, currentPos);
    }
    
    if(getTarget == false && getBackground == false){
      break;
    }
    
    targetPileUp.purgePast();
    backgroundPileUp.purgePast();

    double s = 0;

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

  globalOpts.all.reserve(globalOpts.targetBams.size()  + globalOpts.backgroundBams.size());
  globalOpts.all.insert( globalOpts.all.end(), globalOpts.targetBams.begin(), globalOpts.targetBams.end() );
  globalOpts.all.insert( globalOpts.all.end(), globalOpts.backgroundBams.begin(), globalOpts.backgroundBams.end() );

  BamMultiReader targetReader, backgroundReader;
  
  prepBams(targetReader, "target");
  prepBams(backgroundReader, "background");

  if(!getInsertDists()){
    cerr << "FATAL: " << "problem while generating insert lengths dists" << endl;
    exit(1);
  }
  
  RefVector sequences = targetReader.GetReferenceData();

  targetReader.Close();
  backgroundReader.Close();

  printHeader();

  int seqidIndex = 0;

  //#pragma omp parallel
  {
    for(vector< RefData >::iterator sit = sequences.begin(); sit != sequences.end(); sit++){
      //    cerr << (*sit).RefName << endl;
      
      if(globalOpts.region.size() == 2){
	if((*sit).RefName == globalOpts.seqid){
	  //	cerr << "runing region" << endl;
	  if(! runRegion(seqidIndex, globalOpts.region[0], globalOpts.region[1], sequences)){
	    cerr << "WARNING: region failed to run properly: " << (*sit).RefName << endl;
	  }
	  //	cerr << "ran region" << endl;
	}
      }
      else{
	//	cerr << "runing region" << endl;
	if(! runRegion(seqidIndex, 0, (*sit).RefLength, sequences)){
	  cerr << "WARNING: region failed to run properly: " << (*sit).RefName << endl;

	}
	//      cerr << "ran region" << endl;
      }
      seqidIndex += 1;
    }
  }

  
  cerr << "INFO: WHAM-BAM finished normally." << endl;
  return 0;

}
