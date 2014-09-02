#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <cmath>
#include <time.h>
#include <algorithm>
#include "split.h"
#include "JenksBreaks.h"

#include <omp.h>
#include "api/BamMultiReader.h"
#include "readPileUp.h"

using namespace std;
using namespace BamTools;

struct regionDat{
  int seqidIndex ;
  int start      ;
  int end        ;
};

struct indvDat{
  bool   support;
  int    nReads;
  int    mappedPairs;
  int    nAboveAvg;
  int    notMapped;
  int    mateMissing;
  int    sameStrand;
  int    mateCrossChromosome;
  int    maxLength ;
  double insertSum ;
  double insertMean;
  double lengthSum;
  double clipped; 
  vector< BamAlignment > alignments;
  vector<double> inserts;
  vector<double> hInserts;
  map<string, int> badFlag;
  vector<int> MapQ;
  map<int, vector<string> > cluster;
};

struct insertDat{
  map<string, double> mus; // mean of insert length for each indvdual across 1e6 reads
  map<string, double> sds;  // standard deviation
  map<string, double> lq ;  // 25% of data
  map<string, double> up ;  // 75% of the data
} insertDists;

struct global_opts {
  string         targetBam     ;
  int            nthreads      ;
  string         seqid         ;
  vector<int>    region        ; 
} globalOpts;

static const char *optString ="ht:b:r:x:";

omp_lock_t lock;

void initIndv(indvDat * s){
  s->support             = false;
  s->nReads              = 0;
  s->nAboveAvg           = 0;
  s->notMapped           = 0;
  s->mappedPairs         = 0;
  s->mateMissing         = 0;
  s->insertSum           = 0;
  s->sameStrand          = 0; 
  s->maxLength           = 0;
  s->mateCrossChromosome = 0;
  s->lengthSum           = 0;
  s->clipped             = 0;
  s->alignments.clear();
  s->inserts.clear();
  s->hInserts.clear();
  s->badFlag.clear();
  s->MapQ.clear();
}

void printIndv(indvDat * s, int t){
  cout << "#INDV: " << t << endl;
  cout << "# nReads      :" << s->nReads      << endl;
  cout << "# notMapped   :" << s->notMapped   << endl;
  cout << "# mateMissing :" << s->mateMissing << endl;
  cout << "# crossChrom  :" << s->mateCrossChromosome << endl;
  cout << "# n 1sd above :" << s->nAboveAvg           << endl;

  double sum = 0;

  for(int i = 0; i <  s->MapQ.size(); i++){
    sum += double (s->MapQ[i]);
  }

  cout << "# average mapping q : " << sum/double(s->nReads) << endl;
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

    cout << globalOpts.targetBam ;
    cout << endl;  
}

void printVersion(void){
  cerr << "Version 0.0.1 ; Zev Kronenberg; zev.kronenberg@gmail.com " << endl;
  cerr << endl;
}

void printHelp(void){
  cerr << "usage: WHAM-BAM -t <STRING> -b <STRING> -r <STRING>" << endl << endl;
  cerr << "required: t <STRING> -- comma separated list of target bam files"           << endl ;
  cerr << "required: b <STRING> -- comma separated list of background bam files"       << endl ;
  cerr << "option  : r <STRING> -- a genomic region in the format \"seqid:start-end\"" << endl ;
  cerr << "option  : x <INT>    -- set the number of threads, otherwise max          " << endl ; 
  cerr << endl;
  printVersion();
}

void parseOpts(int argc, char** argv){
  int opt = 0;

  opt = getopt(argc, argv, optString);

  while(opt != -1){
    switch(opt){
    case 'x':
      {
	globalOpts.nthreads = atoi(((string)optarg).c_str());
	cerr << "INFO: OpenMP will roughly use " << globalOpts.nthreads << " threads" << endl;
	break;
      }
    case 't':
      {
	globalOpts.targetBam     = optarg;
	cerr << "INFO: target bam:\n" << globalOpts.targetBam ;
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

  double clipped  = 0;
  double naligned = 0;

  bamR.Open(file);

  int i = 1;
  while(i < 100000 && bamR.GetNextAlignment(al) && abs(double(al.InsertSize)) < 10000){
    if(al.IsFirstMate() && al.IsMapped() && al.IsMateMapped()){
      i += 1;
      alIns.push_back(abs(double(al.InsertSize)));
    }
    if(al.IsMapped()){
      naligned += 1;
      vector< CigarOp > cd = al.CigarData;

      if(cd.back().Type == 'H' || cd.back().Type == 'S'){
	clipped += double(cd.back().Length) / double (al.Length);
      }
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
  cerr << "INFO: fraction of total length clipped, file : " << clipped / naligned 
       << ", " 
       << file << endl;


  return true;
  
}


bool getInsertDists(void){

    grabInsertLengths(globalOpts.targetBam);
  
    return true;
    
}

void prepBams(BamReader & bamReader, string group){

  int errorFlag    = 3;
  string errorMessage ;
 
  bool attempt = true;
  int  tried   = 0   ;
  
  while(attempt && tried < 500){
    
    if( bamReader.Open(globalOpts.targetBam) && bamReader.LocateIndex() ){
      attempt = false;
    }
    else{
      tried += 1;
    }
  }
  if(attempt == true){
    cerr << "FATAL: unable to open BAMs or indices after: " << tried << " attempts"  << endl;  
    cerr << "Bamtools error message:\n" <<  bamReader.GetErrorString()  << endl;
    cerr << "INFO : try using less CPUs in the -x option" << endl;
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


int processGenotype(indvDat * idat, vector<double> & gls, string * geno, double * nr, double * na, double * ga, double * gb, double * aa, double * ab, double * gc, insertDat & localDists, int * ss, int * insert ){

  int alt = 0;
  
  double aal = 0;
  double abl = 0;
  double bbl = 0;
  (*geno) = "./.";
  
  int nreads = idat->badFlag.size();

  if(nreads < 3){
    gls.push_back(-255.0);
    gls.push_back(-255.0);
    gls.push_back(-255.0);

    return alt;
    
  }

  double nref = 0;
  double nalt = 0;

  int ri = 0;

  for(map< string, int >::iterator rit = idat->badFlag.begin(); rit != idat->badFlag.end(); rit++){
    
    double mappingP = unphred(idat->MapQ[ri]);

    if( idat->badFlag[rit->first] == 1 && idat->support == true){ 
      nalt += 1;
      alt  += 1;
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
    ri++;
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

  //   cerr << (*geno) << "\t" << aal << "\t" << abl << "\t" << bbl << "\t" << nref << "\t"  << endl;

  return nalt;

}

void pseudoCounts(vector< vector <double> > &dat, double *alpha, double * beta){
  
  for(int i = 0; i < dat.size(); i++){

    (*alpha) += 2 * exp(dat[i][0]);
    (*alpha) +=     exp(dat[i][1]);
    (*beta ) +=     exp(dat[i][1]);
    (*beta ) += 2 * exp(dat[i][2]);

  }
}

double permute(vector< vector<double> > & dat, int tsize, int bsize, double score){
  
  double nrep = 0;
  double nsuc = 0;
  
  for(int i = 0; i < 1000000; i++){
    nrep += 1;
    
    if(nsuc > 1){
      break;
    }

    double tr = 0.000001;
    double ta = 0.000001;
    double br = 0.000001;
    double ba = 0.000001;
    double ar = 0.000001;
    double aa = 0.000001;

    random_shuffle(dat.begin(), dat.end());

    //    cerr << "TEST" << endl;

    int index = 0;
    
    for(; index < tsize; index++){
      tr += 2 * exp(dat[index][0]);
      tr +=     exp(dat[index][1]);
      ta +=     exp(dat[index][1]);
      ta += 2 * exp(dat[index][2]);
    }
    for(; index < dat.size(); index++){
      br += 2 * exp(dat[index][0]);
      br +=     exp(dat[index][1]);
      ba +=     exp(dat[index][1]);
      ba += 2 * exp(dat[index][2]);
    }
    
    ar += tr + br;
    aa += ta + ba;

    double taf = bound(ta / (ta+tr));
    double baf = bound(ba / (ba+br));
    
    double aaf = bound((ta+ba)/ (ta+ba+tr+br));
    
    double alt  = logLbinomial(ta, (tr+ta), taf) + logLbinomial(ba, (br+ba), baf);
    double null = logLbinomial(ta, (tr+ta), aaf) + logLbinomial(ba, (br+ba), aaf);

    double lrt = 2 * (alt - null);
    
    if(lrt > score){
      nsuc += 1;
    }

  }

  return nsuc / nrep;

}
    
bool loadIndv(map<string, indvDat*> & ti, readPileUp & pileup, global_opts localOpts, insertDat & localDists, long int pos){

  for(list<BamAlignment>::iterator r = pileup.currentData.begin(); r != pileup.currentData.end(); r++){

    if( (*r).IsMapped() && (pos > ( (*r).Position) ) && ( (*r).MapQuality > 0 )){

      string fname = (*r).Filename;

      int bad = 0;      

      ti[fname]->alignments.push_back(*r);

      ti[fname]->nReads += 1;

      vector< CigarOp > cd = (*r).CigarData;

      ti[fname]->lengthSum += (*r).Length;
      
      if(cd.back().Type == 'H' || cd.back().Type == 'S'){
	ti[fname]->clipped += cd.back().Length;
	int location = (*r).GetEndPosition();
	ti[fname]->cluster[location].push_back((*r).Name);
      }
      
      if(!(*r).IsMateMapped()){
	ti[fname]->mateMissing  += (!(*r).IsMateMapped());
      	bad = 1;
      }

      if((*r).IsMapped() && (*r).IsMateMapped()){
	
	ti[fname]->insertSum    += abs(double((*r).InsertSize));	
	ti[fname]->mappedPairs  += 1;

	if(( (*r).IsReverseStrand() && (*r).IsMateReverseStrand() ) || ( !(*r).IsReverseStrand() && !(*r).IsMateReverseStrand() )){
	  bad = 1;
	  ti[fname]->sameStrand += 1;
	}
	
	double ilength = abs ( double ( (*r).InsertSize ));

	double iDiff = abs ( ilength - localDists.mus[(*r).Filename] );

	ti[fname]->inserts.push_back(ilength);

	if(iDiff > (3.0 * insertDists.sds[(*r).Filename]) ){
	  bad = 1;
	  ti[fname]->nAboveAvg += 1;
	  ti[fname]->hInserts.push_back(ilength);
	}
      }
      
      ti[fname]->badFlag[(*r).Name] = bad;
      ti[fname]->MapQ.push_back((*r).MapQuality);      
    }
  }


  // looping over indviduals
  for(map < string, indvDat*>::iterator indvs = ti.begin(); indvs != ti.end(); indvs++){
    // looping over clusters
    for( map< int, vector<string> >::iterator ci = ti[indvs->first]->cluster.begin(); ci != ti[indvs->first]->cluster.end(); ci++){
      if(ci->second.size() > 1){
	// setting support
	ti[indvs->first]->support = true;
	// looping over reads
	for(vector<string>::iterator readName = ti[indvs->first]->cluster[ci->first].begin(); readName != ti[indvs->first]->cluster[ci->first].end(); readName++ ){
	  ti[indvs->first]->badFlag[(*readName)] = 1;
	}
      }
    }
  }

  return true;
}

  bool score(string seqid, long int pos, readPileUp & targDat, long int cp,  insertDat & localDists, string & results, global_opts localOpts){

    map < string, indvDat*> ti, bi;

    indvDat * i;
    i = new indvDat;
    initIndv(i);
    ti[globalOpts.targetBam] = i;

    stringstream oddInsert;

    loadIndv(ti, targDat, localOpts, localDists, pos);
    if(ti[globalOpts.targetBam]->nReads < 3){
      return true;
    }

    for( map< int, vector<string> >::iterator ci = ti[globalOpts.targetBam]->cluster.begin(); ci != ti[globalOpts.targetBam]->cluster.end(); ci++){
      oddInsert << "cluster:";
      oddInsert << ci->first << "->" << ci->second.size() << ",";
    }

//  CJenksBreaks TwoAlleles   (&ti[globalOpts.targetBam]->hInserts, 2);
//  CJenksBreaks ThreeAlleles (&ti[globalOpts.targetBam]->hInserts, 3);
//    
//  vector <long> * TwoBreaks     = TwoAlleles.get_Results();
//  vector <long> * ThreeBreaks   = ThreeAlleles.get_Results();


  vector< vector < double > > targetGls, backgroundGls, totalGls;
  vector<string> genotypes;

  vector<int> depth;
  vector<int> sameStrand;
  vector<int> inserts;
  vector<int> nMissingMate;

  int nAlt = 0;

  string genotype;
  int ss     = 0;
  int insert = 0;

  double ta        = 0;
  double tb        = 0;
  double aa        = 0;
  double ab        = 0;
  double tgc       = 0;
  double targetAlt = 0;
  double targetRef = 0;

  vector<double> Targ;

  nMissingMate.push_back(ti[globalOpts.targetBam]->mateMissing);
  depth.push_back(ti[globalOpts.targetBam]->nReads);
  nAlt += processGenotype(ti[globalOpts.targetBam], Targ, &genotype, &targetRef, &targetAlt, &ta, &tb, &aa, &ab, &tgc, localDists, &ss, &insert);
  genotypes.push_back(genotype);
      

  double clipper =   double(ti[globalOpts.targetBam]->clipped) / double(ti[globalOpts.targetBam]->lengthSum);
  
  oddInsert << clipper << ";";
  

  stringstream tmpOutput;

  tmpOutput  << seqid   << "\t" ;       // CHROM
  tmpOutput  << pos +1  << "\t" ;       // POS
  tmpOutput  << "."     << "\t" ;       // ID
  tmpOutput  << "NA"    << "\t" ;       // REF
  tmpOutput  << "SV"    << "\t" ;       // ALT
  tmpOutput  << "."     << "\t" ;       // QUAL
  tmpOutput  << "."     << "\t" ;       // FILTER
  tmpOutput  << nAlt     << "\t" ;       // INFO
  tmpOutput  << "GT:GL:MM:DP" << "\t" ; // FORMAT
  tmpOutput  << genotype << "\t";
  tmpOutput  << ti[globalOpts.targetBam]->notMapped   << "\t"
	     << ti[globalOpts.targetBam]->nAboveAvg   << "\t"
	     << ti[globalOpts.targetBam]->sameStrand  << "\t"
	     << ti[globalOpts.targetBam]->mateMissing << "\t"
    //         << clipper << "\t" 
    //             << (*TwoBreaks)[1] << ":" << (*ThreeBreaks)[2] <<  ":" << ti[globalOpts.targetBam]->hInserts.size() << "\t"
         << oddInsert.str() << "\t"
    // 	     << ti[globalOpts.targetBam]->maxLength  << "\t"
	     << ti[globalOpts.targetBam]->insertSum  / double (ti[globalOpts.targetBam]->mappedPairs)  << "\t"
	     << ti[globalOpts.targetBam]->nReads     << endl;
  
  

  results.append(tmpOutput.str());
      
  delete ti[globalOpts.targetBam];
  return true;
}

 
bool runRegion(int seqidIndex, int start, int end, vector< RefData > seqNames){
  
  string regionResults;


  omp_set_lock(&lock);

  global_opts localOpts = globalOpts;
  insertDat localDists = insertDists;

  omp_unset_lock(&lock);

  BamReader targetReader;
  
  prepBams(targetReader, "target");

  
  int setRegionErrorFlag = 0;

  setRegionErrorFlag += targetReader.SetRegion(seqidIndex, start, seqidIndex, end);

  if(setRegionErrorFlag != 1){
    return false;
  }
    
  BamAlignment alt;
  
  readPileUp targetPileUp;

  if(! targetReader.GetNextAlignment(alt)){
    targetReader.Close();
    return false;
  }

  long int currentPos = -1;

  targetPileUp.processAlignment(alt,     currentPos);


  bool getTarget     = true;

  while(1){
    
    while(currentPos > targetPileUp.CurrentStart && getTarget){
      getTarget = targetReader.GetNextAlignment(alt);
      targetPileUp.processAlignment(alt, currentPos);
    }
    
    if(getTarget == false ){
      break;
    }
    
    targetPileUp.purgePast();

    if(! score(seqNames[seqidIndex].RefName, currentPos, targetPileUp, currentPos, localDists, regionResults, localOpts )){
      cerr << "FATAL: problem during scoring" << endl;
      cerr << "FATAL: wham exiting"           << endl;
      exit(1);
    }

    currentPos += 1;

    if(regionResults.size() > 100000){
      omp_set_lock(&lock);
      cout << regionResults;
      omp_unset_lock(&lock);
      regionResults.clear();
      
    }

  }

  omp_set_lock(&lock);

  cout << regionResults;

  omp_unset_lock(&lock);
  
  regionResults.clear();
  
  targetReader.Close();
  
  return true;
}


int main(int argc, char** argv) {

  omp_init_lock(&lock);

  srand((unsigned)time(NULL));

  globalOpts.nthreads = -1;

  parseOpts(argc, argv);
  
  if(globalOpts.nthreads == -1){
  
  }
  else{
    omp_set_num_threads(globalOpts.nthreads);
  }
 

  BamReader targetReader;
  
  prepBams(targetReader, "target");

  if(!getInsertDists()){
    cerr << "FATAL: " << "problem while generating insert lengths dists" << endl;
    exit(1);
  }
  
  RefVector sequences = targetReader.GetReferenceData();

  targetReader.Close();

  printHeader();

  int seqidIndex = 0;

  if(globalOpts.region.size() == 2){
    for(vector< RefData >::iterator sit = sequences.begin(); sit != sequences.end(); sit++){      
      if((*sit).RefName == globalOpts.seqid){
	break;
      }
      seqidIndex += 1;
    }
  }

  if(seqidIndex != 0){
    if(! runRegion(seqidIndex, globalOpts.region[0], globalOpts.region[1], sequences)){
      cerr << "WARNING: region failed to run properly." << endl;
    }
    cerr << "INFO: WHAM-BAM finished normally." << endl;
    return 0;
  }
  
  vector< regionDat* > regions; 

  for(vector< RefData >::iterator sit = sequences.begin(); sit != sequences.end(); sit++){
    int start = 0;
    if((*sit).RefLength < 1000){
      continue;
    }

    for(;start < ((*sit).RefLength) ; start += 10000000){
      regionDat * chunk = new regionDat;
      chunk->seqidIndex = seqidIndex;
      chunk->start      = start;
      chunk->end        = start + 10000000 ;
      regions.push_back(chunk);
    }
    regionDat * lastChunk = new regionDat;
    lastChunk->seqidIndex = seqidIndex;
    lastChunk->start = start;
    lastChunk->end   = (*sit).RefLength;
    seqidIndex += 1;
    if(start < (*sit).RefLength){
      regions.push_back(lastChunk);
    }
  }

 #pragma omp parallel for
  
  for(int re = 0; re < regions.size(); re++){

    //    cerr << regions[re]->seqidIndex << "\t" << regions[re]->start << "\t" << regions[re]->end << endl;

    if(! runRegion( regions[re]->seqidIndex, regions[re]->start, regions[re]->end, sequences)){
      omp_set_lock(&lock);
      cerr << "WARNING: region failed to run properly: " 
	   << sequences[regions[re]->seqidIndex].RefName 
	   << ":"  << regions[re]->start << "-" 
	   << regions[re]->end 
	   <<  endl;
      omp_unset_lock(&lock);
    }

  }

    //    (*chunk) = NULL;
    //    delete (*chunk);


  cerr << "INFO: WHAM-BAM finished normally." << endl;
  return 0;

}
