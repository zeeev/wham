#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <cmath>
#include <time.h>
#include <algorithm>
#include "split.h"
#include "KMERUTILS.h"

// openMP - swing that hammer
#include <omp.h>

// bamtools and my headers
#include "api/BamMultiReader.h"
#include "readPileUp.h"

// msa headers
#include <seqan/align.h>
#include <seqan/graph_msa.h>

using namespace std;
using namespace BamTools;

struct cigar{
  char type;
  unsigned int length;
};

struct regionDat{
  int seqidIndex ;
  int start      ;
  int end        ;
};

struct indvDat{
  bool   support       ;
  string genotype      ;
  int    genotypeIndex ; 
  int    nReads        ;
  int    mappedPairs   ;
  int    nAboveAvg     ;
  int    notMapped     ;
  int    mateMissing   ;
  int    sameStrand    ;
  int    nClipping     ; 
  int    mateCrossChromosome ;
  int    maxLength      ;
  double    nBad        ;
  double    nGood       ;
  double insertSum ;
  double insertMean ;
  double lengthSum ;
  double clipped ;
  vector<double> inserts;
  vector<double> hInserts;
  vector<long double> gls;
  vector< BamAlignment > alignments;
  vector<int> badFlag;
  vector<int> MapQ;
  map<int, vector<string> > cluster;
};

struct insertDat{
  map<string, double> mus; // mean of insert length for each indvdual across 1e6 reads
  map<string, double> sds;  // standard deviation
  map<string, double> lq ;  // 25% of data
  map<string, double> up ;  // 75% of the data
  map<string, double> avgD; 
  double overallDepth;
} insertDists;

struct global_opts {
  vector<string> targetBams    ;
  vector<string> backgroundBams;
  vector<string> all           ;
  int            nthreads      ;
  string         seqid         ;
  string         bed           ; 
  string         mask          ;
  vector<int>    region        ; 
} globalOpts;


struct info_field{
  double lrt;

  double taf;
  double baf;
  double aaf;

  double nat;
  double nbt;

  double nab;
  double nbb;

  double tgc;
  double bgc;

};

template<typename T>
inline bool aminan(T value)
{
  return value != value;

}

static const char *optString ="ht:b:r:x:e:m:";

// this lock prevents threads from printing on top of each other

omp_lock_t lock;

bool sortStringSize(string i, string j) {return (i.size() < j.size());}

//ripped from EKG @ github Shannon Entropy

double entropy(string& st) {
  vector<char> stvec(st.begin(), st.end());
  set<char> alphabet(stvec.begin(), stvec.end());
  vector<double> freqs;
  for (set<char>::iterator c = alphabet.begin(); c != alphabet.end(); ++c) {
    int ctr = 0;
    for (vector<char>::iterator s = stvec.begin(); s != stvec.end(); ++s) {
      if (*s == *c) {
	++ctr;
      }
    }
    freqs.push_back((double)ctr / (double)stvec.size());
  }
  double ent = 0;
  double ln2 = log(2);
  for (vector<double>::iterator f = freqs.begin(); f != freqs.end(); ++f) {
    ent += *f * log(*f)/ln2;
  }
  ent = -ent;
  return ent;
}

void initInfo(info_field * s){
  s->lrt = 0;
  s->taf = 0.000001;
  s->baf = 0.000001;
  s->aaf = 0.000001;
  s->nat = 0;
  s->nbt = 0;
  s->nab = 0;
  s->nbb = 0;
  s->tgc = 0;
  s->bgc = 0;
}

void initIndv(indvDat * s){
  s->support             = false;
  s->genotype            = "./.";
  s->genotypeIndex       = -1;
  s->nBad                = 0;
  s->nClipping           = 0;
  s->nGood               = 0;
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

string collapseCigar(vector<CigarOp> & v){
  
  stringstream ss;

  for(vector<CigarOp>::iterator it = v.begin(); it != v.end(); it++){
    ss << (*it).Length;
    ss << (*it).Type;
    
  }
  
  return ss.str();

}

void printHeader(void){
  cout << "##fileformat=VCFv4.1"                                                                << endl;
  cout << "##INFO=<ID=LRT,Number=1,Type=Float,Description=\"Likelihood Ratio Test Statistic\">" << endl;
  cout << "##INFO=<ID=WAF,Number=3,Type=Float,Description=\"Allele frequency of: background,target,combined\">" << endl;
  cout << "##INFO=<ID=GC,Number=2,Type=Integer,Description=\"Number of called genotypes in: background,target\">"  << endl;
  cout << "##INFO=<ID=AT,Number=15,Type=Float,Description=\"Attributes for classification\">"                      << endl;
  cout << "##INFO=<ID=KM,Number=3,Type=Float,Description=\"Kmer filters. The number of 17bp kmers that do not contain N, the number of hits to the masking DB, the fraction of hits\">" << endl;
  cout << "##INFO=<ID=PU,Number=1,Type=Integer,Description=\"Number of reads read supporting position\">" << endl;
  cout << "##INFO=<ID=SU,Number=1,Type=Integer,Description=\"Number of supplement read supporting position\">" << endl;
  cout << "##INFO=<ID=CU,Number=1,Type=Integer,Description=\"Number of neighboring all soft clip clusters across all individuals at pileup position \">" << endl;
  cout << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Number of reads at pileup position across individuals passing filters\">" << endl;
  cout << "##INFO=<ID=SP,Number=1,Type=String,Description=\"Support for endpoint;  none:., mp:mate pair, sr:split read, al:alternative alignment\">" << endl;
  cout << "##INFO=<ID=BE,Number=3,Type=String,Description=\"Best end position: chr,position,count or none:.\">"                << endl;
  cout << "##INFO=<ID=DI,Number=1,Type=Character,Description=\"Consensus is from front or back of pileup : f,b\">"          << endl;
  cout << "##INFO=<ID=NC,Number=1,Type=String,Description=\"Number of soft clipped sequences collapsed into consensus\">"   << endl;
 cout << "##INFO=<ID=MQ,Number=1,Type=String,Description=\"Average mapping quality\">"   << endl;
 cout << "##INFO=<ID=MQF,Number=1,Type=String,Description=\"Fraction of reads with MQ less than 50\">"   << endl;
  cout << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">"      << endl;
  cout << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT allele\">"         << endl;
  cout << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"                                                  << endl;
  cout << "##FORMAT=<ID=GL,Number=A,Type=Float,Description=\"Genotype likelihood\">"                                        << endl;
  cout << "##FORMAT=<ID=NR,Number=1,Type=Integer,Description=\"Number of reads that do not support a SV\">"                 << endl;
  cout << "##FORMAT=<ID=NA,Number=1,Type=Integer,Description=\"Number of reads supporting a SV\">"                          << endl;
  cout << "##FORMAT=<ID=NS,Number=1,Type=Integer,Description=\"Number of reads with a softclip at POS for individual\">"    << endl;
  cout << "##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Number of reads passing filters\">"                          << endl;
  cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" << "\t";

  for(unsigned int b = 0; b < globalOpts.all.size(); b++){
    cout << globalOpts.all[b] ;
    if(b < globalOpts.all.size() - 1){
      cout << "\t";
    }
  }
  cout << endl;  
}

string joinComma(vector<string> & strings){

  stringstream ss;

  for(unsigned int i = 0; i < strings.size() -1 ; i++){
    ss << strings[i] << "," ;
  }

  ss << strings.back();

  return ss.str();

}

string printIndvDat(  indvDat * d ){

  stringstream ss;

  ss << "Genotype       : .......... " << d->genotype          << endl;
  ss << "Genotype index : .......... " << d->genotypeIndex     << endl;
  ss << "Number of reads: .......... " << d->badFlag.size()    << endl;
  ss << "Number of mapped mates: ... " << d->mappedPairs       << endl;
  ss << "Number of odd insert size:  " << d->nAboveAvg         << endl;
  ss << "Number of reads not mapped: " << d->notMapped         << endl;
  ss << "Number of mates missing: .. " << d->mateMissing       << endl;
  ss << "Number of mates same strand:" << d->sameStrand        << endl;
  ss << "Number of reads clipped: .. " << d->nClipping         << endl;
  ss << "Number of alternative reads:" << d->nBad              << endl;
  ss << "Number of reference reads:  " << d->nGood             << endl;
  ss << "Genotype likelihoods: ..... " << d->gls[0] 
     << ":" << d->gls[1] 
     << ":" << d->gls[2] 
     << endl;
  ss << endl;
  ss << "Reads: " << endl;

  int z = 0;

  for(vector< BamAlignment >::iterator it = d->alignments.begin(); it != d->alignments.end(); it++){
   
    ss   << " " << (*it).Name << " " 
         << d->badFlag[z] << " "
         << (*it).RefID    << " " 
	 << (*it).Position << " " 
	 << (*it).GetEndPosition(false,true) << " "
	 << (*it).Position + (*it).Length << " "
         << (*it).MapQuality << " " 
         << (*it).MateRefID    << " " 
         << (*it).MatePosition << " " 
	 << collapseCigar((*it).CigarData) << " " 
         << (*it).AlignmentFlag << " " 
	 << (*it).QueryBases 
	 << endl;
    z++;
  }

  return ss.str();

}

string join(vector<string> strings){

  string joined = "";

  for(vector<string>::iterator sit = strings.begin(); sit != strings.end(); sit++){
    joined = joined + " " + (*sit) + "\n";
  }
  return joined;
}

void printVersion(void){
  cerr << "Version: " << VERSION << endl;
  cerr << "Contact: zev.kronenberg [at] gmail.com " << endl;
  cerr << "Notes  : -If you find a bug, please open a report on github!" << endl;
  cerr << endl;
}

void printHelp(void){
  cerr << "usage  : WHAM-BAM -m <STRING> -x <INT> -r <STRING>     -e <STRING>  -t <STRING>    -b <STRING>   " << endl << endl;
  cerr << "example: WHAM-BAM -m microSat_and_simpleRep_hg19.wham.masking.txt -x 20 -r chr1:0-10000 -e genes.bed -t a.bam,b.bam -b c.bam,d.bam" << endl << endl; 

  cerr << "required   : t <STRING> -- comma separated list of target bam files"           << endl ;
  cerr << "recommended: m <STRING> -- kmer database for downstream filtering"             << endl ; 
  cerr << "option     : b <STRING> -- comma separated list of background bam files"       << endl ;
  cerr << "option     : r <STRING> -- a genomic region in the format \"seqid:start-end\"" << endl ;
  cerr << "option     : x <INT>    -- set the number of threads, otherwise max          " << endl ; 
  cerr << "option     : e <STRING> -- a bedfile that defines regions to score           " << endl ; 
  cerr << endl;
  printVersion();
}

void parseOpts(int argc, char** argv){
  int opt = 0;

  globalOpts.mask = "NA";
  globalOpts.bed  = "NA";

  opt = getopt(argc, argv, optString);

  while(opt != -1){
    switch(opt){
    case 'm':
      {
	globalOpts.mask = optarg;
	cerr << "INFO: WHAM-BAM will screen breakpoints for simple repeats and microstats: " << globalOpts.mask << endl;
	break;
      }
    case 'e':
      {
	globalOpts.bed = optarg;
	cerr << "INFO: WHAM-BAM will only score within bed coordiates provided: " << globalOpts.bed << endl;
	break;
      }
    case 'x':
      {
	globalOpts.nthreads = atoi(((string)optarg).c_str());
	cerr << "INFO: OpenMP will roughly use " << globalOpts.nthreads << " threads" << endl;
	break;
      }
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
  if( globalOpts.targetBams.empty() && globalOpts.backgroundBams.empty() ){
    cerr << "FATAL: Failure to specify target and/or background bam files." << endl;
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

// gerates per bamfile statistics 

void grabInsertLengths(string targetfile){

  vector<string> target_group ;
  target_group.push_back(targetfile);
  
  vector<double> alIns;
  vector<double> nReads;

  BamMultiReader bamR;
  if(!bamR.Open(target_group) && bamR.LocateIndexes() ){
    cerr << "FATAL: cannot read - or find index for: " << targetfile << endl;
    exit(0);
  }

  SamHeader SH = bamR.GetHeader();
  if(!SH.HasSortOrder()){
    cerr << "FATAL: sorted bams must have the @HD SO: tag in each SAM header." << endl;
    exit(1);
  }

  RefVector sequences = bamR.GetReferenceData();

  bamR.Close();

  int i = 0; // index for while loop
  int n = 0; // number of reads
  
  BamAlignment al;

  while(i < 3 || n < 20000){

    const int randomChr = rand() % (sequences.size() -1);
    const int randomPos = rand() % (sequences[randomChr].RefLength -1);
    const int randomEnd = randomPos + 10000;

    if(randomEnd > sequences[randomChr].RefLength){
      continue; 
    }
    if(sequences[randomChr].RefLength < 10000){
      cerr << "FATAL: Trying to randomly sample depth from the first seqid: " << sequences[randomChr].RefName << endl;
      cerr << "      " << sequences[randomChr].RefName << " is " << sequences[randomChr].RefLength << " bp "; 
      cerr << "       Current wham needs the seqid to be longer than 10kb, please contact zev kronenberg if you get this error" << endl << endl;                 
    }
    
    BamMultiReader bamR;
    if(!bamR.Open(target_group) && bamR.LocateIndexes() ){
      cerr << "FATAL: cannot read - or find index for: " << targetfile << endl;
      exit(0);
    }
    if(! bamR.SetRegion(0, randomPos, 0, randomEnd)){      
      continue;
    }
        
    if(!bamR.GetNextAlignmentCore(al)){
      continue;
    }

    i++;
    
    long int cp = al.GetEndPosition();
    
    readPileUp allPileUp;
    
    while(bamR.GetNextAlignmentCore(al)){
      if(al.Position > cp){
	allPileUp.purgePast(&cp);
	cp = al.Position;
	nReads.push_back(allPileUp.currentData.size());
      }
      if(al.IsMapped()
	 && al.IsMateMapped() 
	 && abs(double(al.InsertSize)) < 10000 
	 && al.RefID == al.MateRefID
	 ){	
	allPileUp.processAlignment(al);
	alIns.push_back(abs(double(al.InsertSize)));
	n++;
      }
    }
    bamR.Close();
  }

  double mu       = mean(alIns        );
  double mud      = mean(nReads       );
  double variance = var(alIns, mu     );
  double sd       = sqrt(variance     );
  double sdd      = var(nReads, mud   );
  
  omp_set_lock(&lock);

  insertDists.mus[target_group[0]]  = mu;
  insertDists.sds[target_group[0]]  = sd;
  insertDists.avgD[target_group[0]] = mud;

  cerr << "INFO: for file:" << target_group[0] << endl            
       << "     " << target_group[0] << ": mean depth: " << mud << endl
       << "     " << target_group[0] << ": sd   depth: " << sdd << endl
       << "     " << target_group[0] << ": mean insert length: " << insertDists.mus[target_group[0]] << endl
       << "     " << target_group[0] << ": sd   insert length: " << insertDists.sds[target_group[0]] << endl
       << "     " << target_group[0] << ": number of reads used: " << n  << endl << endl;
       
  omp_unset_lock(&lock);
}

void prepBams(BamMultiReader & bamMreader, string group){

  string errorMessage ;

  vector<string> files;

  if(group == "target"){
    files = globalOpts.targetBams;
  }
  if(group == "background"){
    files = globalOpts.backgroundBams;
  }
  if(group == "all"){
        files = globalOpts.all;
  }
  if(files.empty()){
    cerr << "FATAL: no files ?" << endl;
    exit(1);
  }

  bool attempt = true;
  int  tried   = 0   ;
  
  while(attempt && tried < 500){
    
    //    sleep(int(rand() % 10)+1);

    if( bamMreader.Open(files) && bamMreader.LocateIndexes() ){
      attempt = false;
    }
    else{
      tried += 1;
    }
  }
  if(attempt == true){
    cerr << "FATAL: unable to open BAMs or indices after: " << tried << " attempts"  << endl;  
    cerr << "Bamtools error message:\n" <<  bamMreader.GetErrorString()  << endl;
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


double unphred(double p){
  return pow(10, (-p/10));
}

bool processGenotype(string & fname,
		     indvDat * idat, 
		     double * totalAlt, 
		     double * totalAltGeno,
		     double * relativeDepth,
		     insertDat * stats){

  string genotype = "./.";

  long double aal = 0;
  long double abl = 0;
  long double bbl = 0;

  if(idat->badFlag.size() < 3){
    idat->gls.push_back(-255.0);
    idat->gls.push_back(-255.0);
    idat->gls.push_back(-255.0);
    return true;
  }

  double nref = 0.0;
  double nalt = 0.0;

  if(idat->badFlag.size() > 1000){
    std::random_shuffle ( idat->badFlag.begin(), idat->badFlag.end() );
  }

  int ri = 0;
 
  for(vector< int >::iterator rit = idat->badFlag.begin(); rit != idat->badFlag.end(); rit++){
    if(ri > 1000){
      break;
    }

    double mappingP = unphred(idat->MapQ[ri]);

    // will this improve stuff 1-> 2?

    if( (*rit == 1 && idat->nClipping > 1) ){
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
    ri++;
  }

  idat->nBad  = nalt;
  idat->nGood = nref;

  double nreads = idat->nReads;

  if(nreads > 1000){
    nreads = 1000;
  }

  aal = aal - log(pow(2,nreads)); // the normalization of the genotype likelihood
  abl = abl - log(pow(2,nreads)); // this is causing underflow for really high depth.
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
  genotype     = "0/0";
  idat->genotypeIndex = 0;

  if(abl > max){
    genotype = "0/1";
    idat->genotypeIndex = 1;
    (*totalAlt) += 1;
    max = abl;
  }
  if(bbl > max){
    genotype = "1/1";
    idat->genotypeIndex = 2;
    (*totalAlt) += 2;
  }
  if(idat->genotypeIndex != 0){
    *relativeDepth += (idat->badFlag.size() / stats->avgD[fname]);
    *totalAltGeno += 1;
  }

  idat->genotype = genotype;
  idat->gls.push_back(aal);
  idat->gls.push_back(abl);
  idat->gls.push_back(bbl);

#ifdef DEBUG
   cerr << genotype << "\t" << aal << "\t" << abl << "\t" << bbl << "\t" << nref << "\t"  << endl;
#endif
  return true;

 
}

bool checkN(string & s){

  int sl = s.size();
  int n  = 0;
  int i  = 0;

  if(sl < 20){
    return true;
  }
  for(i = 0 ; i < 10; i++){
    if(s[i] == 'N'){
      n += 1;
    }
  }
  for(i = sl - 1; i < sl - 10; i--){
    if(s[i] == 'N'){
      n += 1;
    }
  }

  if(n > 5){
    return true;
  }

  return false;
}

bool loadIndv(map<string, indvDat*> & ti, 
	      readPileUp & pileup, 
	      global_opts & localOpts, 
	      insertDat & localDists, 
	      long int * pos
	      
	      ){    

  
  for(list<BamAlignment>::iterator r = pileup.currentData.begin(); r != pileup.currentData.end(); r++){
   
    if((*r).Position > *pos){
      continue;
    }

    if( ((*r).AlignmentFlag & 0x0800) != 0 || ! (*r).IsPrimaryAlignment()){
      continue;
    }

    string fname = (*r).Filename;

    int bad = 0;

    vector< CigarOp > cd = (*r).CigarData;

    if( ((pileup.primary[(*r).Position].size() > 1) || (pileup.primary[(*r).GetEndPosition(false,true)].size() > 1))
	&& (cd.front().Type == 'S' || cd.back().Type == 'S') ){
      bad = 1;
      ti[fname]->nClipping++;
    }
    
    if((*r).IsMapped() && (*r).IsMateMapped() && ((*r).RefID == (*r).MateRefID)){

      ti[fname]->insertSum    += abs(double((*r).InsertSize));
      ti[fname]->mappedPairs  += 1;
      
      if(( (*r).IsReverseStrand() && (*r).IsMateReverseStrand() ) || ( !(*r).IsReverseStrand() && !(*r).IsMateReverseStrand() )){
	bad = 1;
	ti[fname]->sameStrand += 1;
      }
      
      if((*r).IsReverseStrand() && ! (*r).IsMateReverseStrand() && (*r).Position < (*r).MatePosition){
	pileup.evert++;
      }

      double ilength = abs ( double ( (*r).InsertSize ));
      
      double iDiff = abs ( ilength - localDists.mus[(*r).Filename] );
      
      ti[fname]->inserts.push_back(ilength);
      
      if(iDiff > (3.0 * insertDists.sds[(*r).Filename]) ){
	bad = 1;
	ti[fname]->nAboveAvg += 1;
	ti[fname]->hInserts.push_back(ilength);
      
	if( ilength - localDists.mus[(*r).Filename] < 0){
	  pileup.mateTooClose++;
	}
	if( ilength - localDists.mus[(*r).Filename] > 0){
	  pileup.mateTooFar++;
	}
      }
    }

    ti[fname]->nReads++;

    map<string, int>::iterator os = pileup.odd.find((*r).Name);

    if(os !=  pileup.odd.end()){
      bad = 1;
    }

    if(bad == 1){
      ti[fname]->nBad += 1;
    }
    else{
      ti[fname]->nGood += 1;
    }
    ti[fname]->badFlag.push_back( bad );
#ifdef DEBUG
    ti[fname]->alignments.push_back(*r);
#endif
    ti[fname]->MapQ.push_back((*r).MapQuality);
  }
  return true;
}


bool cleanUp( map < string, indvDat*> & ti, global_opts localOpts){
  for(vector<string>::iterator all = localOpts.all.begin(); all != localOpts.all.end(); all++ ){
    delete ti[*all];
  }
  return true;
}

bool loadInfoField(map<string, indvDat*> dat, info_field * info, global_opts & opts){
  
  for(unsigned int b = 0; b < opts.backgroundBams.size(); b++){

    if( dat[opts.backgroundBams[b]]->genotypeIndex == -1){
      continue;
    }
    info->nat += 2 - dat[opts.backgroundBams[b]]->genotypeIndex;
    info->nbt +=     dat[opts.backgroundBams[b]]->genotypeIndex;
    info->tgc += 1;
  }

  for(unsigned int t = 0; t < opts.targetBams.size(); t++){
    if(dat[opts.targetBams[t]]->genotypeIndex == -1){
      continue;
    }
    info->nab += 2 - dat[opts.targetBams[t]]->genotypeIndex;
    info->nbb +=     dat[opts.targetBams[t]]->genotypeIndex;
    info->bgc += 1;
  }

  info->taf += info->nbt / (info->nat + info->nbt);
  info->baf += info->nbb / (info->nab + info->nbb);
  info->aaf += (info->nbt + info->nbb) / (info->nat + info->nbt + info->nab + info->nbb);

  double alt  = logLbinomial(info->nbt, (info->tgc * 2), info->taf) + logLbinomial(info->nbb, (2* info->bgc), info->baf);
  double null = logLbinomial(info->nbt, (info->tgc * 2), info->aaf) + logLbinomial(info->nbb, (2* info->bgc), info->aaf);

  info->lrt = 2 * (alt - null);
  
  if(info->lrt < 0 ){
    info->lrt = 0;
  }
  if(aminan(info->lrt)){
    info->lrt = 0;
  }
  return true;
}


string infoText(info_field * info){
  
  stringstream ss;

  ss << "LRT=" << info->lrt << ";";
  if(aminan(info->baf)){
    ss << "WAF=" << info->taf << ",.," << info->aaf << ";";
  }
  else if(aminan(info->taf)){  
    ss << "WAF=" << ".," << info->baf << "," << info->aaf << ";";
  }
  else{
    ss << "WAF=" << info->taf << "," << info->baf << "," << info->aaf << ";";
  }
 
  ss << "GC="  << info->tgc << "," << info->bgc << ";";

  return ss.str();

}

bool burnCigar(string s, vector<cigar> & cigs){
  
  unsigned int offset = 0;
  unsigned int index  = 0; 

  for(string::iterator it = s.begin() ; it != s.end(); it++){

    switch(*it){
    case 'M':
      {
	cigar c;
	c.type   = 'M';
	c.length = atoi(s.substr(offset, index-offset).c_str());
	cigs.push_back(c);
	offset = index + 1;
	break;
      }
    case 'I':
      {
        cigar c;
	c.type   = 'I';
        c.length = atoi(s.substr(offset, index-offset).c_str());
	cigs.push_back(c);
	offset = index + 1;
	break;
      }
    case 'D':
      {
        cigar c;
	c.type   = 'D';
        c.length = atoi(s.substr(offset, index-offset).c_str());
	cigs.push_back(c);
	offset = index + 1;
	break;
      }
    case 'N':
      {
        cigar c;
	c.type   = 'N';
        c.length = atoi(s.substr(offset, index-offset).c_str());
	cigs.push_back(c);
	offset = index + 1;
	break;
      }
    case 'S':
      {
        cigar c;
	c.type   = 'S';
        c.length = atoi(s.substr(offset, index-offset).c_str());
	cigs.push_back(c);
	offset = index + 1;
	break;
      }
    case 'H':
      {
        cigar c;
	c.type   = 'H';
        c.length = atoi(s.substr(offset, index-offset).c_str());
	cigs.push_back(c);
	offset = index + 1;
	break;
      }
    case 'P':
      {
        cigar c;
	c.type   = 'P';
        c.length = atoi(s.substr(offset, index-offset).c_str());
	cigs.push_back(c);
	offset = index + 1;
	break;
      }
    case 'X':
      {
        cigar c;
	c.type   = 'X';
        c.length = atoi(s.substr(offset, index-offset).c_str());
	cigs.push_back(c);
	offset = index + 1;
	break;
      } 
    case '=':
      {
        cigar c;
	c.type   = '=';
        c.length = atoi(s.substr(offset, index-offset).c_str());
	cigs.push_back(c);
	offset = index + 1;
	break;
      }
    default:
      break;
    }
    index += 1;
  }
  //cerr << "done burning" << endl;
  return true;
}

void endPos(vector<cigar> & cigs, long int * pos){
  
  for(vector<cigar>::iterator it = cigs.begin(); 
      it != cigs.end(); it++){
    
    switch( (*it).type ){
    case 'I':
      {
	*pos += (*it).length;
	break;
      }
    case 'M':
      {
	*pos += (*it).length;
	break;
      }
    default:
      break;
    }
  }
}



int SplitReadEndFinder(long int * pos,
	       map<long int, vector < BamAlignment > > & supliment, 
	       string & bestEnd,
	       string & bestSeqid,
	       int * support,
	       long int * otherPos,
	       string & currentSeqid
	       ){
  
#ifdef DEBUG
  cerr << "N secondary:" << supliment.size() << endl;
#endif

  if(supliment[*pos].empty()){
    return 0;
  }

  // seqid, pos, number of reads supporting 
  map<string, map< long int, int > > otherPositions;

#ifdef DEBUG
  cerr << "Chimeric mapping:" << endl;
#endif

  for(vector<BamAlignment>::iterator it = supliment[*pos].begin(); 
      it != supliment[*pos].end(); it++){

    string saTag;
    if(! (*it).GetTag("SA", saTag)){
      cerr << "FATAL::INTERNAL: no sa\n";
      exit(1);
    }
    
#ifdef DEBUG
    cerr << saTag << endl;
#endif

    vector<string> saInfo = split(saTag, ";");

    for(vector<string>::iterator ch = saInfo.begin(); 
	ch != saInfo.end(); ch++){
      
      if((*ch).empty()){
	break;
      }
      
      vector<string> chimera  = split((*ch), ",");

      vector<cigar> c;

      burnCigar(chimera[3], c);

      if(c.front().type == 'S'){
	if(otherPositions[chimera[0]].find(atoi(chimera[1].c_str())) 
	   == otherPositions[chimera[0]].end()){
	  otherPositions[chimera[0]][atoi(chimera[1].c_str())] = 1;
	}
	else{
	  otherPositions[chimera[0]][atoi(chimera[1].c_str())]++;
	}
      } 
      if(c.back().type  == 'S'){
	long int endP = atoi(chimera[1].c_str());
	endPos(c, &endP);
	if(otherPositions[chimera[0]].find(endP) ==
	   otherPositions[chimera[0]].end()
	   ){
	  otherPositions[chimera[0]][endP] = 1;	
	}
	else{
	  otherPositions[chimera[0]][endP]++;	
	}
      } 
    }
  }
  
  string bestOpt; 
  int    nbe = 0;
  
  // checking within chromosome

  for(map<long int, int>::iterator pos = otherPositions[currentSeqid].begin();
      pos != otherPositions[currentSeqid].end(); pos++){
    if(pos->second > nbe ){
      nbe =  pos->second;
      stringstream best;
      best << currentSeqid << "," << pos->first << "," << pos->second;
      *otherPos = pos->first;
      bestSeqid = currentSeqid;
      *support = pos->second;
      bestOpt = best.str();
    }
  }

  // checking across chromosomes

  if(nbe < 1){
    for(map< string, map <long int, int> >::iterator ab = otherPositions.begin();
	ab != otherPositions.end(); ab++){
      
      for(map<long int, int>::iterator pos = ab->second.begin();
	  pos != ab->second.end(); pos++){
	
	if(pos->second > nbe && nbe < 2){
	  nbe =  pos->second;
	  stringstream best;
	  best << ab->first << "," << pos->first << "," << pos->second;
	  *otherPos = pos->first;
	  bestSeqid = ab->first;
	  *support = pos->second;
	  bestOpt = best.str();
	}
      }
    }
  }
  if(!bestOpt.empty()){
    bestEnd = bestOpt;
  }
  
  return otherPositions.size();

}

//XATAG
int otherBreakAlternative(long int * pos,
			  map<long int, vector < BamAlignment > > & supliment, 
			  string & bestEnd,
			  string & bestSeqid,
			  int * support,
			  long int * otherPos,
			  string & currentSeqid			  
	       ){
  
#ifdef DEBUG
  cerr << "N alternative:" << supliment[*pos].size() << endl;
#endif

 
  map<string, map< long int, int > > otherPositions;

#ifdef DEBUG
  cerr << "Alternative mapping:" << endl;
#endif

  for(vector<BamAlignment>::iterator it = supliment[*pos].begin(); 
      it != supliment[*pos].end(); it++){

#ifdef DEBUG
    cerr << (*it).Name << endl;
#endif
    
    string saTag;
    if(!(*it).GetTag("XA", saTag)){
      continue;
    }
    
    vector<string> xaInfo = split(saTag, ";");
    
    for(vector<string>::iterator ch = xaInfo.begin(); 
	ch != xaInfo.end(); ch++){
      
      if((*ch).empty()){
	break;
      }
      
#ifdef DEBUG
      cerr << saTag << endl;
#endif

      vector<string> chimera  = split((*ch), ",");

      // remove the strand from XA
      chimera[1].erase(0,1); 

      vector<cigar> c;

      burnCigar(chimera[2], c);

      if(c.front().type == 'S'){
	if(otherPositions[chimera[0]].find(atoi(chimera[1].c_str())) 
	   == otherPositions[chimera[0]].end()){
	  otherPositions[chimera[0]][atoi(chimera[1].c_str())] = 1;
	}
	else{
	  otherPositions[chimera[0]][atoi(chimera[1].c_str())]++;
	}
      } 
      if(c.back().type  == 'S'){
	long int endP = atoi(chimera[1].c_str());
	endPos(c, &endP);
	if(otherPositions[chimera[0]].find(endP) ==
	   otherPositions[chimera[0]].end()
	   ){
	  otherPositions[chimera[0]][endP] = 1;
	}
	else{
	  otherPositions[chimera[0]][endP]++;
	}
      } 
    }
  }
  
  string bestOpt; 
  int    nbe = 0;

  // checking within chromosome

  for(map<long int, int>::iterator pos = otherPositions[currentSeqid].begin();
      pos != otherPositions[currentSeqid].end(); pos++){
    if(pos->second > nbe ){
      nbe =  pos->second;
      stringstream best;
      best << currentSeqid << "," << pos->first << "," << pos->second;
      *otherPos = pos->first;
      bestSeqid = currentSeqid;
      *support = pos->second;
      bestOpt = best.str();
    }
  }

  // checking across chromosome


  if(nbe < 1){
    for(map< string, map <long int, int> >::iterator ab = otherPositions.begin();
	ab != otherPositions.end(); ab++){
      
      for(map<long int, int>::iterator pos = ab->second.begin();
	  pos != ab->second.end(); pos++){
	
	if(pos->second > nbe && nbe < 1){
	  nbe =  pos->second;
	  stringstream best;
	  best << ab->first << "," << pos->first << "," << pos->second;
	  *otherPos = pos->first;
	  bestSeqid = ab->first;
	  *support = pos->second;
	  bestOpt = best.str();
	}
      }
    }
  }
  
  if(!bestOpt.empty()){
    bestEnd = bestOpt;
  }

  return otherPositions.size();

}

bool uniqClips(long int * pos, 
	       map<long int, vector < BamAlignment > > & clusters, 
	       vector<string> & alts, string & direction){

  map<string, vector<string> >  clippedSeqs;

  int bcount = 0;
  int fcount = 0;

  for( vector < BamAlignment >::iterator it = clusters[(*pos)].begin(); 
       it != clusters[(*pos)].end(); it++){
    
    if(((*it).AlignmentFlag & 0x0800) != 0 || !(*it).IsPrimaryAlignment()){
      continue;
    }
    
    vector< CigarOp > cd = (*it).CigarData;
  
    if((*it).Position == (*pos)){
      string clip = (*it).QueryBases.substr(0, cd.front().Length);
      if(clip.size() < 10){
	continue;
      }
      clippedSeqs["f"].push_back(clip);
      fcount += 1;
    }
    if((*it).GetEndPosition(false,true) == (*pos)){
      string clip = (*it).QueryBases.substr( (*it).Length - cd.back().Length );
      if(clip.size() < 10){
	continue;
      }
      clippedSeqs["b"].push_back(clip);    
      bcount += 1;
    }
  }

  string key = "f";
  if(bcount > fcount){
    key = "b";
  }

  direction = key;

  for(vector<string>::iterator seqs = clippedSeqs[key].begin();
      seqs != clippedSeqs[key].end(); seqs++
      ){
    alts.push_back(*seqs);
  }

  return true;
}



string consensus(vector<string> & s, double * nn, string & direction){

  if(s.empty()){
    return ".";
  }
  
  if(s.size() == 1){
    return s[0];
  }
  
  stringstream con;
  
  {
    using namespace seqan;
   
    typedef String< Dna > TSequence;
    typedef Graph<Alignment<StringSet<TSequence, Dependent<> > > > TGraph;

    StringSet<TSequence> seq;
    
    int index = 0;

    if(direction.compare("f") == 0){
      for(uint clips = 0; clips < s.size(); clips++){
	appendValue(seq, s[clips]);
	if(clips > 19){
	  break;
	}
      }
    }
    else{
      for(int clips = s.size() -1; clips > -1; clips--){
	appendValue(seq, s[clips]);          
	index+=1;
	if(index > 19){
	  break;
	}
      }
    }
    TGraph aliG(seq);
    globalMsaAlignment(aliG, Blosum62(-1, -11));

    String<char> align;
    convertAlignment(aliG, align);
    
    unsigned int nseq   = length(seq);
    unsigned int colLen = length(align) / nseq;
    
    for(unsigned int z = 0 ; z < colLen; z++){
      map<char, int> columnBases;
      for(unsigned int s = 0; s < nseq; s++){
	if(align[z + (s*colLen)] != gapValue<char>()){
	  columnBases[align[z + (s*colLen)]]++;
	}
      }
      if(columnBases.size() == 1){
	con << columnBases.begin()->first;
      }
      else{
	*nn += 1;
	con << "N";
      }
    }
    
#ifdef DEBUG
  cerr << "seqAn alignment" << endl;
  cerr << aliG << endl;
#endif 
  }
  
  return con.str(); 
}

bool clusterMatePos(string & seqid, 
		    long int * pos, 
		    map<long int, vector < BamAlignment > > & primary,
		    string & bestEnd, 
		    int * count,
		    long int * breakpoint
		    ){
  
  int otherSeqids = 0;  
  
  map<long int, int> otherPos;
 
  map<long int, int>::iterator fm;

  for(vector<BamAlignment>::iterator it = primary[*pos].begin();
      it != primary[*pos].end(); it++){

    if(!(*it).IsMateMapped()){
      continue;
    }

    if((*it).MateRefID != (*it).MateRefID){
      otherSeqids++;
      continue;
    }
    
    if(otherPos.find((*it).MatePosition) == otherPos.end()){
      otherPos[(*it).MatePosition] = 1;
    }
    else{
      otherPos[(*it).MatePosition]++;
    }
  }
 
  int maxCount = 0;
  long int bestPos = 0;

  for(map<long int, int>::iterator pp = otherPos.begin(); 
      pp != otherPos.end(); pp++){    
    if(pp->second > maxCount){
      maxCount = pp->second;
      bestPos = pp->first;
    }
  }
  if(maxCount < 2){
    return false;
  }
  else{
    stringstream ss ;
    ss << seqid << "," << bestPos << "," << maxCount;
    *breakpoint = bestPos;
    *count = maxCount;
    bestEnd = ss.str();

    return true;
  }   
}

bool score(string seqid, 
	   long int * pos, 
	   readPileUp & totalDat, 
	   insertDat & localDists, 
	   string & results, 
	   global_opts localOpts,
	   vector<uint64_t> & kmerDB
	   ){

  
  totalDat.processPileup(pos);
  
  if(totalDat.primary[*pos].size() < 3){
    return true;
  }


  stringstream attributes;

  attributes << "AT="
             << double(totalDat.nPaired) / double(totalDat.numberOfReads)
             << ","
             << double(totalDat.nMatesMissing) / double(totalDat.numberOfReads)
             << ","
             << double(totalDat.nSameStrand) / double(totalDat.numberOfReads)
             << ","
             << double(totalDat.nCrossChr) / double(totalDat.numberOfReads)
             << ","
             << double(totalDat.nsplitRead) / double(totalDat.numberOfReads)
             << ","
             << double(totalDat.nf1SameStrand) / double(totalDat.numberOfReads)
             << ","
             << double(totalDat.nf2SameStrand) / double(totalDat.numberOfReads)
             << ","
             << double(totalDat.nf1f2SameStrand) / double(totalDat.numberOfReads)
             << ","
             << double(totalDat.nsplitReadCrossChr) / double(totalDat.numberOfReads)
             << ","
             << double(totalDat.nsplitMissingMates) / double(totalDat.numberOfReads)
             << ","
             << double(totalDat.nDiscordant) / double(totalDat.numberOfReads)
             << ","
             << double(totalDat.nsameStrandDiscordant) / double(totalDat.numberOfReads)
             << ","
             << double(totalDat.ndiscordantCrossChr) / double(totalDat.numberOfReads)
             << ","
             << double(totalDat.internalInsertion) / double(totalDat.numberOfReads)
             << ","
             << double(totalDat.internalDeletion) / double(totalDat.numberOfReads)
             << ","
             << double(totalDat.mateTooClose) / double(totalDat.numberOfReads)
             << ","
             << double(totalDat.mateTooFar) / double(totalDat.numberOfReads)
             << ","
             << double(totalDat.evert) / double(totalDat.numberOfReads);


  // if there are only good mate pairs skip this site (small indels and things that don't generate split reads
  //   are  going to be missed (NOT MANY))

  if(attributes.str().compare("AT=1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0") == 0){
    return true;
  }


  // finding consensus sequence 
  vector<string> alts ; // pairBreaks;

  string direction ;

  uniqClips(pos, totalDat.primary, alts, direction);

  if(alts.size() < 3){
    return true;
  }

  double nn   = 0;

  string altSeq = consensus(alts, &nn, direction);

  if(altSeq.size() < 10){
    return true;
  }

  if(nn / double(altSeq.size()) > 0.30 || nn > 18){
    return true;
  }

  // searchign for repeats 

  stringstream kfilter;

  double nReps  = 0;
  double nAssay = 0;

  if(altSeq.size() > 17){
    for(uint16_t l = 0; l < (altSeq.size() - 17); l++){
      string conKmer = altSeq.substr(l,17);
      std::size_t found = conKmer.find("N");
      if (found!=std::string::npos){
	continue;
      }
      nAssay += 1;
      char * con = new char[18];
      memcpy(con, conKmer.c_str(), 18);
      con[17] = '\0';
      uint64_t front =  charArrayToBin(con, 0);
      if( binary_search(kmerDB.begin(), kmerDB.end(), front) ){
	nReps+=1;
      }
      delete con;
    }
  }

  double kmHitFrac = double(nReps) / double(nAssay) ; 

  if(aminan(kmHitFrac)){
    kmHitFrac = 0;
  }
  kfilter << nAssay << "," << nReps << "," << kmHitFrac ;

  // trying to locate best end

  string bestSeqid ="" ;
  string bestEnd   ="" ;
  string esupport  ="" ;

  int otherBreakPointCount    = 0;
  long int otherBreakPointPos = 0;
  long int SVLEN              = -1;

  // trying to find mate breakpoint using splitread support
  int otherSeqids = SplitReadEndFinder(pos, totalDat.supplement, bestEnd, bestSeqid, &otherBreakPointCount, &otherBreakPointPos, seqid);

  if(!bestEnd.empty()){
    esupport = "sr";
  }

  // trying to find alternative mappings

  if(bestEnd.empty() || seqid.compare(bestSeqid) != 0 ){
    otherSeqids = otherBreakAlternative(pos, totalDat.primary, bestEnd, bestSeqid, &otherBreakPointCount, &otherBreakPointPos, seqid);
    if(!bestEnd.empty()){
      esupport = "al";
    }
  }

  // trying to find mate breakpoint using mate mapping postion.
  if(bestEnd.empty() || seqid.compare(bestSeqid) != 0){
    if(clusterMatePos(seqid, pos, totalDat.primary, bestEnd, &otherBreakPointCount, &otherBreakPointPos)){    
      bestSeqid = seqid;
      esupport = "mp";
    }
  }

  // bestSeqid is the other breakpoint
  // you need at least two reads for a translocation for a split read supported SV
  if(seqid.compare(bestSeqid) != 0 && ! bestSeqid.empty()){
    if(otherBreakPointCount < 2){
      return true;
    }
  }

  // SVs over a megabase require additional support 
  if(seqid.compare(bestSeqid) == 0 && ! bestSeqid.empty()){
    if(abs(*pos - otherBreakPointPos) > 1000000 && otherBreakPointCount < 2){
      return true;
    }
    else{
      SVLEN = abs((*pos+1) - otherBreakPointPos);
    }
  }

  // preparing data structure to load genotypes
  map < string, indvDat*> ti;

  for(unsigned int t = 0; t < localOpts.all.size(); t++){
    indvDat * i;
    i = new indvDat;
    initIndv(i);
    ti[localOpts.all[t]] = i;
  }

  loadIndv(ti, totalDat, localOpts, localDists, pos);

  double nAlt = 0;
  double nAltGeno = 0;

  double alternative_relative_depth_sum = 0;

  for(unsigned int t = 0; t < localOpts.all.size(); t++){
    processGenotype(localOpts.all[t], 
		    ti[localOpts.all[t]], 
		    &nAlt, 
		    &nAltGeno,
		    &alternative_relative_depth_sum , 
		    &localDists);
    #ifdef DEBUG
    cerr << "position: " << *pos << endl; 
    cerr << printIndvDat(ti[localOpts.all[t]]) << endl;
    #endif 
  }
  
  if(nAlt == 0 ){
    cleanUp(ti, localOpts);
    return true;
  }

  attributes << "," << (alternative_relative_depth_sum/nAltGeno) << ";";
  
  
  info_field * info = new info_field; 

  initInfo(info);

  loadInfoField(ti, info, localOpts);

  string infoToPrint = infoText(info);
  
  infoToPrint.append(attributes.str());

  stringstream tmpOutput;

  tmpOutput  << seqid           << "\t"  ;       // CHROM
  tmpOutput  << (*pos) +1       << "\t"  ;       // POS
  tmpOutput  << "."             << "\t"  ;       // ID
  tmpOutput  << "N"             << "\t"  ;       // REF
  tmpOutput  << altSeq          << "\t"  ;       // ALT
  tmpOutput  << "."             << "\t"  ;       // QUAL
  tmpOutput  << "."             << "\t"  ;       // FILTER
  tmpOutput  << infoToPrint                                                          ;
  tmpOutput  << "KM=" << kfilter.str()                    << ";"                     ;
  tmpOutput  << "PU=" << totalDat.primary[*pos].size()    << ";"                     ;
  tmpOutput  << "SU=" << totalDat.supplement[*pos].size() << ";"                     ;
  tmpOutput  << "CU=" << totalDat.primary.size() + totalDat.supplement.size() << ";" ; 
  tmpOutput  << "RD=" << totalDat.numberOfReads << ";"                               ;
  tmpOutput  << "NC=" << alts.size()     << ";"                                      ;
  tmpOutput  << "MQ=" << double(totalDat.mapQsum) / double(totalDat.numberOfReads)  << ";";
  tmpOutput  << "MQF=" << double(totalDat.nLowMapQ) / double(totalDat.numberOfReads)  << ";";

  if(esupport.empty() || bestEnd.empty() || (unsigned int)bestEnd.c_str()[0]  == 0){
    bestEnd  =  "." ;
    esupport =  "." ;
  }
  tmpOutput  << "SP=" << esupport  << ";";
  tmpOutput  << "BE=" << bestEnd   << ";";
  tmpOutput  << "DI="   << direction << ";";
  if(otherBreakPointPos == 0 || SVLEN == -1 ){
    tmpOutput << "END=.;SVLEN=.\t";
  }
  else{
    tmpOutput << "END=" << otherBreakPointPos << ";" << "SVLEN=" << SVLEN << "\t";
  }
  tmpOutput  << "GT:GL:NR:NA:NS:RD" << "\t" ;

  int enrichment = 0;
        
  for(unsigned int t = 0; t < localOpts.all.size(); t++){

    if(ti[localOpts.all[t]]->nBad > 2){
      enrichment = 1;
    }

    tmpOutput << ti[localOpts.all[t]]->genotype 
	      << ":" << ti[localOpts.all[t]]->gls[0]
	      << "," << ti[localOpts.all[t]]->gls[1]
	      << "," << ti[localOpts.all[t]]->gls[2]
	      << ":" << ti[localOpts.all[t]]->nGood
	      << ":" << ti[localOpts.all[t]]->nBad
	      << ":" << ti[localOpts.all[t]]->nClipping
	      << ":" << ti[localOpts.all[t]]->nReads    ;
    if(t < localOpts.all.size() - 1){
      tmpOutput << "\t";
    }
  }
  
  tmpOutput << endl;

  if(enrichment == 0 ){
    cleanUp(ti, localOpts);
    return true;
  }


  #ifdef DEBUG
  cerr << "line: " << tmpOutput.str();
  #endif 
  
  results.append(tmpOutput.str());
  
  cleanUp(ti, localOpts);
  
  delete info;
  
  return true;
}


bool filter(BamAlignment & al){

  if(!al.IsMapped()){
    return false;
  }
  if(al.IsDuplicate()){
    return false;
  }
  if(! al.IsPrimaryAlignment()
     && ((al.AlignmentFlag & 0x0800) == 0)){
    return false;
  }


  string saTag;
  if(al.GetTag("SA", saTag)){ 
  }
  else{
    if(al.MapQuality < 21){
      return false;
    }
  }

  string xaTag;
  
  if(al.GetTag("XA", xaTag)){
    vector<string> xas = split(xaTag, ";");
      if(xas.size() > 2){
	#ifdef DEBUG
	cerr << "failed xa filter" << al.Name << " " << xaTag << endl;
	#endif
	
	return false;
      }
  }

  if(checkN(al.QueryBases)){
    return false;
  }

  return true;
}
 
bool runRegion(int seqidIndex, 
	       int start, 
	       int end, 
	       vector< RefData > seqNames, 
	       vector<uint64_t> kmerDB){
  
  string regionResults;

  omp_set_lock(&lock);

  global_opts localOpts = globalOpts;
  insertDat localDists  = insertDists;

  omp_unset_lock(&lock);

  BamMultiReader All;
  
  prepBams(All, "all");

  if(!All.SetRegion(seqidIndex, start, seqidIndex, end)){
    return false;
  }


  BamAlignment al     ;
  readPileUp allPileUp;
  bool hasNextAlignment = true;

  hasNextAlignment = All.GetNextAlignment(al);
  if(!hasNextAlignment){
    All.Close();
    return false;
  }

  list <long int> clippedBuffer;
  long int currentPos  = -1;
  
  while(hasNextAlignment){    
    while(currentPos >= clippedBuffer.front()){
      if(clippedBuffer.empty()){
	break;
      }
      clippedBuffer.pop_front();
    }
    while(clippedBuffer.empty()){
      hasNextAlignment = All.GetNextAlignment(al);
      if(!hasNextAlignment){
	break;
      }
      if(!filter(al)){
	continue;
      }
      vector< CigarOp > cd = al.CigarData;
      if(cd.front().Type == 'S'){
	clippedBuffer.push_back(al.Position);
      }
      if(cd.back().Type  == 'S'){
	clippedBuffer.push_back(al.GetEndPosition(false,true));
      }
      allPileUp.processAlignment(al);
    }
    
    clippedBuffer.sort();
    
    while(al.Position <= clippedBuffer.front()){
      hasNextAlignment = All.GetNextAlignment(al);
      if(!hasNextAlignment){
        break;
      }
      if(!filter(al)){
        continue;
      }
      vector< CigarOp > cd = al.CigarData;
      if(cd.front().Type == 'S'){
        clippedBuffer.push_back(al.Position);
      }
      if(cd.back().Type  == 'S'){
        clippedBuffer.push_back(al.GetEndPosition(false,true));
      }
      allPileUp.processAlignment(al);
    }
    
    clippedBuffer.sort();

    #ifdef DEBUG
    cerr << "About to score : " << currentPos << endl;
    #endif

    allPileUp.purgePast( &currentPos );    

    if(! score(seqNames[seqidIndex].RefName, 
	       &currentPos, 
	       allPileUp,
	       localDists, 
	       regionResults, 
	       localOpts,
	       kmerDB)){
      cerr << "FATAL: problem during scoring" << endl;
      cerr << "FATAL: wham exiting"           << endl;
      exit(1);
    }

    currentPos = clippedBuffer.front();

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
  
  All.Close();

  return true;
}

bool loadKmerDB(vector<uint64_t> & DB){

  ifstream kmerDB (globalOpts.mask);
  string line;

  uint64_t kmer;

  if(kmerDB.is_open()){
    while(getline(kmerDB, line)){
      kmer = std::stoul(line);
      DB.push_back(kmer);
    }
  }
  else{
    return false;
  }

  kmerDB.close();
  return true;
}

bool loadBed(vector<regionDat*> & features, RefVector seqs){

  map<string, int> seqidToInt;

  int index = 0;

  for(vector< RefData >::iterator sit = seqs.begin(); sit != seqs.end(); sit++){

    seqidToInt[ (*sit).RefName ] = index;

    index+=1;
  }

  ifstream featureFile (globalOpts.bed);

  string line;

  if(featureFile.is_open()){

    while(getline(featureFile, line)){

      vector<string> region = split(line, "\t");

      int start = atoi(region[1].c_str()) ;
      int end   = atoi(region[2].c_str()) ;
  
      regionDat * r = new regionDat;
      r->seqidIndex = seqidToInt[region[0]];
      r->start      = start;
      r->end        = end  ;
      features.push_back(r);
    }
  }
  else{
    return false;
  }

  featureFile.close();
  
  return true;
}

int main(int argc, char** argv) {

#ifdef DEBUG
  cerr << "INFO: WHAM is in debug mode" << endl;
#endif

  omp_init_lock(&lock);

  srand((unsigned)time(NULL));

  globalOpts.nthreads = -1;

  parseOpts(argc, argv);
  
  if(globalOpts.nthreads == -1){
  }
  else{
    omp_set_num_threads(globalOpts.nthreads);
  }
 
  //loading up filenames into a vector

  globalOpts.all.reserve(globalOpts.targetBams.size()                          + globalOpts.backgroundBams.size() );
  globalOpts.all.insert( globalOpts.all.end(), globalOpts.targetBams.begin(),         globalOpts.targetBams.end() );
  globalOpts.all.insert( globalOpts.all.end(), globalOpts.backgroundBams.begin(), globalOpts.backgroundBams.end() );
  
  // loading kmer database
  vector<uint64_t> kmerDB;
  if(! (globalOpts.mask.compare("NA") == 0)){
    if(!loadKmerDB(kmerDB)){
      cerr << "FATAL: masking file was specified, but could not be opened or read." << endl;
      exit(1);
    }
  }

  cerr << "INFO: gathering stats for each bam file." << endl;
  cerr << "INFO: this step can take a few minutes." << endl;

 #pragma omp parallel for
  for(unsigned int i = 0; i < globalOpts.all.size(); i++){
    grabInsertLengths(globalOpts.all[i]);
  }

  // the pooled reader

  BamMultiReader allReader;

  // grabbing sam header and checking for sotrted bams
  prepBams(allReader, "all");
  SamHeader SH = allReader.GetHeader();
  if(!SH.HasSortOrder()){
    cerr << "FATAL: sorted bams must have the @HD SO: tag in each SAM header." << endl;
    exit(1);
  }
  allReader.Close();

  prepBams(allReader, "all");  
  RefVector sequences = allReader.GetReferenceData();
  allReader.Close();

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

  if(seqidIndex != 0 || globalOpts.region.size() == 2 ){
    if(! runRegion(seqidIndex, 
		   globalOpts.region[0], 
		   globalOpts.region[1], 
		   sequences, 
		   kmerDB)){
      cerr << "WARNING: region failed to run properly." << endl;
    }
    cerr << "INFO: WHAM-BAM finished normally." << endl;
    return 0;
  }
  
  vector< regionDat* > regions; 
  if(globalOpts.bed == "NA"){
    for(vector< RefData >::iterator sit = sequences.begin(); sit != sequences.end(); sit++){
      int start = 500;
      if((*sit).RefLength < 2000){
	cerr << "WARNING: " << (*sit).RefName << " is too short for WHAM-BAM: " << (*sit).RefLength << endl;
	continue;
      }
      
      for(;start < ( (*sit).RefLength - 500) ; start += 1000000){
	regionDat * chunk = new regionDat;
	chunk->seqidIndex = seqidIndex;
	chunk->start      = start;
	chunk->end        = start + 1000000 ;
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
  }
  else{
    loadBed(regions, sequences);
  }

 #pragma omp parallel for
  
  for(unsigned int re = 0; re < regions.size(); re++){

    omp_set_lock(&lock);
    cerr << "INFO: running region: " << sequences[regions[re]->seqidIndex].RefName << ":" << regions[re]->start << "-" << regions[re]->end << endl;
    omp_unset_lock(&lock);

    if(! runRegion( regions[re]->seqidIndex, regions[re]->start, regions[re]->end, sequences, kmerDB)){
      omp_set_lock(&lock);
      cerr << "WARNING: region failed to run properly: " 
	   << sequences[regions[re]->seqidIndex].RefName 
	   << ":"  << regions[re]->start << "-" 
	   << regions[re]->end 
	   <<  endl;
      omp_unset_lock(&lock);
    }

  }

  cerr << "INFO: WHAM-BAM finished normally." << endl;
  return 0;
}
