#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <cmath>
#include <time.h>
#include <algorithm>
#include "split.h"

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
  vector<double> gls;
  vector< BamAlignment > alignments;
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
  vector<string> targetBams    ;
  vector<string> backgroundBams;
  vector<string> all           ;
  int            nthreads      ;
  string         seqid         ;
  string         bed           ; 
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

static const char *optString ="ht:b:r:x:e:";

// this lock prevents threads from printing on top of each other

omp_lock_t lock;

bool sortStringSize(string i, string j) {return (i.size() < j.size());}

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
  cout << "##fileformat=VCFv4.1"                                                                                                                  << endl;
  cout << "##INFO=<ID=LRT,Number=1,Type=Float,Description=\"Likelihood Ratio Test Statistic\">"                                                       << endl;
  cout << "##INFO=<ID=AF,Number=3,Type=Float,Description=\"Allele frequency of: background,target,combined\">" << endl;
  cout << "##INFO=<ID=GC,Number=2,Type=Integer,Description=\"Number of called genotypes in: background,target\">"  << endl;
  cout << "##INFO=<ID=CU,Number=1,Type=Integer,Description=\"Number of neighboring soft clip clusters across all individuals at pileup position \">" << endl;
  cout << "##INFO=<ID=ED,Number=.,Type=String,Description=\"Colon separated list of potenial paired breakpoints, in the format: seqid,pos,count\">" << endl;
  cout << "##INFO=<ID=BE,Number=2,Type=String,Description=\"Best end position: chr,position,count\">"                  << endl;
  cout << "##INFO=<ID=NC,Number=1,Type=String,Description=\"Number of soft clipped sequences collapsed into consensus\">"                  << endl;
  cout << "##INFO=<ID=AT,Number=13,Type=Float,Description=\"Attributes for classification\">"                                              << endl;
  cout << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Pseudo genotype\">"                                                                 << endl;
  cout << "##FORMAT=<ID=GL,Number=A,Type=Float,Description=\"Genotype likelihood \">"                                                                 << endl;
  cout << "##FORMAT=<ID=FR,Number=1,Type=Float,Description=\"Fraction of reads with soft or hard clipping\">"                                << endl;
  cout << "##FORMAT=<ID=NR,Number=1,Type=Integer,Description=\"Number of reads supporting a SV\">"                                                      << endl;
  cout << "##FORMAT=<ID=NA,Number=1,Type=Integer,Description=\"Number of reads that do not support a SV\">"                                                      << endl;
  cout << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Number of reads with mapping quality greater than 0\">"                                << endl;
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
  for(vector< BamAlignment >::iterator it = d->alignments.begin(); it != d->alignments.end(); it++){
   
    ss   << " " << (*it).Name << " " 
         << (*it).RefID    << " " 
	 << (*it).Position << " " 
	 << (*it).GetEndPosition() << " "
	 << (*it).Position + (*it).Length << " "
         << (*it).MapQuality << " " 
         << (*it).MateRefID    << " " 
         << (*it).MatePosition << " " 
	 << collapseCigar((*it).CigarData) << " " 
         << (*it).AlignmentFlag << " " 
	 << (*it).QueryBases 
	 << endl;
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
  cerr << "Version 1.2.0 ; Zev Kronenberg; zev.kronenberg@gmail.com " << endl;
  cerr << "Github version: " << VERSION << endl;
  cerr << endl;
}

void printHelp(void){
  cerr << "usage  : WHAM-BAM -x <INT> -r <STRING>     -e <STRING>  -t <STRING>    -b <STRING>   " << endl << endl;
  cerr << "example: WHAM-BAM -x 20    -r chr1:0-10000 -e genes.bed -t a.bam,b.bam -b c.bam,d.bam" << endl << endl; 

  cerr << "required: t <STRING> -- comma separated list of target bam files"           << endl ;
  cerr << "option  : b <STRING> -- comma separated list of background bam files"       << endl ;
  cerr << "option  : r <STRING> -- a genomic region in the format \"seqid:start-end\"" << endl ;
  cerr << "option  : x <INT>    -- set the number of threads, otherwise max          " << endl ; 
  cerr << "option  : e <STRING> -- a bedfile that defines regions to score           " << endl ; 
  cerr << endl;
  printVersion();
}

void parseOpts(int argc, char** argv){
  int opt = 0;

  globalOpts.bed = "NA";

  opt = getopt(argc, argv, optString);

  while(opt != -1){
    switch(opt){
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

bool grabInsertLengths(string file){

  BamReader bamR;
  BamAlignment al;

  vector<double> alIns;

  double clipped  = 0;
  double naligned = 0;

  bamR.Open(file);

  int i = 1;
  while(i < 100000 && bamR.GetNextAlignment(al)){
    if(al.IsFirstMate() && al.IsMapped() && al.IsMateMapped() && abs(double(al.InsertSize)) < 10000){
      i += 1;
      alIns.push_back(abs(double(al.InsertSize)));
    }
    if(al.IsMapped()){
      naligned += 1;
      vector< CigarOp > cd = al.CigarData;

      if(cd.back().Type == 'S' || cd.back().Type == 'H' ){
        clipped += double(cd.back().Length) / double (al.Length);
      }
      if(cd.front().Type == 'S' || cd.front().Type == 'H'){
	clipped += double(cd.front().Length) / double (al.Length);
      }

    }
  }

  bamR.Close();

  double mu       = mean(alIns        );
  double variance = var(alIns, mu     );
  double sd       = sqrt(variance     );

  insertDists.mus[file] = mu;
  insertDists.sds[file] = sd;

  cerr << "INFO: mean insert length, fragment standard deviation, number of reads, file : "
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

  bool flag = true;

  for(unsigned int i = 0; i < globalOpts.all.size(); i++){
    flag = grabInsertLengths(globalOpts.all[i]);
  }
  
  return true;
  
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


bool processGenotype(indvDat * idat, double * totalAlt){

  string genotype = "./.";

  double aal = 0;
  double abl = 0;
  double bbl = 0;

  if(idat->badFlag.size() < 3){
    idat->gls.push_back(-255.0);
    idat->gls.push_back(-255.0);
    idat->gls.push_back(-255.0);
    return true;
  }

  double nref = 0.0;
  double nalt = 0.0;

  int ri = 0;

  for(map< string, int >::iterator rit = idat->badFlag.begin(); rit != idat->badFlag.end(); rit++){

    double mappingP = unphred(idat->MapQ[ri]);

    if( idat->badFlag[rit->first] == 1 ){
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

  aal = aal - log(pow(2,idat->nReads));
  abl = abl - log(pow(2,idat->nReads));
  bbl = bbl - log(pow(2,idat->nReads));

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

  if(n > 4){
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


    if( (pileup.allCount[(*r).Position] > 1) || (pileup.allCount[(*r).GetEndPosition()] > 1)){
      bad = 1;
      ti[fname]->nClipping++;
    }
    
    if((*r).IsMapped() && (*r).IsMateMapped() && (*r).IsProperPair()){

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
    ti[fname]->badFlag[(*r).Name] = bad;
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
  
  return true;

}


string infoText(info_field * info){
  
  stringstream ss;

  ss << "LRT=" << info->lrt << ";";
  ss << "AF="  << info->taf << "," << info->baf << "," << info->aaf << ";";
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



int otherBreak(long int * pos,
	       map<long int, vector < BamAlignment > > & supliment, 
	       string & otherside,
	       string & bestEnd,
	       string & bestSeqid,
	       int * support,
	       long int * otherPos
	      
	       ){
  
#ifdef DEBUG
  cerr << "N secondary:" << supliment.size() << endl;
#endif

  if(supliment[*pos].empty()){
    return 0;
  }

  map<string, map< long int, int > > otherPositions;

#ifdef DEBUG
  cerr << "Chimeric mapping:" << endl;
#endif

  for(vector<BamAlignment>::iterator it = supliment[*pos].begin(); 
      it != supliment[*pos].end(); it++){

    string saTag;
    if(! (*it).GetTag("SA", saTag)){
      cerr << "no sa\n";
      return false;
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
	otherPositions[chimera[0]][atoi(chimera[1].c_str())]++;
      }

      if(c.back().type  == 'S'){
	long int endP = atoi(chimera[1].c_str());
	endPos(c, &endP);
	otherPositions[chimera[0]][endP]++;	
      } 
    }
  }
    
  stringstream ends;
  string bestOpt; 
  int    nbe = 0;

  for(map< string, map <long int, int> >::iterator ab = otherPositions.begin();
      ab != otherPositions.end(); ab++){

    for(map<long int, int>::iterator pos = ab->second.begin();
	pos != ab->second.end(); pos++){

      ends << ab->first << "," << pos->first << "," << pos->second << ":";

      if(pos->second > nbe){
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

  bestEnd.append(bestOpt);
  
  otherside.append(ends.str());



  return otherPositions.size();

}

bool uniqClips(long int * pos, 
	       map<long int, vector < BamAlignment > > & clusters, 
	       vector<string> & alts){

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
      if(clip.size() < 4){
	continue;
      }
      clippedSeqs["f"].push_back(clip);
      fcount += 1;
    }
    if((*it).GetEndPosition() == (*pos)){
      string clip = (*it).QueryBases.substr( (*it).Length - cd.back().Length );
      if(clip.size() < 4){
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

  for(vector<string>::iterator seqs = clippedSeqs[key].begin();
      seqs != clippedSeqs[key].end(); seqs++
      ){
    alts.push_back(*seqs);
  }

  return true;
}



string consensus(vector<string> & s, double * nn){

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

    for(int clips = s.size() -1; clips > -1; clips--){
      appendValue(seq, s[clips]);          
      index+=1;
      if(index > 19){
	break;
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

bool score(string seqid, 
	   long int * pos, 
	   readPileUp & totalDat, 
	   insertDat & localDists, 
	   string & results, 
	   global_opts localOpts){


  totalDat.processPileup(pos);

  if(double(totalDat.nPaired) / double(totalDat.numberOfReads) > 0.9999){
    return true;
  }

  if(totalDat.primaryCount[*pos] < 2 && totalDat.supplementCount[*pos] < 2){
    return true;
  }

  #ifdef DEBUG
  cerr << "Passed Cluster filters: " << totalDat.primaryCount[*pos] << " " << totalDat.supplementCount[*pos] << endl;
  #endif

  string ends, bestEnd, bestSeqid;

  int otherBreakPointCount    = 0;
  long int otherBreakPointPos = 0;

  int otherSeqids = otherBreak(pos, totalDat.supplement, ends, bestEnd, bestSeqid, &otherBreakPointCount, &otherBreakPointPos);

  if(otherSeqids > 3){
    return true;
  }
  if(seqid.compare(bestSeqid) != 0 && ! bestSeqid.empty()){
    if(otherBreakPointCount < 2){
      return true;
    }
  }

  if(seqid.compare(bestSeqid) == 0 && ! bestSeqid.empty()){
    if(abs(*pos - otherBreakPointPos) > 1000000 && otherBreakPointCount < 2){
      return true;
    }
  }
  
  

  if(ends.empty()){
    ends = "nan";
  }
  if(bestEnd.empty()){
    bestEnd = "nan";
  }
   
  vector<string> alts ; // pairBreaks;

  uniqClips(pos, totalDat.primary, alts);

  if(alts.size() < 1){
    return true;
  }

  double nn   = 0;

  string altSeq = consensus(alts, &nn);

  if(altSeq.size() < 11){
    return true;
  }

  if(nn / double(altSeq.size()) > 0.30){
    return true;
  }

  map < string, indvDat*> ti;

  for(unsigned int t = 0; t < localOpts.all.size(); t++){
    indvDat * i;
    i = new indvDat;
    initIndv(i);
    ti[localOpts.all[t]] = i;
  }

  loadIndv(ti, totalDat, localOpts, localDists, pos);

  double nAlt = 0;

  for(unsigned int t = 0; t < localOpts.all.size(); t++){
    processGenotype(ti[localOpts.all[t]], &nAlt);
    #ifdef DEBUG
    cerr << "position: " << *pos << endl; 
    cerr << printIndvDat(ti[localOpts.all[t]]) << endl;
    #endif 
  }
  
  if(nAlt == 0 ){
    cleanUp(ti, localOpts);
    return true;
  }
  
  info_field * info = new info_field; 

  initInfo(info);

  loadInfoField(ti, info, localOpts);

  string infoToPrint = infoText(info);
  
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
	     << ";";

  infoToPrint.append(attributes.str());


  stringstream tmpOutput;

  tmpOutput  << seqid           << "\t"  ;       // CHROM
  tmpOutput  << (*pos)          << "\t"  ;       // POS
  tmpOutput  << "."             << "\t"  ;       // ID
  tmpOutput  << "N"             << "\t" ;       // REF
  tmpOutput  << altSeq          << "\t"  ;       // ALT
  tmpOutput  << "."             << "\t"  ;       // QUAL
  tmpOutput  << "."             << "\t"  ;       // FILTER
  tmpOutput  << infoToPrint << ""  ;
  tmpOutput  << "CU=" << totalDat.allCount.size() << ";"  ;
  tmpOutput  << "NC=" << alts.size()     << ";"  ;
  tmpOutput  << "ED=" << ends   << ";";
  tmpOutput  << "BE=" << bestEnd << "\t";
  tmpOutput  << "GT:GL:NR:NA:DP:FR" << "\t" ;
        
  for(unsigned int t = 0; t < localOpts.all.size(); t++){

    double fr = 0;
    if(ti[localOpts.all[t]]->nClipping > 0 && ti[localOpts.all[t]]->nReads > 0){
      fr = double(ti[localOpts.all[t]]->nClipping) / double(ti[localOpts.all[t]]->nReads);
    }

    tmpOutput << ti[localOpts.all[t]]->genotype 
	      << ":" << ti[localOpts.all[t]]->gls[0]
	      << "," << ti[localOpts.all[t]]->gls[1]
	      << "," << ti[localOpts.all[t]]->gls[2]
	      << ":" << ti[localOpts.all[t]]->nGood
	      << ":" << ti[localOpts.all[t]]->nBad
	      << ":" << ti[localOpts.all[t]]->nReads
              << ":" << fr ;
    if(t < localOpts.all.size() - 1){
      tmpOutput << "\t";
    }
  }
  
  tmpOutput << endl;

  #ifdef DEBUG
  cerr << "line: " << tmpOutput.str();
  #endif 
  
  results.append(tmpOutput.str());
  
  cleanUp(ti, localOpts);
  
  delete info;
  
  return true;
}


bool filt(BamAlignment & al){

  if(!al.IsMapped()){
    return false;
  }
  if(al.MapQuality == 0){
    return false;
  }
  if(al.IsDuplicate()){
    return false;
  }
  if(! al.IsPrimaryAlignment()
     && ((al.AlignmentFlag & 0x0800) == 0)){
    return false;
  }

  string xaTag;
  
  if(al.GetTag("XA", xaTag)){
    vector<string> xas = split(xaTag, ";");
      if(xas.size() > 2){
	return false;
      }
  }

  if(checkN(al.QueryBases)){
    return false;
  }


  return true;
}
 
bool runRegion(int seqidIndex, int start, int end, vector< RefData > seqNames){
  
  string regionResults;

  omp_set_lock(&lock);

  global_opts localOpts = globalOpts;
  insertDat localDists = insertDists;

  omp_unset_lock(&lock);

  BamMultiReader All;
  
  prepBams(All, "all");

  if(!All.SetRegion(seqidIndex, start, seqidIndex, end)){
    return false;
  }

  BamAlignment al     ;
  readPileUp allPileUp;

  if(! All.GetNextAlignment(al)){
    All.Close();
    return false;
  }

  long int currentPos = 0;
  bool getNextAl      = true;

  while(1){  
  
    bool clipped = false;
    
    if(getNextAl == false){
      break;
    }
    
    while(clipped == false && getNextAl){
      
      getNextAl = All.GetNextAlignment(al);
      
      if(filt(al)){
	
	allPileUp.processAlignment(al, currentPos);
	
	vector< CigarOp > cd = al.CigarData;
	
	if(cd.front().Type == 'S' ){
	  currentPos = al.Position;
	  clipped = true;
	}

	if(cd.back().Type == 'S' && ! clipped){
          currentPos = al.GetEndPosition();
          clipped = true;
	}

       	while(al.Position <= currentPos && getNextAl && clipped){
	  getNextAl = All.GetNextAlignment(al);
	  
	  if(al.IsMapped() &&  al.MapQuality > 0 && ! al.IsDuplicate() && getNextAl){
	    
	    allPileUp.processAlignment(al, currentPos);
	  }
	}
      }	
    }
    
    allPileUp.purgePast();    

    if(! score(seqNames[seqidIndex].RefName, 
	       &currentPos, 
	       allPileUp,
	       localDists, 
	       regionResults, 
	       localOpts )){
      cerr << "FATAL: problem during scoring" << endl;
      cerr << "FATAL: wham exiting"           << endl;
      exit(1);
    }
    
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
 
  globalOpts.all.reserve(globalOpts.targetBams.size()  + globalOpts.backgroundBams.size());
  globalOpts.all.insert( globalOpts.all.end(), globalOpts.targetBams.begin(), globalOpts.targetBams.end() );
  globalOpts.all.insert( globalOpts.all.end(), globalOpts.backgroundBams.begin(), globalOpts.backgroundBams.end() );

  BamMultiReader allReader;
  

  // checking for bam 
  prepBams(allReader, "all");
  SamHeader SH = allReader.GetHeader();
  if(!SH.HasSortOrder()){
    cerr << "FATAL: sorted bams must have the @HD SO: tag in each SAM header." << endl;
    exit(1);
  }
  


  allReader.Close();

  if(!getInsertDists()){
    cerr << "FATAL: " << "problem while generating insert lengths dists" << endl;
    exit(1);
  }

  cerr << "INFO: generated distributions" << endl;
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
    if(! runRegion(seqidIndex, globalOpts.region[0], globalOpts.region[1], sequences)){
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
