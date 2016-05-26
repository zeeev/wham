#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <cmath>
#include <time.h>
#include <algorithm>
#include "split.h"
#include "KMERUTILS.h"
#include "fastahack/Fasta.h"
#include "ssw_cpp.h"

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
  vector<string> targetBams      ;
  vector<string> backgroundBams  ;
  vector<string> all             ;
  int            nthreads        ;
  string         fasta           ;
  string         seqid           ;
  string         bed             ; 
  string         mask            ;
  int            qualLookup[126] ;
  int            qualCut         ; // average base pair quality for 5nt left or right of the breakpoint
  int            MQ              ; // required mapping quality for soft-clipped reads
  unsigned int   min             ; // number of soft-clips with the same start or end
  vector<int>    region          ; 
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

static const char *optString ="ht:f:b:r:x:e:m:q:p:i";


// this lookup is good for sanger and illumina 1.8+

int SangerLookup[126] =    {-1,-1,-1,-1,-1, -1,-1,-1,-1,-1, // 0-9     1-10
			    -1,-1,-1,-1,-1, -1,-1,-1,-1,-1, // 10-19   11-20
			    -1,-1,-1,-1,-1, -1,-1,-1,-1,-1, // 20-29   21-30
			    -1,-1, 0, 1, 2,  3, 4, 5, 6, 7, // 30-39   31-40
			     8, 9,10,11,12, 13,14,15,16,17, // 40-49   41-50
			    18,19,20,21,22, 23,24,25,26,27, // 50-59   51-60
			    28,29,30,31,32, 33,34,35,36,37, // 60-69   61-70
			    38,39,40,41,42, 43,44,45,46,47, // 70-79   71-80
			    48,49,50,51,52, 53,54,55,56,57, // 80-89   81-90
			    58,59,60,61,62, 63,64,65,66,67, // 90-99   91-100
			    68,69,70,71,72, 73,74,75,76,77, // 100-109 101-110
			    78,79,80,81,82, 83,84,85,86,87, // 110-119 111-120
			    88,89,90,91,92, 93           }; // 120-119 121-130

int IlluminaOneThree[126] = {-1,-1,-1,-1,-1, -1,-1,-1,-1,-1, // 0-9     1-10
                             -1,-1,-1,-1,-1, -1,-1,-1,-1,-1, // 10-19   11-20
                             -1,-1,-1,-1,-1, -1,-1,-1,-1,-1, // 20-29   21-30
                             -1,-1,-1,-1,-1, -1,-1,-1,-1,-1, // 30-39   31-40
                             -1,-1,-1,-1,-1  -1,-1,-1,-1,-1, // 40-49   41-50
                             -1,-1,-1,-1,-1, -1,-1,-1,-1,-1, // 50-59   51-60
                             -1,-1,-1,-1,-1,  0, 1, 2, 3, 4, // 60-69   61-70
                              5, 6, 7, 8, 9, 10,11,12,13,14, // 70-79   71-80
                             15,16,17,18,19, 20,21,22,23,24, // 80-89   81-90
                             25,26,27,28,29, 30,31,32,33,34, // 90-99   91-100
                             35,36,37,38,39, 40,-1,-1,-1,-1, // 100-109 101-110
                             -1,-1,-1,-1,-1, -1,-1,-1,-1,-1, // 110-119 111-120
                             -1,-1,-1,-1,-1, -1,-1        }; // 120-119 121-130


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

// rounding ints 

int RoundNum(int num)
{
  int rem = num % 10;
  return rem >= 5 ? (num - rem + 10) : (num - rem);
}

int roundUp(int numToRound, int multiple) 
{ 
  if(multiple == 0) 
    { 
      return numToRound; 
    } 

  int remainder = numToRound % multiple;
  if (remainder == 0)
    return numToRound;
  return numToRound + multiple - remainder;
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
  cout << "##fileformat=VCFv4.1"                                                                                            << endl;
  cout << "##INFO=<ID=LRT,Number=1,Type=Float,Description=\"Likelihood Ratio Test statistic\">"                             << endl;
  cout << "##INFO=<ID=WAF,Number=3,Type=Float,Description=\"Allele frequency of: background,target,combined\">"             << endl;
  cout << "##INFO=<ID=GC,Number=2,Type=Integer,Description=\"Number of called genotypes in: background,target\">"           << endl;
  cout << "##INFO=<ID=AT,Number=15,Type=Float,Description=\"Pileup attributes\">"                                           << endl;
  cout << "##INFO=<ID=CF,Number=1,Type=Float,Description=\"Fraction of reads with more than three cigar operations\">"      << endl;
  cout << "##INFO=<ID=CISTART,Number=2,Type=Integer,Description=\"PE confidence interval around POS\">"                     << endl;
  cout << "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"PE confidence interval around END\">"                       << endl;
  cout << "##INFO=<ID=PU,Number=1,Type=Integer,Description=\"Number of reads supporting position\">"                        << endl;
  cout << "##INFO=<ID=SU,Number=1,Type=Integer,Description=\"Number of supplemental reads supporting position\">"           << endl;
  cout << "##INFO=<ID=CU,Number=1,Type=Integer,Description=\"Number of neighboring all soft clip clusters across all individuals at pileup position\">" << endl;
  cout << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Number of reads at pileup position across individuals passing filters\">" << endl;
  cout << "##INFO=<ID=NC,Number=1,Type=Float,Description=\"Number of soft clipped sequences collapsed into consensus\">"   << endl;
  cout << "##INFO=<ID=MQ,Number=1,Type=Float,Description=\"Average mapping quality\">"                                      << endl;
  cout << "##INFO=<ID=MQF,Number=1,Type=Float,Description=\"Fraction of reads with MQ less than 50\">"                      << endl;
  cout << "##INFO=<ID=SP,Number=3,Type=Integer,Description=\"Number of reads supporting endpoint: mate-position,split-read,alternative-mapping\">" << endl;
  cout << "##INFO=<ID=CHR2,Number=3,Type=String,Description=\"Other seqid\">"                                               << endl;
  cout << "##INFO=<ID=DI,Number=1,Type=Character,Description=\"Consensus is from front or back of pileup: f,b\">"           << endl;
  cout << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">"      << endl;
  cout << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between POS and END\">"                << endl;
  cout << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"                                                  << endl;
  cout << "##FORMAT=<ID=GL,Number=A,Type=Float,Description=\"Genotype Likelihood\">"                                        << endl;
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
         << (*it).InsertSize   << " " 
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
  cerr << "         -The options -m,-q, and -p, control the sensitivity and specificity" << endl;
  cerr << "         -If you have exome data use the -e option for best performance     " << endl;
  cerr << endl;
}

void printHelp(void){
  cerr << "usage  : WHAM-BAM -f <STRING> -m <INT> -q <INT> -p <INT> -x <INT> -r <STRING> -e <STRING> -t <STRING> -b <STRING> " << endl << endl;
  cerr << "example: WHAM-BAM if my.fasta -m 2 -q 15 -p 10 -x 20 -r chr1:0-10000 -e genes.bed -t a.bam,b.bam -b c.bam,d.bam" << endl << endl; 

  cerr << "required   : t <STRING> -- comma separated list of target bam files          " << endl ;
  cerr << "required   : f <STRING> -- reference sequence reads were aligned to          " << endl ;

  cerr << "option     : b <STRING> -- comma separated list of background bam files      " << endl ;
  cerr << "option     : r <STRING> -- a genomic region in the format \"seqid:start-end\"" << endl ;
  cerr << "option     : x <INT>    -- set the number of threads, otherwise max [all]    " << endl ; 
  cerr << "option     : e <STRING> -- a bedfile that defines regions to score  [none]   " << endl ; 
  cerr << "option     : m <INT>    -- minimum number of soft-clips supporting           " << endl ;
  cerr << "                           START [3]                                         " << endl ;                
  cerr << "option     : q <INT>    -- exclude soft-cliped sequences with average base   " << endl ;
  cerr << "                           quality below phred scaled value (0-41) [20]      " << endl ; 
  cerr << "option     : p <INT>    -- exclude soft-clipped reads with mapping quality   " << endl ;
  cerr << "                           below value [15]                                  " << endl ; 
  cerr << "option     : i          -- base quality is Illumina 1.3+ Phred+64            " << endl ; 
  cerr << endl;
  printVersion();
}

void parseOpts(int argc, char** argv){
  int opt = 0;

  globalOpts.mask  = "NA";
  globalOpts.bed   = "NA";

  globalOpts.qualCut = 20 ;
  globalOpts.MQ      = 15 ;
  globalOpts.min     = 3  ;

  opt = getopt(argc, argv, optString);

  while(opt != -1){
    switch(opt){
    case 'i':
      {
	cerr << "INFO: base quality is Illumina 1.3+ Phred+64 NOT sanger Phred+33" << endl;
	for(unsigned int z = 0; z < 126; z++){

	  cerr <<   z << " " << SangerLookup[z] << " " << IlluminaOneThree[z] << endl;

	  SangerLookup[z] = IlluminaOneThree[z];
	}
	break;
      }
    case 'p':
      {
	globalOpts.MQ = atoi(((string)optarg).c_str());
	cerr << "INFO: WHAM-BAM skip soft-clips with mapping quaity below: " << globalOpts.MQ << endl;
	break;
      }
    case 'q':
    {
      globalOpts.qualCut = atoi(((string)optarg).c_str());
      cerr << "INFO: WHAM-BAM skip soft-clips with average base quality below: " << globalOpts.qualCut << endl;
      break; 
    }
    
    case 'f':
      {
	globalOpts.fasta =  optarg;
	cerr << "INFO: WHAM-BAM will using the following fasta: " << globalOpts.fasta << endl;
	break;
      }
      
    case 'm':
      {
	globalOpts.min = atoi(((string)optarg).c_str());
	cerr << "INFO: WHAM-BAM requires : " << globalOpts.min << " soft-clips to score breakpoint " << endl;
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

	if(tmp_region.size() != 2 || tmp_region[1].empty() || tmp_region[0].empty()){
	  cerr << "FATAL: region was not set correctly" << endl;
	  cerr << "INFO:  region format: seqid:start-end" << endl; 
	  exit(1);
	}

	vector<string> start_end = split(tmp_region[1], "-");
        globalOpts.seqid = tmp_region[0];
        globalOpts.region.push_back(atoi(start_end[0].c_str()));
        globalOpts.region.push_back(atoi(start_end[1].c_str()));


	if(start_end.size() !=2 || start_end[0].empty() || start_end[1].empty()){
	  cerr << "FATAL: region was not set correctly" << endl;
          cerr << "INFO:  region format: seqid:start-end" << endl;
          exit(1);
	}
		
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
  if(globalOpts.fasta.empty()){
    cerr << "FATAL: Failure to specify fasta file." << endl;
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

double mean(vector<int> & data){

  double sum = 0;

  for(vector<int>::iterator it = data.begin(); it != data.end(); it++){
    sum += (*it);
  }
  return sum / data.size();
}


double mean(vector<double> & data){

  double sum = 0;

  for(vector<double>::iterator it = data.begin(); it != data.end(); it++){
    sum += (*it);
  }
  return sum / data.size();
}

double var(vector<double> & data, double mu){
  double variance = 0;

  for(vector<double>::iterator it = data.begin(); it != data.end(); it++){
    variance += pow((*it) - mu,2);
  }

  return variance / (data.size() - 1);
}

// gerates per bamfile statistics 

void grabInsertLengths(string & targetfile){


  omp_set_lock(&lock);
  int quals [126];
  memcpy(quals, SangerLookup, 126*sizeof(int));
  omp_unset_lock(&lock);

  vector<double> alIns;
  vector<double> nReads;

  BamReader bamR;
  if(!bamR.Open(targetfile)   ){
    cerr << "FATAL: cannot find - or - read : " << targetfile << endl;
    exit(1);
  }

  if(! bamR.LocateIndex()){
    cerr << "FATAL: cannot find - or - open index for : " << targetfile << endl;
    exit(1);
  }

  SamHeader SH = bamR.GetHeader();
  if(!SH.HasSortOrder()){
    cerr << "FATAL: sorted bams must have the @HD SO: tag in each SAM header." << endl;
    exit(1);
  }

  RefVector sequences = bamR.GetReferenceData();

  int i = 0; // index for while loop
  int n = 0; // number of reads
  
  BamAlignment al;

  int qsum = 0;
  int qnum = 0;

  int fail = 0;


  while(i < 5 || n < 100000){

    fail += 1;
    if(fail > 1000000){
      cerr << "FATAL: was not able to gather stats on bamfile: " << targetfile << endl;
      exit(1);
    }

    unsigned int max = 20;
    
    if(sequences.size() < max){
      max = sequences.size() ;
    }
    
    int randomChr = rand() % (max -1);
    int randomPos = rand() % (sequences[randomChr].RefLength -1);
    int randomEnd = randomPos + 2000;

    if(randomEnd > sequences[randomChr].RefLength){
      continue; 
    }

    if(! bamR.SetRegion(randomChr, randomPos, randomChr, randomEnd)){      
      cerr << "FATAL: Cannot set random region";
      exit(1);
    }
        
    if(!bamR.GetNextAlignmentCore(al)){
      continue;
    }

    i++;
    
    long int cp = al.GetEndPosition(false,true);
    
    readPileUp allPileUp;
    
    while(bamR.GetNextAlignment(al)){
      if(!al.IsMapped() || ! al.IsProperPair()){
	continue;
      }
      string any;
      if(al.GetTag("XA", any)){
	continue;
      }
      if(al.GetTag("SA", any)){
	continue;
      }
      if(al.IsDuplicate()){
	continue;
      }

      string squals = al.Qualities;
      
      // summing base qualities (quals is the lookup)

      for(unsigned int q = 0 ; q < squals.size(); q++){
	qsum += quals[ int(squals[q]) ];
	qnum += 1;

	if(quals[int(squals[q])] < 0){
	  omp_set_lock(&lock);
	  cerr << endl;
          cerr << "FATAL: base quality is not sanger or illumina 1.8+ (0,41) in file : " << targetfile << endl;
          cerr << "INFO : offending qual string   : " << squals << endl;
          cerr << "INFO : offending qual char     : " << squals[q] << endl;
	  cerr << "INFO : -1 qual ; qual ; qual +1: " << quals[ int(squals[q]) -1 ] << " " << quals[ int(squals[q])  ] << " " << quals[ int(squals[q]) +1 ] << endl;
          cerr << "INFO : rescale qualities or contact author for additional quality ranges" << endl;
	  cerr << endl;
	  omp_unset_lock(&lock);
	  exit(1);
	}
      }


      if(al.Position > cp){
	allPileUp.purgePast(&cp);
	cp = al.GetEndPosition(false,true);
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
  }
  bamR.Close();

  double mu       = mean(alIns        );
  double mud      = mean(nReads       );
  double variance = var(alIns, mu     );
  double sd       = sqrt(variance     );
  double sdd      = sqrt(var(nReads, mud ));
  
  omp_set_lock(&lock);

  insertDists.mus[  targetfile ] = mu;
  insertDists.sds[  targetfile ] = sd;
  insertDists.avgD[ targetfile ] = mud;

  cerr << "INFO: for file:" << targetfile << endl            
       << "      " << targetfile << ": mean depth: ......... " << mud << endl
       << "      " << targetfile << ": sd   depth: ......... " << sdd << endl
       << "      " << targetfile << ": mean insert length: . " << insertDists.mus[targetfile] << endl
       << "      " << targetfile << ": sd   insert length: . " << insertDists.sds[targetfile] << endl
       << "      " << targetfile << ": lower insert length:  " << insertDists.mus[targetfile] -(2.5*insertDists.sds[targetfile]) << endl
       << "      " << targetfile << ": upper insert length:  " << insertDists.mus[targetfile] +(2.5*insertDists.sds[targetfile])   << endl 
       << "      " << targetfile << ": average base quality: " << double(qsum)/double(qnum) << " " << qsum << " " << qnum << endl
       << "      " << targetfile << ": number of reads used: " << n  << endl << endl;

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
    cerr << "INFO : check bam names " << endl;
    cerr << "INFO : check that the indicies are *.bam.bai "      << endl;
    cerr << "INFO : check that the ulimit is not being exceeded" << endl;
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

bool processGenotype(string    & fname        ,
		     indvDat   * idat         , 
		     double    * totalAlt     , 
		     double    * totalAltGeno ,
		     double    * relativeDepth,
		     insertDat * stats        ){

  string genotype = "./.";

  long double aal = 0;
  long double abl = 0;
  long double bbl = 0;

  if(idat->badFlag.size() < 2){
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

  double within_sample_clip_frq = double( idat->nClipping )  /  double( idat->badFlag.size() ) ;
  
  for(vector< int >::iterator rit = idat->badFlag.begin(); rit != idat->badFlag.end(); rit++){
    if(ri > 1000){
      break;
    }

    if(idat->MapQ[ri] == 0){
      idat->MapQ[ri] = 1;
    }

    double mappingP = unphred(idat->MapQ[ri]);

  #ifdef DEBUG
    cerr << "Geno: clipping must reach 0.01 frq: " << within_sample_clip_frq << endl;
  #endif

    if( (*rit == 1 && ( within_sample_clip_frq > 0.01 )) && idat->nClipping > 1){
#ifdef DEBUG
      cerr << "geno mq unphred: " << mappingP << " " << idat->MapQ[ri] << endl;
      cerr << "geno aal abl bbl: " << aal << " " << abl << " " << bbl <<  endl;
#endif
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


  #ifdef DEBUG
  cerr << "Geno: genotype likelihoods before corrections: 0/0 , 0/1, 1/1 " << aal << " " << abl << " " << bbl << endl;
  #endif

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

bool loadReadsIntoIndvs(map<string, indvDat*> & ti   , 
			readPileUp  & pileup         , 
			global_opts & localOpts      , 
			insertDat   & localDists     , 
			long int    * pos            , 
			vector<int> & insertLength   ,
			vector<int> & outOfBounds    ){    

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
    
    if((( abs( (*r).Position - *pos ) < 5 )
	|| ( abs( (*r).GetEndPosition(false,true) - *pos) < 5 ))
       && (cd.front().Type == 'S' || cd.back().Type == 'S')  ){
      
      if( ((pileup.primary[(*r).Position].size() > 1) || (pileup.primary[(*r).GetEndPosition(false,true)].size() > 1))){
	
	bad = 1;
	ti[fname]->nClipping++;
      }    
    }
    


//    if( ((pileup.primary[(*r).Position].size() > 1) || (pileup.primary[(*r).GetEndPosition(false,true)].size() > 1))
//	&& (cd.front().Type == 'S' || cd.back().Type == 'S') ){
//      bad = 1;
//      ti[fname]->nClipping++;
//    }
    
    if((*r).IsMapped() && (*r).IsMateMapped() && ((*r).RefID == (*r).MateRefID)){
      
      ti[fname]->insertSum    += abs(double((*r).InsertSize));
      ti[fname]->mappedPairs  += 1;
      
      if(( (*r).IsReverseStrand() && (*r).IsMateReverseStrand() ) || ( !(*r).IsReverseStrand() && !(*r).IsMateReverseStrand() )){
	bad = 1;
	ti[fname]->sameStrand += 1;
      }
      
      double ilength = abs ( double ( (*r).InsertSize ));
      
      double iDiff = abs ( ilength - localDists.mus[(*r).Filename] );
            
      int test = 0;

      if(iDiff > ( 2.5 * insertDists.sds[(*r).Filename] ) ){
	test = 1;

	outOfBounds.push_back(  (*r).MatePosition );
	insertLength.push_back( (*r).InsertSize  );

	bad = 1;
	ti[fname]->nAboveAvg += 1;
	ti[fname]->hInserts.push_back(ilength);
      
	if( ilength < localDists.mus[(*r).Filename]){
	  pileup.mateTooClose++;
	}
	if( ilength > localDists.mus[(*r).Filename] ){
	  pileup.mateTooFar++;
	}
      }
#ifdef DEBUG
      cerr << "IL: " << ilength << " " << iDiff << " " << 3.0 * insertDists.sds[(*r).Filename] << " " <<  test << endl;
#endif
      
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
    ti[fname]->MapQ.push_back((*r).MapQuality);
#ifdef DEBUG
    ti[fname]->alignments.push_back(*r);
#endif 
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

void endPos(vector<cigar> & cigs, int * pos){
  
  // MXDN= all consume reference

  for(vector<cigar>::iterator it = cigs.begin(); 
      it != cigs.end(); it++){
    
    switch( (*it).type ){
    case 'M':
      {
	*pos += (*it).length;
        break;
      }
    case 'X':
      {
        *pos += (*it).length;
        break;
      }
    case 'D':
      {
	*pos += (*it).length;
        break;
      }
    case '=':
      {
	*pos += (*it).length;
        break;
      }
    case 'N':
      {
	*pos += (*it).length;
        break;
      }
            
    default:
      break;
    }
  }
}


bool uniqClips(long int * pos, 
	       map<long int, vector < BamAlignment > > & clusters, 
	       vector<string> & alts, string & direction, global_opts & localOpts){

  map<string, vector<string> >  clippedSeqs;

  int bcount = 0;
  int fcount = 0;

  for( vector < BamAlignment >::iterator it = clusters[(*pos)].begin(); 
       it != clusters[(*pos)].end(); it++){
        
    vector< CigarOp > cd = (*it).CigarData;
    
    if((*it).Position == (*pos)){
      string clip = (*it).QueryBases.substr(0, cd.front().Length);
      if(clip.size() < 5){
	continue;
      }

      int start = cd.front().Length - 5;      

      if(start < 0){
	start = 0 ;
      }

      int sum = 0;
      int num = 0;

      string quals = (*it).Qualities.substr(start, 5);

      #ifdef DEBUG
      cerr << "clip seq : " << clip << endl;
      cerr << "clip qual: " << quals << endl;

      #endif

      for(unsigned int q = 0 ; q < quals.size() ; q++){
        sum += localOpts.qualLookup[ int(quals[q]) ] ;
	num += 1;
      }

      #ifdef DEBUG
      cerr << "clip avgQ: " << double(sum)/double(num) << endl;
      #endif 

      if((double(sum)/double(num)) < localOpts.qualCut){
      
#ifdef DEBUG
	cerr << "clip fail due to low qual" << endl;
#endif
	continue;
      }


      clippedSeqs["f"].push_back(clip);
      fcount += 1;
    }
    if((*it).GetEndPosition(false,true) == (*pos)){

      string clip = (*it).QueryBases.substr( (*it).Length - cd.back().Length );
      if(clip.size() < 5){
	continue;
      }

      int start = (*it).Length - cd.back().Length;

      int sum = 0;
      int num = 0;

      string quals = (*it).Qualities.substr(start, 5);


      #ifdef DEBUG
      cerr << "clip seq : " << clip << endl;
      cerr << "clip qual: " << quals << endl;

      #endif

      for(unsigned int q = 0 ; q < quals.size() ; q++){
        sum += localOpts.qualLookup[ int(quals[q]) ] ;
        num += 1;
      }

#ifdef DEBUG
      cerr << "clip avgQ: " << double(sum)/double(num) <<endl;
#endif

      if((double(sum)/double(num)) < localOpts.qualCut){

#ifdef DEBUG
	cerr << "clip fail due to low qual" << endl;
#endif
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


void gatherAlternativeAlignments(vector <int> & outBounds,
				 long int     * lb       ,
				 long int     * ub       ,
				 readPileUp   & pileup   ,
				 string       & seqid    ,
				 long int     * pos      ,
				 vector<string> & supports ){

  string xaTag;
  vector<string> XAhits;

  for(list<BamAlignment>::iterator r = pileup.currentData.begin();
      r != pileup.currentData.end(); r++){

    if((*r).Position > *pos){
      continue;
    }

    if(! (*r).IsPrimaryAlignment()){
      continue;
    }

    if((*r).GetTag("XA", xaTag)){
      vector<string> tmpHit = split(xaTag, ";");
      for(vector<string>::iterator s = tmpHit.begin();
          s != tmpHit.end(); ++s){
        XAhits.push_back(*s);
      }
    }
  }

  for(vector<string>::iterator xa = XAhits.begin();
      xa != XAhits.end(); ++xa){
    // he uses a semi-colon on the last one :(
    if((*xa).empty()){
      continue;
    }

    vector<string> chimera  = split((*xa), ",");

    // we are looking for intra chrom
    if(chimera[0].compare(seqid) != 0){
      continue;
    }
    
    // remove the strand from XA
    chimera[1].erase(0,1); 

    int chimPos = atoi( chimera[1].c_str() ) -1;

    vector<cigar> c;
    
    burnCigar(chimera[2], c);

    if(c.front().type == 'S' && c.back().type == 'S'){
      continue;
    }
    if(c.back().type == 'S'){
      endPos(c, &chimPos);
      chimPos = chimPos - 1;
    }
    supports.push_back("xa");
    outBounds.push_back( chimPos );
    
#ifdef DEBUG
    cerr << "gathering XA within CHR support cigar: " << *xa << endl;
#endif
  }
}


void gatherSplitReadSupport(vector <int> & outBounds,
			    long int     * lb       ,
			    long int     * ub       ,
			    readPileUp   & pileup   ,
			    string       & seqid    ,
			    long int     * pos      ,
			    vector<string> & supports ){

  vector<string> SAhits;
  
  for(list<BamAlignment>::iterator r = pileup.currentData.begin(); 
      r != pileup.currentData.end(); r++){
    
    if((*r).Position > *pos){
      continue;
    }

    if(! (*r).IsPrimaryAlignment()){
      continue;
    }
    string saTag;   
    if((*r).GetTag("SA", saTag)){      
      vector<string> tmpHit = split(saTag, ";");
      for(vector<string>::iterator s = tmpHit.begin();
	  s != tmpHit.end(); ++s){
	SAhits.push_back(*s);
      }
    }    
  }

  for(vector<string>::iterator sa = SAhits.begin(); 
      sa != SAhits.end(); sa++){
    // he uses a semi-colon on the last one :(
    if((*sa).empty()){
      continue;
    }

    vector<string> chimera  = split((*sa), ",");

    if(chimera[0].compare(seqid) != 0){
      continue;
    }

#ifdef DEBUG
    cerr << "gathering within CHR support cigar: " << *sa << endl;
#endif
    
    vector<cigar> c;

    burnCigar(chimera[3], c);

    int chimPos = atoi( chimera[1].c_str() ) -1;

    if(c.front().type == 'S' && c.back().type == 'S'){
      continue;
    }
    if(c.back().type == 'S'){
      endPos(c, &chimPos);
      chimPos = chimPos -1;
    } 
    supports.push_back("sr");
    outBounds.push_back( chimPos );
  }
}

bool intraChromosomeSvEnd( vector <int>      & inBounds   ,
			   vector <int>      & outOfBounds,
			   int supports []                ,
			   readPileUp        & totalDat   ,
			   string            & seqid      ,
			   string            & chr2       , 
			   long int          * pos	  ,
			   long int          & setPos     ,
			   long int          & hi         ,
			   long int          & lo         ,
			   FastaReference    & RefSeq     ,
			   string            & altSeq     ){
  
#ifdef DEBUG
  cerr << "hunting for intra chromosomal end" << endl;
  cerr << "  n out: " << outOfBounds.size()   << endl; 
  cerr << "  n inb: " << inBounds.size()      << endl; 
#endif

  vector<string> supportType;
  
  for(unsigned int st = 0; st < outOfBounds.size(); st++){
    supportType.push_back("mp");
  }

  unsigned int nz = 0;
  for(vector<int>::iterator nzit = inBounds.begin(); 
      nzit != inBounds.end(); ++nzit){

    #ifdef DEBUG
    cerr << "ILOB: " << *nzit << endl; 
    #endif 

    if(*nzit == 0){
      nz += 1;
    }
  }

#ifdef DEBUG
  cerr << "NZ: " << nz << endl;
#endif

#ifdef DEBUG
  cerr << "About to gather target zone: " << endl;
  for(vector<int>::iterator pi = outOfBounds.begin();
      pi != outOfBounds.end(); pi++){
    cerr << "OB: " << *pi << endl;
  }  

#endif

  int targetZone = *pos ;
  if(outOfBounds.size() > 1){
    targetZone = mean(outOfBounds);
  }
  
  string endChunk   = RefSeq.getSubSequence(chr2, targetZone-300, 600);

  // Declares a default Aligner
  StripedSmithWaterman::Aligner aligner;
  // Declares a default filter
  StripedSmithWaterman::Filter filter;
  // Declares an alignment that stores the result
  StripedSmithWaterman::Alignment alignment;
  // Aligns the query to the ref

  //  aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment);

  aligner.Align(altSeq.c_str(), endChunk.c_str(), endChunk.size(), filter, &alignment);

  if(abs((targetZone - 500 + alignment.ref_begin)-targetZone) < 500){
    

    
  #ifdef DEBUG
    cerr << "targetZone before switch: " << targetZone  << endl;
  #endif 
  targetZone = (targetZone - 500 + alignment.ref_begin);
  
  #ifdef DEBUG
  cerr << "targetZone after switch: " << targetZone  << endl;
  #endif

}


  long int lowerBoundOut = -1;
  long int upperBoundOut = -1;


#ifdef DEBUG
    cerr << "hunting for intra chromosomal end out of bounds" << endl;
#endif
    
    lowerBoundOut = 0 ;
    upperBoundOut = 0 ;   
    
    gatherSplitReadSupport(outOfBounds, &lowerBoundOut,      
 			   &upperBoundOut, totalDat, chr2, pos, supportType); 
    gatherAlternativeAlignments(outOfBounds,&lowerBoundOut, 
				&upperBoundOut, totalDat, chr2, pos, supportType);
    
    if(outOfBounds.empty()){
      return false;
    }
 
     
#ifdef DEBUG
    for(vector<int>::iterator test = outOfBounds.begin();
        test != outOfBounds.end(); ++test){
      cerr << "PP pos:" << *test << endl;
    }
#endif

    // loading all possible breakpoints into a map.

    map<int, int> posCounts;
   
    for(vector<int>::iterator pi = outOfBounds.begin();
	pi != outOfBounds.end(); pi++){

      if(posCounts.find(RoundNum(*pi)) == posCounts.end()){
	posCounts[RoundNum(*pi)] =  1;
      }
      else{
	posCounts[RoundNum(*pi)] += 1;
      }
    }

    map<int, int> RawPosCounts;

    for(vector<int>::iterator pi = outOfBounds.begin();
	pi != outOfBounds.end(); pi++){

      if(RawPosCounts.find(*pi) == RawPosCounts.end()){
	RawPosCounts[*pi] =  1;
      }
      else{
        RawPosCounts[*pi] += 1;
      }
    }

    #ifdef DEBUG
    cerr << "targetZone: " <<  targetZone  << endl;
    #endif 

    int bestPos  = 0; 
    int maxCount = 0;
    
    for(map<int, int>::iterator ci = posCounts.begin();
	ci != posCounts.end(); ci++){      

      #ifdef DEBUG
      cerr << "breakpoint: " << ci->first << " " << ci->second << " " << (int(targetZone) - ci->first) << endl;
      #endif
      
      if(ci->second > maxCount || (ci->second == maxCount
				  && abs(targetZone - ci->first) < abs(bestPos - targetZone)) ){

	bestPos  = ci->first;
	maxCount = ci->second; 
      }
    }
    
    if(maxCount < 2){
      for(map<int, int>::iterator ci = posCounts.begin();
	  ci != posCounts.end(); ci++){
	if(abs(ci->first - targetZone) < abs(bestPos - targetZone)){
	  bestPos = ci->first;
	}
      }
    }

    vector<double>   trimmed;

    #ifdef DEBUG
    cerr << "bestPos: " << bestPos << endl;
    #endif 
    

    if(bestPos != 0){

      map<string, int> supportStrings;
      
      supportStrings["mp"] = 0;
      supportStrings["xa"] = 0;
      supportStrings["sr"] = 0;
      
      int supportVectorIndex = 0;
      
      int tmpNewBestPos = 0;

      for(vector<int>::iterator tr = outOfBounds.begin(); 
	  tr != outOfBounds.end(); ++tr){
	
	if(abs(bestPos - tmpNewBestPos) > abs(bestPos - *tr) ){
	  tmpNewBestPos = *tr;
	}
	
	if(abs(bestPos - *tr) < 150){
	  trimmed.push_back(double(*tr));
	  supportStrings[supportType[supportVectorIndex]] += 1;
	}
	supportVectorIndex += 1;
      }

      bestPos = tmpNewBestPos;
      
      supports[0] = supportStrings["mp"];
      supports[1] = supportStrings["sr"];
      supports[2] = supportStrings["xa"];


      int bh  = 0;
      int nbp = bestPos;

      for(map<int,int>::iterator rc = RawPosCounts.begin(); rc != RawPosCounts.end(); rc++ ){

	#ifdef DEBUG
	cerr << "hone pos:" << rc->first << " " << rc->second << endl;
	#endif

	if(abs(rc->first - bestPos) > 20){
	  continue;
	}
	if(rc->second > bh){
	  bh = rc->second;
	  nbp = rc->first;
	}
      }
      
      bestPos = nbp;

      double mu = mean(trimmed);
      double sd = sqrt(var(trimmed, mu));
      
      if(trimmed.size() > 1){
	hi = bestPos + int(1.96*sd);
	lo = bestPos - int(1.96*sd);	
      }
      else{
	hi  = bestPos;
	lo  = bestPos;
      }
      setPos = int(bestPos) + 1;  
    }
    
    if(nz > trimmed.size()){
      setPos = int(*pos) +1;
      hi     = int(*pos) +1;
      lo     = int(*pos) +1;
    }
    return true;
}

void interChromosomeSvEnd(readPileUp   &     totalDat, 
			  string       &        seqid, 
			  string       &         chr2, 
			  long int     *          pos,	
			  vector<int>     & positions,
			  vector<RefData> & seqnames ){

  map<string, int> alternativeChrCounts;
  
  alternativeChrCounts[seqid] = positions.size();

  for(list<BamAlignment>::iterator it = totalDat.currentData.begin();
	 it != totalDat.currentData.end(); it++){
    
    if(! (*it).IsPrimaryAlignment()
       && (((*it).AlignmentFlag & 0x0800) == 0)){
      continue;
    }
    if(! (*it).IsMateMapped()){
      continue; 
    }
    if((*it).Position > *pos){
      continue;
    }

    string chr2seqid = seqnames[(*it).MateRefID].RefName;
    
    if(chr2seqid.compare(seqid) != 0){
      
      if(alternativeChrCounts.find(chr2seqid) != alternativeChrCounts.end()){
	alternativeChrCounts[chr2seqid] += 1;
      }
      else{
	alternativeChrCounts[chr2seqid] = 1;
      }
    }

    vector<string> SAhits;
    string saTag;

    if((*it).GetTag("SA", saTag)){
      vector<string> tmpHit = split(saTag, ";");
      for(vector<string>::iterator s = tmpHit.begin();
          s != tmpHit.end(); ++s){
        SAhits.push_back(*s);
      }
    }
  
    for(vector<string>::iterator sah = SAhits.begin();
	sah != SAhits.end(); sah++){

      if((*sah).empty()){
	continue;
      }
      vector<string> saDat = split(*sah, ",");
      if(alternativeChrCounts.find(saDat[0]) != alternativeChrCounts.end()){
	alternativeChrCounts[saDat[0]] += 1;
      }
      else{
	alternativeChrCounts[saDat[0]] = 1;
      }
    }
  
//    vector<string> XAhits;
//    string xaTag;
//
//    if((*it).GetTag("XA", xaTag)){
//      vector<string> tmpHit = split(xaTag, ";");
//      for(vector<string>::iterator s = tmpHit.begin();
//          s != tmpHit.end(); ++s){
//        XAhits.push_back(*s);
//      }
//    }
//
//    for(vector<string>::iterator xah = XAhits.begin();
//        xah != XAhits.end(); xah++){
//
//      if((*xah).empty()){
//        continue;
//      }
//      vector<string> xaDat = split(*xah, ",");
//      if(alternativeChrCounts.find(xaDat[0]) != alternativeChrCounts.end()){
//        alternativeChrCounts[xaDat[0]] += 1;
//      }
//      else{
//        alternativeChrCounts[xaDat[0]] = 1;
//      }
//    }
//  
  }
  string bestChr2;
  int Chr2Count = 0;

  for(map<string,int>::iterator ac = alternativeChrCounts.begin();
      ac != alternativeChrCounts.end(); ac++){
  
#ifdef DEBUG
    cerr << "otherseqids: " << ac->first << " " << ac->second << endl;
#endif 

  
    if(ac->second > Chr2Count){
      bestChr2 = ac->first;
      Chr2Count = ac->second;
    }
  }

  if(alternativeChrCounts[seqid] == Chr2Count){
    bestChr2 = seqid;
  }

  if(seqid.compare(bestChr2) != 0){

    chr2 = bestChr2;


    int    seqidIndx = 0;

    while(bestChr2.compare(seqnames[seqidIndx].RefName) != 0){
      seqidIndx+=1;
    }
    
    positions.clear();


    for(vector<BamAlignment>::iterator it = totalDat.primary[*pos].begin();
	it != totalDat.primary[*pos].end(); it++){

      if(! (*it).IsPrimaryAlignment()
	 && (((*it).AlignmentFlag & 0x0800) == 0)){
	continue;
      }
      if(! (*it).IsMateMapped()){
	continue;
      }

      if((*it).MateRefID != seqidIndx ){
	continue;
      }
      else{
	positions.push_back((*it).MatePosition);
      }
    }
  }
}

bool score(string seqid                 , 
	   long int         * pos       , 
	   readPileUp       & totalDat  , 
	   insertDat        & localDists, 
	   string           & results   , 
	   global_opts      & localOpts ,
	   vector<uint64_t> & kmerDB    ,
	   vector<RefData>  & seqnames  ,
	   int              & offset    ,
	   string           & RefChunk  ,
	   FastaReference   & RefSeq    ){
  
  totalDat.processPileup(pos);
  
  if(totalDat.primary[*pos].size() < localOpts.min && 
     totalDat.supplement[*pos].size() < 2
     ){
#ifdef DEBUG
    cerr << "left scoring because there were too few soft clips" << endl;
#endif
    return true;
  }

  if(totalDat.nDiscordant   == 0 
     && totalDat.nsplitRead == 0 
     && totalDat.evert      == 0 
     && totalDat.primary[*pos].size() < globalOpts.min){
#ifdef DEBUG
    cerr << "left scoring because there was no discordant, no splits, no everts" << endl;
#endif
    return true;
  }
  
//  if((double(totalDat.nLowMapQ) / double(totalDat.numberOfReads)) == 1){
//#ifdef DEBUG
//    cerr << "left scoring because mapping quality was way too low" << endl;
//#endif
//    return true;
//  }

  
  double cf = double(totalDat.tooManyCigs) / double(totalDat.numberOfReads);

  stringstream attributes;
  
  attributes << "AT="
	     << double(totalDat.nPaired)           / double(totalDat.numberOfReads)
	     << ","
	     << double(totalDat.nDiscordant)       / double(totalDat.numberOfReads)
	     << ","             
	     << double(totalDat.nMatesMissing)     / double(totalDat.numberOfReads)
	     << ","
	     << double(totalDat.nSameStrand)       / double(totalDat.numberOfReads)
	     << ","
	     << double(totalDat.nCrossChr)         / double(totalDat.numberOfReads)
	     << ","
	     << double(totalDat.nsplitRead)        / double(totalDat.numberOfReads)
	     << "," 
	     << double(totalDat.nf1SameStrand)     / double(totalDat.numberOfReads)
	     << ","
	     << double(totalDat.nf2SameStrand)     / double(totalDat.numberOfReads)
	     << ","
	     << double(totalDat.nf1f2SameStrand)   / double(totalDat.numberOfReads)
	    << ","
	     << double(totalDat.internalInsertion) / double(totalDat.numberOfReads)
	     << ","
	     << double(totalDat.internalDeletion)  / double(totalDat.numberOfReads);
  
  // finding consensus sequence 
  vector<string> alts ; // pairBreaks;
  
  string direction ;
  
  uniqClips(pos, totalDat.primary, alts, direction, localOpts);
  
  if(alts.size() < 2){
    
    #ifdef DEBUG
    cerr << "returned too few alts" << endl;
    #endif
    
    return true;
  }
  
  double nn   = 0;
  
  string altSeq = consensus(alts, &nn, direction);

  if(altSeq.size() < 10){

    
#ifdef DEBUG
    cerr << "alt size too short" << endl;
#endif
    return true;
  }
  
  if(nn / double(altSeq.size()) > 0.50 ){
    #ifdef DEBUG
    cerr << "too many N" << endl;
#endif
    return true;
  }

  map < string, indvDat*> ti;

  for(unsigned int t = 0; t < localOpts.all.size(); t++){
    indvDat * i;
    i = new indvDat;
    initIndv(i);
    ti[localOpts.all[t]] = i;
  }

  vector<int> outOfbounds  ;
  vector<int> insertLength ; 

  loadReadsIntoIndvs(ti, totalDat, 
		     localOpts, localDists, 
		     pos, insertLength, outOfbounds);
  
#ifdef DEBUG
  cerr << "loaded reads into indvs" << endl;
  cerr << "  inbound : " << outOfbounds.size() << endl;
  cerr << "  outbound: " << insertLength.size() << endl;
#endif

  string esupport=".";
  string bestEnd= ".";

  string chr2 = seqid;

  long int SVLEN              = 0;
  long int otherBreakPointPos = 0;
  long int ehigh              = 0;
  long int elow               = 0;
  long int shigh              = *pos;
  long int slow               = *pos;

  // mp,sr,sa
  int supports[3] = {0,0,0};

  interChromosomeSvEnd(totalDat, seqid, chr2, pos, outOfbounds, seqnames);

  if(intraChromosomeSvEnd(insertLength, 
			  outOfbounds, 
			  supports,
			  totalDat, 
			  seqid, 
			  chr2,
			  pos, 
			  otherBreakPointPos,
			  ehigh,
			  elow,
			  RefSeq,
			  altSeq
			  )){
    SVLEN = abs( (*pos) - otherBreakPointPos  ) ;
  
    if(seqid.compare(chr2) != 0){
      SVLEN = 0;
    }
}


  //  string localChunk = RefSeq.getSubSequence(seqid, (*pos-200),400);

  #ifdef DEBUG
  cerr << "substr: " << seqid << " " << *pos << " " << otherBreakPointPos << endl;
  cerr << "offset: " << offset - 500 << " " << RefChunk.size() << endl;
#endif

  string localChunk = RefChunk.substr(offset-200, 400);
  string endChunk   = RefSeq.getSubSequence(chr2, otherBreakPointPos-200, 400);

  #ifdef DEBUG
  cerr << "substrok: " << seqid << " " << *pos << " " << otherBreakPointPos << endl;
  #endif


  // Declares a default Aligner
  StripedSmithWaterman::Aligner aligner;
  // Declares a default filter
  StripedSmithWaterman::Filter filter;
  // Declares an alignment that stores the result
  StripedSmithWaterman::Alignment alignmentEnd ;
  StripedSmithWaterman::Alignment alignment;
  // Aligns the query to the ref

  //  aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment);

  aligner.Align(altSeq.c_str(), localChunk.c_str(), localChunk.size(), filter, &alignment);
  aligner.Align(altSeq.c_str(), endChunk.c_str(), endChunk.size(), filter, &alignmentEnd);

#ifdef DEBUG
  cerr << "local Pos SSW: "  << (*pos - 200 + alignment.ref_begin) << " "  << alignment.sw_score << endl;
  cerr << "end Pos SSW  : "  << (otherBreakPointPos - 200 + alignmentEnd.ref_begin) << " " << alignmentEnd.sw_score << endl;
#endif 




//  cerr << "===== SSW result =====" << endl;
//  cerr << "Best Smith-Waterman score:\t" << alignment.sw_score << endl
//       << "Next-best Smith-Waterman score:\t" << alignment.sw_score_next_best << endl
//       << "Reference start:\t" << alignment.ref_begin << endl
//       << "Reference end:\t" << alignment.ref_end << endl
//       << "Query start:\t" << alignment.query_begin << endl
//       << "Query end:\t" << alignment.query_end << endl
//       << "Next-best reference end:\t" << alignment.ref_end_next_best << endl
//       << "Number of mismatches:\t" << alignment.mismatches << endl
//       << "Cigar: " << alignment.cigar_string << endl
//       << "current seqid: " << seqid << endl
//       << "best seqid2  : " << chr2  << endl
//       << "best Pos2    : " << otherBreakPointPos << endl
//       << "con          : " << altSeq << endl
//       << "Pos1: " << *pos << " " << (*pos - 200 + alignment.ref_begin) << endl;
//  cerr << "======================" << endl;


  string reportRef = "N";

  if(alignment.mismatches < 5 && 
     totalDat.internalDeletion > 0 && 
     totalDat.internalDeletion > totalDat.internalInsertion){
    otherBreakPointPos = (*pos - 200 + alignment.ref_begin);
    elow  = otherBreakPointPos -10;
    ehigh = otherBreakPointPos +10;

    if(seqid.compare(chr2) == 0){
      
      if(otherBreakPointPos >= *pos){
	SVLEN = abs(*pos - otherBreakPointPos)-1;
	reportRef = RefChunk.substr(offset, SVLEN+1);
	altSeq    = RefChunk.substr(offset, 1); 
      }
      else{
	otherBreakPointPos = 0;
      }
    }
  }

  
  if(alignment.mismatches < 5 &&
     alignment.query_begin > 0 &&
     totalDat.internalInsertion > 0 &&
				    totalDat.internalDeletion  < totalDat.internalInsertion){
  otherBreakPointPos = (*pos - 200 + alignment.ref_begin);
  elow  = otherBreakPointPos -10;
  ehigh = otherBreakPointPos +10;

  if(seqid.compare(chr2) == 0){
      
      if(otherBreakPointPos >= *pos){
	SVLEN = alignment.query_begin ;
	reportRef = RefChunk.substr(offset-1, 2);
        altSeq    = reportRef + altSeq.substr(0, alignment.query_begin);
      }
      else{
        otherBreakPointPos = 0;
      }
    }
  }

//  if(seqid.compare(chr2) == 0){
//    if(*pos > otherBreakPointPos){
//
//    #ifdef DEBUG
//      cerr << "left scoring: start is upstream" << endl;
//#endif
//
//      return true;
//    }
//  }
 
  if(SVLEN > 1000000 && supports[0] == 0 && supports[1] == 0 ){
#ifdef DEBUG
    cerr << "SV too large for zero support" << endl;
    cerr << "mp: " << supports[0] << " sr: " << supports[1] << " xa: " << supports[2] << endl;
#endif

    return true;
  }

  


  if(otherBreakPointPos == 0 ){
#ifdef DEBUG
    cerr << "other breakpoint was not found" << endl;
#endif
    
    return true;
  }

  vector<double> allPrimaryClippingPos;
  
  for(std::map <long int, std::vector<BamTools::BamAlignment> >::iterator iz = totalDat.primary.begin();
	iz != totalDat.primary.end(); ++iz){
    

    if(abs(int(*pos) - int(iz->first)) > 50){
      continue;
    }

    for(unsigned int i = 0; i < iz->second.size(); ++i){
      allPrimaryClippingPos.push_back(double(iz->first));
    }
  }

  if(allPrimaryClippingPos.size() > 1){
    double mu = mean(allPrimaryClippingPos);
    double sd = sqrt(var(allPrimaryClippingPos, mu)); 
    
    shigh += int(1.96 * sd) + 1;
    slow  -= int(1.96 * sd) + 1;

  }


  if(supports[0] > 0 || supports[1] > 0 || supports[2] > 0){
    stringstream sp;
    sp << supports[0] << "," << supports[1] << "," << supports[2] ;
    esupport = sp.str();
  }


  double nAlt     = 0;
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
  

  if(altSeq.size() < 5){
    cleanUp(ti, localOpts);
    return true;
  #ifdef DEBUG
    cerr << "left scoring because alt seq is less than 5" << endl;
  #endif

  }

  if(nAlt == 0 ){
    cleanUp(ti, localOpts);

  #ifdef DEBUG
    cerr << "left scoring because no alt genotypes" << endl;
  #endif
    return true;
  }


  attributes << "," 
	     << double(totalDat.mateTooClose)      / double(totalDat.numberOfReads)
	     << ","
	     << double(totalDat.mateTooFar)        / double(totalDat.numberOfReads)
	     << ","
	     << double(totalDat.evert)             / double(totalDat.numberOfReads)
	     << ","
	     << (alternative_relative_depth_sum/nAltGeno) << ";";
  
  info_field * info = new info_field; 

  initInfo(info);

  loadInfoField(ti, info, localOpts);

  string infoToPrint = infoText(info);
  
  infoToPrint.append(attributes.str());

  stringstream tmpOutput;
  
  tmpOutput  << seqid           << "\t"  ;       // CHROM
  tmpOutput  << (*pos) +1       << "\t"  ;       // POS
  tmpOutput  << "."             << "\t"  ;       // ID
  tmpOutput  << reportRef       << "\t"  ;       // REF
  tmpOutput  << altSeq          << "\t"  ;       // ALT
  tmpOutput  << "."             << "\t"  ;       // QUAL
  tmpOutput  << "."             << "\t"  ;       // FILTER
  tmpOutput  << infoToPrint                                                          ;
  tmpOutput  << "CF=" << cf                               << ";"                     ;
  tmpOutput  << "CISTART=" << slow << "," << shigh                              << ";"                     ;
  tmpOutput  << "CIEND=" << elow << "," << ehigh                              << ";"                     ;
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
  tmpOutput  << "SP="   << esupport  << ";";
  tmpOutput  << "CHR2=" << chr2     << ";";
  tmpOutput  << "DI="   << direction << ";";
  if(otherBreakPointPos == 0  ){
    tmpOutput << "END=.;SVLEN=.\t";
  }
  else{
    tmpOutput << "END=" << otherBreakPointPos << ";" << "SVLEN=" << SVLEN << "\t";
  }
  tmpOutput  << "GT:GL:NR:NA:NS:RD" << "\t" ;

  int enrichment = 0;
        
  for(unsigned int t = 0; t < localOpts.all.size(); t++){

    if(ti[localOpts.all[t]]->nBad > 1){
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
  #ifdef DEBUG
    cerr << "left scoring: no enrichment" << endl;
  #endif

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
#ifdef DEBUG
    //    cerr << "failed not mapped " << al.Name << " " << endl;
#endif
  }
  if(al.IsDuplicate()){
    return false;
#ifdef DEBUG
    //    cerr << "failed duplicate read " << al.Name << " " << endl;
#endif

  }
  if(! al.IsPrimaryAlignment()
     && ((al.AlignmentFlag & 0x0800) == 0)){
    
    #ifdef DEBUG
    //    cerr << "failed not primary not split " << al.Name << " " << endl;
#endif


    return false;
  }

  string saTag;
  if(al.GetTag("SA", saTag)){ 
  }
  else{
    if(al.MapQuality < 1){
#ifdef DEBUG
      //      cerr << "failed mapping quality too low" << al.Name << " " << endl;
#endif
      return false;
    }
  }

  string xaTag;  
  if(al.GetTag("XA", xaTag)){
    vector<string> xas = split(xaTag, ";");
      if(xas.size() > 5){
	#ifdef DEBUG
	//	cerr << "failed xa filter" << al.Name << " " << xaTag << endl;
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
  
  if(seqNames[seqidIndex].RefLength < 2000){
    cerr << "WARNING: " << seqNames[seqidIndex].RefName << " is too short to assay: " << seqNames[seqidIndex].RefLength << endl; 
    omp_unset_lock(&lock);
    return false;
  }


  if((start - 500) < 0){
    start = 500;
  }
  if((end + 500) > seqNames[seqidIndex].RefLength){
    end = seqNames[seqidIndex].RefLength -500; 
  }

  global_opts localOpts = globalOpts;

  memcpy(localOpts.qualLookup, SangerLookup, 126*sizeof(int));

  insertDat localDists  = insertDists;

  FastaReference RefSeq;

  RefSeq.open(localOpts.fasta);
  
  string refChunk = RefSeq.getSubSequence(seqNames[seqidIndex].RefName, start-500, (end-start)+500);
  
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
  
  while(hasNextAlignment || ! clippedBuffer.empty()){    
    while(currentPos >= clippedBuffer.front()){
      if(clippedBuffer.empty()){
	break;
      }
      clippedBuffer.pop_front();
    }
    while(clippedBuffer.empty()){
      hasNextAlignment = All.GetNextAlignment(al);
#ifdef DEBUG
      //      cerr << "clipping buffer empty " << al.Name << " " << al.Position << endl;
#endif
      
      if(!hasNextAlignment){
	break;
      }
      if(!filter(al)){
	continue;
      }
      vector< CigarOp > cd = al.CigarData;
      if(cd.front().Type == 'S' && cd.front().Length > 9 && al.MapQuality  > localOpts.MQ){
	clippedBuffer.push_back(al.Position);
      }
      if(cd.back().Type  == 'S' && cd.back().Length  > 9 && al.MapQuality  > localOpts.MQ){
	clippedBuffer.push_back(al.GetEndPosition(false,true));
      }
      allPileUp.processAlignment(al);
    }
    
    clippedBuffer.sort();
    
#ifdef DEBUG
    //    cerr << "clipping buffer not empty; front: " << clippedBuffer.front() << "back: " << clippedBuffer.back() << endl;
#endif

    while(al.Position <= clippedBuffer.front()){
      hasNextAlignment = All.GetNextAlignment(al);

#ifdef DEBUG
      //      cerr << "clipping not buffer empty " << al.Name << " " << al.Position << endl;
#endif
      if(!hasNextAlignment){
        break;
      }
      if(!filter(al)){
        continue;
      }
      vector< CigarOp > cd = al.CigarData;
      if(cd.front().Type == 'S' && cd.front().Length > 9 && al.MapQuality > 5){
        clippedBuffer.push_back(al.Position);
      }
      if(cd.back().Type  == 'S' && cd.back().Length > 9 && al.MapQuality > 5){
        clippedBuffer.push_back(al.GetEndPosition(false,true));
      }
      allPileUp.processAlignment(al);
    }
    
    clippedBuffer.sort();

    #ifdef DEBUG
    cerr << "About to score : " << currentPos << endl;
    #endif

    allPileUp.purgePast( &currentPos );    

    int offset = currentPos - start + 500;

    #ifdef DEBUG
    cerr << "cp: " << currentPos << " st: " << start << " st+: " << start + 200 << " of: " << offset << endl;
    #endif

    if(! score(seqNames[seqidIndex].RefName, 
	       &currentPos, 
	       allPileUp,
	       localDists, 
	       regionResults, 
	       localOpts,
	       kmerDB,
	       seqNames,
	       offset,
	       refChunk,
	       RefSeq
	       )){
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

//bool loadKmerDB(vector<uint64_t> & DB){
//
//  ifstream kmerDB (globalOpts.mask.c_str());
//  string line;
//
//  uint64_t kmer;
//
//  if(kmerDB.is_open()){
//    while(getline(kmerDB, line)){
//      kmer = strtoull(line.c_str());
//      DB.push_back(kmer);
//    }
//  }
//  else{
//    return false;
//  }
//
//  kmerDB.close();
//  return true;
//}

bool loadBed(vector<regionDat*> & features, RefVector seqs){

  map<string, int> seqidToInt;

  int index = 0;

  for(vector< RefData >::iterator sit = seqs.begin(); sit != seqs.end(); sit++){

    seqidToInt[ (*sit).RefName ] = index;

    index+=1;
  }

  ifstream featureFile (globalOpts.bed.c_str());

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
//  if(! (globalOpts.mask.compare("NA") == 0)){
//    if(!loadKmerDB(kmerDB)){
//      cerr << "FATAL: masking file was specified, but could not be opened or read." << endl;
//      exit(1);
//    }
//  }

  cerr << "INFO: gathering stats for each bam file." << endl;



#pragma omp parallel for schedule(dynamic, 3)
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

  vector< regionDat* > regions; 


  if(globalOpts.region.size() == 2){
    for(vector< RefData >::iterator sit = sequences.begin(); sit != sequences.end(); sit++){      
      if((*sit).RefName == globalOpts.seqid){
	break;
      }
      seqidIndex += 1;
    }
  }

  if(globalOpts.region.size() == 2){

    int start = globalOpts.region[0];
    int end   = globalOpts.region[1];

    int p = start;
    int e = 0;
    for(; (p+1000000) <= end; p += 1000000){
      regionDat * regionInfo = new regionDat;
      regionInfo->seqidIndex = seqidIndex             ;
      regionInfo->start      = p                      ;
      regionInfo->end        = 1000000 + p            ;
      regions.push_back(regionInfo);
      e = p + 1000000;
    }
    if(e < end){
      regionDat * regionInfo = new regionDat;
      regionInfo->seqidIndex = seqidIndex               ;
      regionInfo->start      = p                        ;
      regionInfo->end        = end                      ;
      regions.push_back(regionInfo)                     ;
    }

  }

  if(globalOpts.bed.compare("NA") == 0 && globalOpts.region.size() != 2){
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

#pragma omp parallel for schedule(dynamic, 3)
  
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
