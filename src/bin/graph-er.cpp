/*
This program was created at:  Thu May  7 12:10:40 2015
This program was created by:  Zev N. Kronenberg

Contact: zev.kronenber@gmail.com

Organization: Unviersity of Utah
    School of Medicine
    Salt Lake City, Utah

The MIT License (MIT)

Copyright (c) <2015> <Zev N. Kronenberg>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*/

#include <string>
#include <iostream>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <unistd.h>

#include "split.h"
#include "fastahack/Fasta.h"
#include "ssw_cpp.h"

// openMP - swing that hammer
#include <omp.h>

// bamtools and my headers
#include "api/BamMultiReader.h"
#include "readPileUp.h"

// gsl header
#include "gauss.h"

// paired-like hmm
#include "phredUtils.h"
#include "alignHMM.h"


#define BAM_CIGAR_SHIFT 4
#define BAM_CIGAR_MASK  ((1 << BAM_CIGAR_SHIFT) - 1)

using namespace std;
using namespace BamTools;

struct options{
  std::vector<string> targetBams;
  bool statsOnly         ;  
  bool skipGeno          ;
  int MQ                 ;
  int nthreads           ;
  string fasta           ;
  string graphOut        ;
  map<string, int> toSkip;
  string seqid           ;
  vector<int> region     ;
  string svs             ; 

}globalOpts;

struct regionDat{
  int seqidIndex ;
  int start      ;
  int end        ;
};

struct readPair{
  int          flag;
  int          count;
  BamAlignment al1;
  BamAlignment al2;  
};


struct cigDat{
  int  Length;
  char Type;
};

struct saTag{
  int seqid;
  int pos;
  bool strand;
  vector<cigDat> cig;
};

struct node;
struct edge;

struct edge{
  node * L;
  node * R;
  map<char,int> support;
};

struct node{
  int   seqid          ;
  int   pos            ;
  int   endSupport     ;
  int   beginSupport   ; 
  bool  collapsed      ;
  vector <edge *> eds  ;
  map<string,int> sm   ;
};

struct graph{
  map< int, map<int, node *> > nodes;
  vector<edge *>   edges;
}globalGraph;

struct insertDat{
  map<string, double> mus ; // mean of insert length for each indvdual across 1e6 reads
  map<string, double> sds ;  // standard deviation
  map<string, double> lq  ;  // 25% of data
  map<string, double> up  ;  // 75% of the data
  map<string, double> swm ;
  map<string, double> sws ;
  map<string, double> avgD;
  double overallDepth;
} insertDists;

// options

struct breakpoints{
  bool fail              ;
  bool two               ;
  char type              ;
  int seqidIndexL        ;
  int seqidIndexR        ;
  string seqid           ;
  int merged             ;
  int refined            ;
  int five               ;
  int three              ;
  int svlen              ;
  int collapsed          ; 
  int totalSupport       ; 
  string id              ; 
  vector<string> alleles ;
  vector<int>    supports;
  vector<vector<double> > genotypeLikelhoods ;
  vector<int>             genotypeIndex      ;
  vector<int>             nref               ;
  vector<int>             nalt               ;
  vector<string>          sml                ; 
  vector<string>          smr                ; 

  double lref;
  double lalt;

};

static const char *optString = "b:m:r:a:g:x:f:e:hsk";



// omp lock

omp_lock_t lock;
omp_lock_t glock;


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : p = prob successs, q prob fail, m = success, n = failure;

 Function does   : adds soft clip length to end of read to see if it overlaps
                   end pos

 Function returns: bool

*/

bool loadExternal(vector<breakpoints *> & br, map<string, int> & inverse_lookup){

  ifstream featureFile (globalOpts.svs.c_str());
  
  string line;
  
  if(featureFile.is_open()){
    
    while(getline(featureFile, line)){
      
      vector<string> SV   = split(line, "\t");
      
      vector<string> info = split(SV[10], ";");

      map<string, string> infoDat;

      for(vector<string>::iterator iz = info.begin(); iz != info.end(); iz++){
	
	vector<string> kv = split(*iz, "=");
	
	infoDat[kv.front()] = kv.back();

      }

	
      if( (infoDat["SVTYPE"].compare("DEL") == 0) || (infoDat["SVTYPE"].compare("DUP") == 0) || (infoDat["SVTYPE"].compare("INV") == 0) ){

	breakpoints * bk = new breakpoints;
	
	bk->fail = false;
	bk->two  = true ;
	bk->type = 'D';
	
	if((infoDat["SVTYPE"].compare("DUP") == 0)){
	  bk->type = 'U';
	}
	if((infoDat["SVTYPE"].compare("INV") == 0)){
          bk->type = 'V';
        }
	if(inverse_lookup.find(SV[0]) == inverse_lookup.end() ){
	  cerr << "FATAL: could not find seqid in inverse lookup: " << SV[0] << endl;
	  exit(1);
	}
	
	vector<string> POS = split(infoDat["POS"], ",");
	vector<string> SUP = split(infoDat["SUPPORT"], ",");

	bk->seqidIndexL = inverse_lookup[SV[0]];
	bk->seqidIndexR = inverse_lookup[SV[0]];
	bk->seqid       = SV[0];
	bk->merged      = false;
	bk->refined     = false; 

	cerr << POS[0] << " " << POS[1] << endl;

	bk->five        = atoi(POS[0].c_str());
	bk->three       = atoi(POS[1].c_str());

	bk->id           = infoDat["ID"];
	bk->collapsed    = atoi(infoDat["COLLAPSED"].c_str());
	bk->lalt = 0;
	bk->lref = 0;

	bk->sml = split(infoDat["LID"], ",");
	bk->smr = split(infoDat["RID"], ",");

	bk->supports.push_back(atoi(SUP[0].c_str()));
	bk->supports.push_back(atoi(SUP[0].c_str()));
	  

	if(bk->five > bk->three){
	  cerr << "FATAL: SV starts before it ends: " << line << endl;
	  exit(1);
	}

	bk->svlen = (bk->three - bk->five);

	
	br.push_back(bk);

      }
      else{
	cerr << "FATAL: loading external breakpoints: unknown type: " << SV[3] << endl;
	exit(1);
      }
      
    }
  }
  featureFile.close();
  return true;
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : p = prob successs, q prob fail, m = success, n = failure;

 Function does   : adds soft clip length to end of read to see if it overlaps
                   end pos

 Function returns: bool

*/


double ldbinomial(double p, double q, int m, int n)
{ 
  double temp = lgamma(m + n + 1.0);
  temp -=  lgamma(n + 1.0) + lgamma(m + 1.0);
  temp += m*log(p) + n*log(q);
  return temp;
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : bam alignemnt, pos

 Function does   : adds soft clip length to end of read to see if it overlaps
                   end pos

 Function returns: bool

*/

bool endsBefore(BamAlignment & r, int & pos, int b){
  if(!r.IsMapped()){
    return false;
  }
  int end = r.GetEndPosition();

  if(r.CigarData.back().Type == 'S'){
    end  += r.CigarData.back().Length;
  }
  if((end) < pos+b){
    return true;
  }
  else{
    return false;
  }
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : bam alignemnt, pos

 Function does   : substracts soft clip length from start to see if it overlaps
                   start pos

 Function returns: bool

*/

bool startsAfter(BamAlignment & r, int & pos, int b){
  if(!r.IsMapped()){
    return false;
  }
  int start = r.Position;

  if(r.CigarData.front().Type == 'S'){
    start  -= r.CigarData.front().Length;
  }
  if((start+b) > pos){
    return true;
  }
  else{
    return false;
  }
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of strings and separator

 Function does   : joins vector with separator

 Function returns: string

*/

string join(vector<string> & strings, string sep){

  string joined = *(strings.begin());

  for(vector<string>::iterator sit = strings.begin()+1; sit != strings.end();
      sit++){

    joined = joined + sep + (*sit) ;
  }
  return joined;
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of doubles and a seperator

 Function does   : joins vector with separator

 Function returns: string

*/

string join(vector<double> & ints, string sep){

  stringstream ss;

  for(vector<double>::iterator sit = ints.begin(); sit != ints.end(); sit++){
    ss << *sit << sep;
  }
  return ss.str();
}



double totalAlignmentScore(vector<BamAlignment> & reads, breakpoints * br){

  int sum = 0;
  int n   = 0;

  // Declares a default Aligner
  StripedSmithWaterman::Aligner aligner;
  // Declares a default filter
  StripedSmithWaterman::Filter filter;
  // Declares an alignment that stores the result
  StripedSmithWaterman::Alignment alignment;
  // Aligns the query to the ref

  for(vector<BamAlignment>::iterator r = reads.begin(); 
      r != reads.end(); r++){
    
    if(endsBefore((*r), br->five,10) || startsAfter((*r), br->three,10)){
      continue;
    }

    
    if((*r).IsMapped()){
      if(((*r).CigarData.front().Type != 'S' && (*r).CigarData.front().Length < 10) &&
	 ((*r).CigarData.back().Type  != 'S' && (*r).CigarData.back().Length  < 10)){
        continue;
      }

    }

    n += 1;

    aligner.Align((*r).QueryBases.c_str(), br->alleles.back().c_str(),
                  br->alleles.back().size(),  filter, &alignment);


    sum +=  alignment.sw_score;



  }
  if(sum > 0){
    return double(sum) / double(n);
  }
  else{
    return 0;
  }
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : breakpoint pointer

 Function does   : finds alignment alt score

 Function returns: double

*/


bool getPopAlignments(vector<string> & bamFiles, breakpoints * br, vector<BamAlignment> & reads, int buffer){

  omp_set_lock(&lock);

  BamMultiReader bamR;
  if(!bamR.Open(bamFiles)){
    cerr << "FATAL: count not open bamfiles " << join(bamFiles, " ");
  }
  if(! bamR.LocateIndexes()){
    vector<string> modifiedIndexName;
    for(vector<string>::iterator files = bamFiles.begin();
	files != bamFiles.end(); files++){
      vector<string> fileName = split(*files, ".");
      fileName.back() = "bai";
      string indexName = join(fileName, ".");
      modifiedIndexName.push_back(indexName);
    }
    if(! bamR.OpenIndexes(modifiedIndexName) ){
      cerr << "FATAL: cannot find bam index." << endl;
      exit(1);
    }
  }
  omp_unset_lock(&lock);

  if(!bamR.SetRegion(br->seqidIndexL, br->five -buffer, br->seqidIndexL, br->five +buffer)){
    cerr << "FATAL: cannot set region for breakpoint refinement." << endl;
    exit(1);
  }

  BamAlignment al;

  while(bamR.GetNextAlignment(al)){
    if((al.AlignmentFlag & 0x0800) != 0 ){
      continue;
    }
    if(! al.IsPaired()){
      continue;
    }
    if(! al.IsMapped() && ! al.IsMateMapped()){
      continue;
    }
    if(al.IsDuplicate()){
      continue;
    }
    if(! al.IsPrimaryAlignment()){
      continue;
    }
    if(al.IsMapped() && al.MapQuality < globalOpts.MQ){
      continue;
    }

    reads.push_back(al);
  }
  if(br->two){
    if(!bamR.SetRegion(br->seqidIndexL, br->three-buffer, br->seqidIndexL, br->three+buffer)){
      cerr << "FATAL: cannot set region for genotyping." << endl;
      exit(1);
    }

    while(bamR.GetNextAlignment(al)){
      if((al.AlignmentFlag & 0x0800) != 0 ){
        continue;
      }
      if(! al.IsPaired()){
        continue;
      }
      if(! al.IsMapped() && ! al.IsMateMapped()){
        continue;
      }
      if(al.IsDuplicate()){
        continue;
      }
      if(! al.IsPrimaryAlignment()){
        continue;
      }
      if(al.IsMapped() && al.MapQuality < globalOpts.MQ){
	continue;
      }

      reads.push_back(al);
    }
  }

  bamR.Close();
  
  return true;
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : string reference

 Function does   : revcomp

 Function returns: string

*/
void Comp(string & seq){
  
  locale loc;

  for (size_t i = 0; i < seq.size(); ++i)
    {
      switch (toupper(seq[i], loc))
	{
	case 'A':
	  {
	    seq[i] = 'T';
	    break;
	  }
	case 'T':
	  {
	    seq[i] = 'A';
	    break;
	  }
	case 'G':
	  {
	    seq[i] = 'C';
	    break;
	  }
	case 'C':
	  {
	    seq[i] = 'G';
	    break;
	  }
	default:
	  {
	    seq[i] = 'N';
	    break;
	  }
	}
    } 
}


int phred(double v){

  double s = (-10 * log10(v) ); 
  
  if(int(s) == 0){
    return 1;
  }
  else{
    return int(s);
  }

}

double unPhred(int v){
  
  return pow(10.0,(-1*double(v)/10));

}



//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : nothing

 Function does   : prints help

 Function returns: NA

*/

void printVersion(void){
  cerr << "Version: " << VERSION << endl;
  cerr << "Contact: zev.kronenberg [at] gmail.com " << endl;
  cerr << "Notes  : -If you find a bug, please open a report on github!" << endl;
  cerr << endl;
}


void printHelp(void){
//------------------------------- XXXXXXXXXX --------------------------------  
  cerr << " Usage:  " << endl;
  cerr << "       WHAM-GRAPHENING -f my.bam -a my.fasta \\ " << endl;
  cerr << "       -g my.graph.out.txt -e M,GL000207.1 2> wham.err > wham.out" << endl;
  cerr << endl;
  cerr << " Required:  " << endl;
//------------------------------- XXXXXXXXXX --------------------------------  

  cerr << "          -f - <STRING> - A sorted and indexed bam file or a list" << endl;
  cerr << "                          of bams: a.bam,b.bam,..." << endl;
  cerr << "          -a - <STRING> - The reference genome (indexed fasta).  " << endl;
  cerr << endl;
  cerr << " Optional:  " << endl;
  cerr << "          -s - <FLAG>   - Exits the program after the stats are      " << endl;
  cerr << "                          gathered.                                  " << endl;
  cerr << "          -k - <FLAG>   - Skip genotyping (much faster).             " << endl;
  cerr << "          -g - <STRING> - File to write graph to (very large output)." << endl;
  cerr << "          -e - <STRING> - Comma sep. list of seqids to skip.         " << endl;
  cerr << "          -r - <STRING> - Region in format: seqid:start-end          " << endl;
  cerr << "          -x - <INT>    - Number of CPUs to use [default: all cores]." << endl;
  cerr << "          -m - <INT>    - Mapping quality filter [20].               " << endl;
  cerr << "          -b - <STRING> - external file to genotype.                 " << endl;

  cerr << endl;
  cerr << " Output:  " << endl;
  cerr << "        STDERR: Run statistics and bam stats                        " << endl;    
  cerr << "        STOUT : SV calls in BEDPE format (VCF soon)                 " << endl;  
  cerr << endl;
  printVersion();
}


int  getSupport(node * N){
  
  int support = 0;

  for(vector<edge *>::iterator ed = N->eds.begin() ;
      ed != N->eds.end(); ed++){
    support += (*ed)->support['L'];
    support += (*ed)->support['H'];
    support += (*ed)->support['S'];
    support += (*ed)->support['I'];
    support += (*ed)->support['D'];
    support += (*ed)->support['V'];
    support += (*ed)->support['M'];
    support += (*ed)->support['R'];
    support += (*ed)->support['X'];
  
}
  return support;
}

bool sortNodesByPos(node * L, node * R){
  return (L->pos < R->pos);
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of breakpoint calls

 Function does   : return true if the left is less than the right

 Function returns: bool

*/

bool sortNodesBySupport(node * L, node * R){

  int totalSL = 0;
  int totalSR = 0;

  for(vector<edge *>::iterator ed = L->eds.begin() ;
      ed != L->eds.end(); ed++){
      totalSL += (*ed)->support['L'];
      totalSL += (*ed)->support['H'];
      totalSL += (*ed)->support['S'];
      totalSL += (*ed)->support['I'];
      totalSL += (*ed)->support['D'];
      totalSL += (*ed)->support['V'];
      totalSL += (*ed)->support['M'];
      totalSL += (*ed)->support['R'];
      totalSL += (*ed)->support['X'];
  }
  
  
  for(vector<edge *>::iterator ed = R->eds.begin() ;
      ed != R->eds.end(); ed++){
    totalSR += (*ed)->support['L'];
    totalSR += (*ed)->support['H'];
    totalSR += (*ed)->support['S'];
    totalSR += (*ed)->support['I'];
    totalSR += (*ed)->support['D'];
    totalSR += (*ed)->support['V'];
    totalSR += (*ed)->support['M'];
    totalSR += (*ed)->support['R'];
    totalSR += (*ed)->support['X'];
  }


  return (totalSL > totalSR) ;  
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of breakpoint calls

 Function does   : return true if the left is less than the right

 Function returns: bool

*/

bool sortBreak(breakpoints * L, breakpoints * R){
  
  if(L->seqidIndexL == R->seqidIndexL ){
    if(L->five < R->five){

      return true;
    }
  }
  else{
    if(L->seqidIndexL < R->seqidIndexL){

      return true;
    }
  }
  
  return false;
  

}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of pointers to breakpoints and a bamtools RefVector

 Function does   : prints a bedbe format

 Function returns: nada

*/


void printBEDPE(vector<breakpoints *> & calls, RefVector & seqs){

  int index = 0;

   sort(calls.begin(), calls.end(), sortBreak);


  for(vector<breakpoints *>::iterator c = calls.begin(); c != calls.end(); c++){
    
    int svlen = (*c)->svlen;

    if(svlen < 5){
      (*c)->fail = true;
    }

    if((*c)->fail){
      continue;
    }

    index += 1;
   

    stringstream ss;

    string type = "NONE";
    switch((*c)->type){
    case 'D':
      type = "DEL";
      svlen = -svlen;
      break;
    case 'U':
      type = "DUP";
      break;
    case 'I':
      type = "INR";
      break;
    case 'V':
      type = "INV";
      break;
    default:
      break;
    }

    ss << seqs[(*c)->seqidIndexL].RefName 
       << "\t"
       << ((*c)->five - 10)
       << "\t"
       << ((*c)->five + 10)
       << "\t"
       << seqs[(*c)->seqidIndexL].RefName
       << "\t"
       << ((*c)->three - 10)
       << "\t"
       << ((*c)->three + 10)
       << "\t"
       << type << ":" << (*c)->id
       << "\t"
       << "."
       << "\t"
       << "."
       << "\t"
       << ".";

    ss << "\tSVTYPE=" << type << ";SVLEN=" << svlen << ";ID=" << (*c)->id << ";"
       << "SUPPORT=" << (*c)->supports[0] << "," << (*c)->supports[1] << ";" 
       << "MERGED=" << (*c)->merged << ";"
       << "REFINED=" << (*c)->refined << ";"
       << "POS=" << (*c)->five << "," << (*c)->three << ";" ;

    string SML = ".";
    string SMR = ".";

    if((*c)->sml.size() > 1){
      SML = join((*c)->sml, ",");
    }
    else{
      SML = (*c)->sml.front();
    }

    if((*c)->smr.size() > 1){
      SMR = join((*c)->smr, ",");
    }
    else{
      SMR = (*c)->smr.front();
    }

    ss << "LID=" << SML << ";" ; 
    ss << "RID=" << SMR ; 


   
    for(unsigned int i = 0; i < (*c)->genotypeIndex.size(); i++){
      if((*c)->genotypeIndex[i] == -1){
	ss << "\t" << "./.:" << "."
           << "," << "."
           << "," << "."
           << ":" << (*c)->nalt[i]
           << ":" << (*c)->nref[i];
      }
      else if((*c)->genotypeIndex[i] == 0){
	ss << "\t" << "0/0:" << (*c)->genotypeLikelhoods[i][0] 
	   << "," << (*c)->genotypeLikelhoods[i][1] 
	   << "," << (*c)->genotypeLikelhoods[i][2] 
	   << ":" << (*c)->nalt[i]
	   << ":" << (*c)->nref[i];
      }
      else if((*c)->genotypeIndex[i] == 1){
	ss << "\t" << "0/1:" << (*c)->genotypeLikelhoods[i][0] 
	   << "," << (*c)->genotypeLikelhoods[i][1] 
	   << "," << (*c)->genotypeLikelhoods[i][2] 
	   << ":" << (*c)->nalt[i]
	   << ":" << (*c)->nref[i];
      }
      else if((*c)->genotypeIndex[i] == 2){
	ss << "\t" << "1/1:" << (*c)->genotypeLikelhoods[i][0] 
	   << "," << (*c)->genotypeLikelhoods[i][1] 
	   << "," << (*c)->genotypeLikelhoods[i][2] 
	   << ":" << (*c)->nalt[i]
           << ":" << (*c)->nref[i];
      }     
      else{
	cerr << "FATAL: printBEDPE: unknown genotype." << endl;
	exit(1);
      }
    }
    ss << endl;
    cout << ss.str();
  }
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of ints

 Function does   : calculates the mean

 Function returns: double

*/

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

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of doubles

 Function does   : calculates the var

 Function returns: double

*/

double var(vector<double> & data, double mu){
  double variance = 0;

  for(vector<double>::iterator it = data.begin(); it != data.end(); it++){
    variance += pow((*it) - mu,2);
  }

  return variance / (data.size() - 1);
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : int pointer, vector cigDat

 Function does   : finds end position

 Function returns: void

*/
void endPos(vector<cigDat> & cigs, int * pos){

  for(vector<cigDat>::iterator it = cigs.begin();
      it != cigs.end(); it++){

    switch( (*it).Type ){
    case 'M':
      {
        *pos += (*it).Length;
        break;
      }
    case 'X':
      {
        *pos += (*it).Length;
        break;
      }
    case 'D':
      {
        *pos += (*it).Length;
        break;
      }
    case '=':
      {
        *pos += (*it).Length;
        break;
      }
    case 'N':
      {
	*pos += (*it).Length;
        break;
      }
    default:
      break;
    }
  }
  // WARNING: this needs to be double checked
   *pos -= 1;
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : node pointer, vector of node pointers

 Function does   : populates a vector for a sub graph

 Function returns: void

*/


void getTree(node * n, vector<node *> & ns){
  
  map<int, int>  seen ;
  vector<edge *> edges;
  
  if(!n->eds.empty()){
    edges.insert(edges.end(), n->eds.begin(), n->eds.end());
  }
  
  seen[n->pos] = 1;

  ns.push_back(n);

  // if something is pushed to the back of the vector it changes the positions ! be warned.

  while(!edges.empty()){
    
     uint hit = 0;
     
     if(seen.find(edges.back()->L->pos) != seen.end() && seen.find(edges.back()->R->pos) != seen.end() ){
       hit = 1;
     }
     else if(seen.find(edges.back()->L->pos) == seen.end() && seen.find(edges.back()->R->pos) == seen.end()){
       seen[edges.back()->L->pos] = 1;
       seen[edges.back()->R->pos] = 1;
       ns.push_back(edges.back()->L);
       ns.push_back(edges.back()->R);
       edges.insert(edges.end(), edges.back()->L->eds.begin(), edges.back()->L->eds.end());       
       edges.insert(edges.end(), edges.back()->R->eds.begin(), edges.back()->R->eds.end());
     }
     else if(seen.find(edges.back()->L->pos) == seen.end()){
       
       seen[edges.back()->L->pos] = 1;

       ns.push_back(edges.back()->L);
       
       edges.insert(edges.end(), edges.back()->L->eds.begin(), edges.back()->L->eds.end());
       
     }
     else{
       seen[edges.back()->R->pos] = 1;
       ns.push_back(edges.back()->R);

       edges.insert(edges.end(), edges.back()->R->eds.begin(), edges.back()->R->eds.end());
     }

    if(hit == 1){
      edges.pop_back();
    }
  }
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : breakpoints *, RefSeq

 Function does   : provides ref and alt

 Function returns: bool

*/

bool genAlleles(breakpoints * bp, string & fasta, RefVector & rv){
    
  FastaReference rs;

  omp_set_lock(&lock);
  rs.open(fasta);
  omp_unset_lock(&lock);

  string ref ;
  string alt ;


  bp->seqid = rv[bp->seqidIndexL].RefName;

  if(bp->type == 'D'){   

    if((bp->five - 500) < 0){
      bp->fail = true;
      return false;
    }
    if((bp->three + 500) > rv[bp->seqidIndexL].RefLength){
      bp->fail = true;
      return false;
    }

    ref = rs.getSubSequence(rv[bp->seqidIndexL].RefName, bp->five - 200, bp->svlen + 400 );
    alt = rs.getSubSequence(rv[bp->seqidIndexL].RefName, bp->five - 200, 200) + 
      // bp->five = first base of deletion -1 last ref base + 1 for fasta 
      rs.getSubSequence(rv[bp->seqidIndexL].RefName, bp->three+1, 200); // start one after deletion ends

    #ifdef DEBUG
    cerr << rs.getSubSequence(rv[bp->seqidIndexL].RefName, bp->five - 5, 5) 
	 << " -- " 
	 << rs.getSubSequence(rv[bp->seqidIndexL].RefName, bp->three+1,  5)
         << endl;
    #endif


  }
    //duplication;
  if(bp->type == 'U'){
    if((bp->five - 500) < 0){
      bp->fail = true;
      return false;
    }
    if((bp->three + 500) > rv[bp->seqidIndexL].RefLength){
      bp->fail = true;
      return false;
    }
    ref = rs.getSubSequence(rv[bp->seqidIndexL].RefName, bp->five -200 , bp->svlen + 400);
    alt = rs.getSubSequence(rv[bp->seqidIndexL].RefName, bp->five - 200 , bp->svlen + 200)
      + rs.getSubSequence(rv[bp->seqidIndexL].RefName, bp->five , bp->svlen)
      + rs.getSubSequence(rv[bp->seqidIndexL].RefName, bp->three , 200);

  }
  if(bp->type == 'V'){

    if((bp->five - 500) < 0){
      bp->fail = true;
      return false;
    }
    if((bp->three + 500) > rv[bp->seqidIndexL].RefLength){
      bp->fail = true;
      return false;
    }
    string inv = rs.getSubSequence(rv[bp->seqidIndexL].RefName, bp->five, (bp->svlen) );    
    inv = string(inv.rbegin(), inv.rend());
    Comp(inv);
    
    ref = rs.getSubSequence(rv[bp->seqidIndexL].RefName, bp->five -200, (bp->svlen + 400)) ;
    alt = rs.getSubSequence(rv[bp->seqidIndexL].RefName, bp->five-200, 200) + inv + rs.getSubSequence(rv[bp->seqidIndexL].RefName, bp->three, 200);

    
  }
    
  if(ref.size() > 800 && bp->type != 'U'){
    ref = ref.substr(0,400) + ref.substr(ref.size()-400,  400);
  }
  if(alt.size() > 800 && bp->type != 'U'){
    alt = alt.substr(0,400) + alt.substr(alt.size() -400, 400);
  }
  
  if(ref.size() > 1200 && bp->type == 'U'){
    ref = ref.substr(0,400) + ref.substr(ref.size()-400,  400);
    
  }
  if(alt.size() > 1200 && bp->type == 'U'){
  
    alt = alt.substr(0,400) + alt.substr(bp->svlen, 400) +  alt.substr(alt.size() -400, 400); 

  }
  

  bp->alleles.push_back(ref) ;
  bp->alleles.push_back(alt) ;

  //  cerr << "RS: " << ref.size() << " AS: " << alt.size() << endl;

  return true;
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : edge pointer 

 Function does   : init

 Function returns: void

*/

void initEdge(edge * e){
  e->L = NULL;
  e->R = NULL;


  e->support['L'] = 0;
  e->support['H'] = 0;
  e->support['S'] = 0;
  e->support['I'] = 0;
  e->support['D'] = 0;
  e->support['V'] = 0;
  e->support['M'] = 0;
  e->support['R'] = 0;
  e->support['X'] = 0;
  e->support['K'] = 0;

}




//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of ints and separator

 Function does   : joins vector with separator

 Function returns: string

*/

string join(vector<int> & ints, string sep){

  stringstream ss;

  for(vector<int>::iterator sit = ints.begin(); sit != ints.end(); sit++){
    ss << *sit << sep;
  }
  return ss.str();
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of strings

 Function does   : joins vector with returns;

 Function returns: string

*/

string joinReturn(vector<string> strings){

  string joined = "";

  for(vector<string>::iterator sit = strings.begin(); sit != strings.end(); 
      sit++){
    joined = joined + " " + (*sit) + "\n";
  }
  return joined;
}



//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : read pair pointer

 Function does   : -->s s<-- tests if pairs suggests insertion

 Function returns: bool

*/

bool isPointIn(readPair * rp){

  if(!rp->al1.IsMapped() || !rp->al2.IsMapped()){
    return false;
  }
  if(rp->al1.RefID != rp->al2.RefID){
    return false;
  }
  if(rp->al1.Position <= rp->al2.Position){
    
    if(rp->al1.CigarData.back().Type == 'S' && 
       rp->al2.CigarData.front().Type == 'S'){
      return true;
    }
  }
  else{
    if(rp->al1.CigarData.front().Type == 'S' && 
       rp->al2.CigarData.back().Type == 'S'){
      return true;
    }
  }
  return false;
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : refID (int), left position (int), graph

 Function does   : determine if a node is in the graph

 Function returns: bool

*/


bool isInGraph(int refID, int pos, graph & lc){

  if(lc.nodes.find(refID) == lc.nodes.end()){
    return false;
  }
 
  if(lc.nodes[refID].find(pos) != lc.nodes[refID].end()){
    return true;
  }
  return false;
}



//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : left position, right position, local graph

 Function does   : adds nodes or builds connections

 Function returns: void

*/

void addIndelToGraph(int refID, int l, int r, char s, string & SM){

  omp_set_lock(&glock);

  if( ! isInGraph(refID, l, globalGraph) 
      &&  ! isInGraph(refID, r, globalGraph) ){
    
    //cerr << "addIndelToGraph: neither node found" << endl;

    node * nodeL;
    node * nodeR;
    edge * ed   ;

    nodeL = new node;
    nodeR = new node;
    ed    = new edge;

    if(nodeL->sm.find("SM") == nodeL->sm.end()){
      nodeL->sm[SM] = 0;
    }
    else{
      nodeL->sm[SM] += 1;
    }

    if(nodeR->sm.find("SM") == nodeR->sm.end()){
      nodeR->sm[SM] = 0;
    }
    else{
      nodeR->sm[SM] += 1;
    }


    nodeL->collapsed = false;
    nodeR->collapsed = false;
    nodeL->beginSupport = 0;
    nodeL->endSupport   = 0;

    nodeR->beginSupport = 0;
    nodeR->endSupport   = 0;

    initEdge(ed);

    ed->support[s] +=1;

    ed->L = nodeL;
    ed->R = nodeR;

    nodeL->eds.push_back(ed);
    nodeR->eds.push_back(ed);

    nodeL->pos = l;
    nodeL->seqid = refID;
    nodeR->pos = r;
    nodeR->seqid = refID;

    globalGraph.nodes[refID][l] = nodeL;
    globalGraph.nodes[refID][r] = nodeR;

  }
 else if(isInGraph(refID, l, globalGraph) 
	 &&  ! isInGraph(refID, r, globalGraph)){

   node * nodeR;
   edge * ed;

   nodeR = new node;
   ed    = new edge;

   if(nodeR->sm.find("SM") == nodeR->sm.end()){
     nodeR->sm[SM] = 0;
   }
   else{
     nodeR->sm[SM] += 1;
   }

   if(globalGraph.nodes[refID][l]->sm.find("SM") == globalGraph.nodes[refID][l]->sm.end()){
     globalGraph.nodes[refID][l]->sm[SM] = 0;
   }
   else{
     globalGraph.nodes[refID][l]->sm[SM] += 1;
   }


   nodeR->collapsed = false;
   
   nodeR->beginSupport = 0;
   nodeR->endSupport = 0;

   initEdge(ed);
   ed->support[s] += 1;

   nodeR->pos   = r;
   nodeR->seqid = refID; 
   ed->L = globalGraph.nodes[refID][l];
   ed->R = nodeR;
   
   nodeR->eds.push_back(ed);

   globalGraph.nodes[refID][l]->eds.push_back(ed);
   globalGraph.nodes[refID][r] = nodeR;

 }
 else if(! isInGraph(refID, l, globalGraph) 
	 &&  isInGraph(refID, r, globalGraph)){

   //cerr << "addIndelToGraph: right node found" << endl;
   
   node * nodeL;
   edge * ed;

   nodeL = new node;
   ed    = new edge;

   if(nodeL->sm.find("SM") == nodeL->sm.end()){
     nodeL->sm[SM] = 0;
   }
   else{
     nodeL->sm[SM] += 1;
   }

   if(globalGraph.nodes[refID][r]->sm.find("SM") == globalGraph.nodes[refID][r]->sm.end()){
     globalGraph.nodes[refID][r]->sm[SM] = 0;
   }
   else{
     globalGraph.nodes[refID][r]->sm[SM] += 1;
   }


   nodeL->collapsed = false;
   nodeL->beginSupport = 0;
   nodeL->endSupport   = 0;


   initEdge(ed);
   ed->support[s] +=1;
   nodeL->pos   = l;
   nodeL->seqid = refID;
   ed->R = globalGraph.nodes[refID][r];
   ed->L = nodeL;

   //   cerr << "LPD RP: " << ed->R->pos << endl;
   
   nodeL->eds.push_back(ed);

   globalGraph.nodes[refID][r]->eds.push_back(ed);
   globalGraph.nodes[refID][l] = nodeL;  
}
 else{
   uint hit = 0;

   if(globalGraph.nodes[refID][l]->sm.find(SM) != globalGraph.nodes[refID][l]->sm.end()){
     globalGraph.nodes[refID][l]->sm[SM] = 0;
   }
   else{
     globalGraph.nodes[refID][l]->sm[SM] += 1;
   }
   if(globalGraph.nodes[refID][r]->sm.find(SM) != globalGraph.nodes[refID][r]->sm.end()){
     globalGraph.nodes[refID][r]->sm[SM] = 0;
   }
   else{
     globalGraph.nodes[refID][r]->sm[SM] += 1;
   }



   for(vector<edge *>::iterator ite 
	 = globalGraph.nodes[refID][l]->eds.begin();
       ite != globalGraph.nodes[refID][l]->eds.end(); ite++){
     if((*ite)->L->pos == l && (*ite)->R->pos == r){
       
       (*ite)->support[s] += 1;

       hit = 1;
     }
   }
   if(hit == 0){
     edge * ne;
     ne = new edge;
     initEdge(ne);
     ne->support[s]+=1;
     ne->L =      globalGraph.nodes[refID][l];
     ne->R =      globalGraph.nodes[refID][r];
     globalGraph.nodes[refID][l]->eds.push_back(ne);
     globalGraph.nodes[refID][r]->eds.push_back(ne);
   }
 }
  omp_unset_lock(&glock);
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : bam alignment, two ints.

 Function does   : finds the positions of indels

 Function returns: string

*/

bool indelToGraph(BamAlignment & ba, string & SM){

  if(!ba.IsMapped()){
    return false;
  }

  if(ba.MapQuality < 30){
    return false;
  }

  bool hit = false;

  int p = ba.Position;

  for(vector<CigarOp>::iterator ci = ba.CigarData.begin();
      ci != ba.CigarData.end(); ci++){
    
    switch(ci->Type){
    case 'M':
      {
	p += ci->Length;
	break;
      }
    case '=':
      {
        p += ci->Length;
        break;
      }
    case 'N':
      {
        p += ci->Length;
        break;
      }
    case 'X':
      {
	p += ci->Length;
	break;
      }
    case 'I':
      {
	hit = true;
	addIndelToGraph(ba.RefID, p, p + ci->Length, 'I', SM);
	break;
      }
    case 'D':
      {
	hit = true;
	addIndelToGraph(ba.RefID, p , (p + ci->Length ), 'D', SM);
	p += ci->Length;
	break;
      }
    default :
      {
	break;
      }
    }
  }
  return hit;
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of BamTools::CigarOp

 Function does   : joins vector " ";

 Function returns: string

*/

string joinCig(vector<CigarOp> strings){

  stringstream joined ;

  for(vector<CigarOp>::iterator sit = strings.begin(); 
      sit != strings.end(); sit++){
    joined  << (*sit).Length << (*sit).Type; 
  }
  return joined.str();
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of BamTools::CigarOp, int length

 Function does   : looks for clips greater than or equal to length

 Function returns: bool

*/

bool areBothClipped(vector<CigarOp> & ci){

  if(ci.front().Type == 'S' && ci.back().Type == 'S'){
    return true;
  }
  return false;
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of BamTools::CigarOp, int length

 Function does   : looks for clips greater than or equal to length

 Function returns: bool

*/

bool IsLongClip(vector<CigarOp> & ci, unsigned int len){

  if(ci.front().Type == 'S' && ci.front().Length >= len){
    return true;
  }
  if(ci.back().Type == 'S' && ci.back().Length >= len){
    return true;
  }
  return false;
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector<CigarOp>

 Function does   : calculates the number of matching bases 

 Function returns: unsigned int, total number of matching bases

*/

int match(vector<CigarOp> co){
  
  int m = 0;
  
  for(vector<CigarOp>::iterator it = co.begin(); 
      it != co.end(); it++){
    if(it->Type == 'M'){
      m += it->Length;
    }
  }
  
  return m;
  
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : read pair to filter

 Function does   : fails poor quality or non-sv read pairs

 Function returns: true = bad, false = good;

*/

bool pairFailed(readPair * rp){
  
  if(rp->al1.IsMapped() && rp->al2.IsMapped()){
    if(rp->al1.Length == rp->al1.CigarData[0].Length 
       && rp->al1.CigarData[0].Type == 'M' &&
       rp->al2.Length == rp->al2.CigarData[0].Length 
       && rp->al2.CigarData[0].Type == 'M' ){
      return true;
    }
    if(rp->al1.MapQuality < globalOpts.MQ && rp->al2.MapQuality < globalOpts.MQ){
      return true;
    }
    if((match(rp->al1.CigarData) + match(rp->al2.CigarData)) < 75){
      return true;
    }
  }
  return false;
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector<cigDat>, string

 Function does   : loads cigar data into a string

 Function returns: NA

*/

void parseCigar(vector<cigDat> & parsedCigar, string cigar){

  unsigned int spot = 0;

  for(unsigned int i = 0; i < cigar.size(); i++){
    if(int(cigar[i] > 57)){
      cigDat tup;
      tup.Length = atoi(cigar.substr(spot, i-spot).c_str());
      tup.Type  = cigar[i];
      parsedCigar.push_back(tup);
    }
  }
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : SA tag string, and vector of saTag

 Function does   : parses string and converts positions to BAM

 Function returns: NA

*/

void parseSA(vector<saTag> & parsed, string tag, map<string, int> & il){

  vector<string> sas = split(tag, ';');

  for(unsigned int i = 0 ; i < sas.size() -1 ; i++){
    
    saTag sDat;
    
    vector<string> sat = split (sas[i], ',');

    if(sat.size() != 6){
      cerr << "FATAL: failure to parse SA optional tag" << endl;
      exit(1);
    }

    if(sat[2].compare("-") == 0){
      sDat.strand = true;
    }
    else{
      sDat.strand = false;
    }

    sDat.seqid = il[sat[0]];
    sDat.pos   = atoi(sat[1].c_str()) - 1;
    parseCigar(sDat.cig, sat[3]);
    parsed.push_back(sDat);

  }
  
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : bam alignment, vector<saTag> 

 Function does   : finds links for splitter to put in graph

 Function returns: NA

*/

void splitToGraph(BamAlignment & al, vector<saTag> & sa, string & SM){

  if(!al.IsMapped()){
    return;
  }

  if(sa.size() > 1){
    return;
  }

  if(sa[0].seqid != al.RefID){
    return;
  }

  char support = 'S';

  if((sa[0].strand && ! al.IsReverseStrand()) || (! sa[0].strand && al.IsReverseStrand() )){
    support = 'V';
  }

  if(sa.front().cig.front().Type == 'S' 
     && sa.front().cig.back().Type == 'S'){
    return;
  }

  if(al.CigarData.front().Type == 'S' 
     && al.CigarData.back().Type == 'S'){
    return;
  }

  if(al.CigarData.front().Type == 'S'){
    
    int start = al.Position; 
    int end   = sa.front().pos  ;

    if(sa.front().cig.back().Type == 'S'){
      endPos(sa[0].cig, &end) ;
    
      if(end > start ){
	support = 'X';
      }
    }
    
    if(start > end){
      int tmp = start;
      start = end;
      end   = tmp;
    }
    addIndelToGraph(al.RefID, start, end, support, SM);
  }
  else{
    int start = al.GetEndPosition(false,true);
    int end   = sa.front().pos                ;
    if(sa[0].cig.back().Type == 'S'){
      endPos(sa.front().cig, &end);      
    }
    else{
      if(start > end){
	support = 'X';
      }
    }
    if(start > end){      
      start = sa[0].pos;
      end   = al.GetEndPosition(false,true);
    }
    addIndelToGraph(al.RefID, start, end, support, SM);
  }
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : pointer to readPair;

 Function does   : adds high and low insert sizes to graph. both pairs are 
                   treated as mapped

 Function returns: NA

*/

void deviantInsertSize(readPair * rp, char supportType, string & SM){

  if(! IsLongClip(rp->al1.CigarData, 1) 
     && ! IsLongClip(rp->al2.CigarData, 1)){
    return;
  }
  if(rp->al1.CigarData.front().Type == 'S' 
     || rp->al1.CigarData.back().Type == 'S'){
    int start = rp->al1.Position;
    int end   = rp->al2.Position;
    if(rp->al2.CigarData.back().Type == 'S'){
      end = rp->al2.GetEndPosition();
    }

    if(rp->al1.CigarData.back().Type == 'S'){
      start = rp->al1.GetEndPosition(false,true);
    }
    if(start > end){
      int tmp = end;
      end = start  ;
      start = tmp  ;
    }
    addIndelToGraph(rp->al1.RefID, start, end, supportType, SM);       
  }
  else{
    int start = rp->al2.Position;
    int end   = rp->al1.Position;

    if(rp->al1.CigarData.back().Type == 'S'){
      end = rp->al1.GetEndPosition();
    }

    if(rp->al2.CigarData.back().Type == 'S'){
      start = rp->al2.GetEndPosition(false,true);
    }    
    if(start > end){
      int tmp = end;
      end = start  ;
      start = tmp  ;
    }
    addIndelToGraph(rp->al2.RefID, start, end, supportType, SM);
    
  }

}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : pointer to readPair ; seqid to index

 Function does   : processes Pair

 Function returns: NA

*/

void processPair(readPair * rp,  map<string, int> & il, 
		 double   * low, double * high, string & SM){

#ifdef DEBUG
  cerr << "processing pair" << endl;
#endif 

  
  string sa1;
  string sa2;
  
  bool sameStrand = false;

  if(pairFailed(rp)){

#ifdef DEBUG
    cerr << "pair failed filter: " << rp->al1.Name  << endl;
    cerr << "                    " << rp->al1.Position << endl;
    cerr << "                    " << joinCig(rp->al1.CigarData) << endl;
    cerr << "                    " << rp->al2.Position << endl;
    cerr << "                    " << joinCig(rp->al2.CigarData) << endl;

#endif 
    return;
  }
 
  // not doing translocations 
 
  if(rp->al1.RefID != rp->al2.RefID){
    return;
  }

  indelToGraph(rp->al1, SM);
  indelToGraph(rp->al2, SM);
  
  if( rp->al1.IsMapped() && rp->al2.IsMapped() ){
    if( ! IsLongClip(rp->al1.CigarData, 5) 
	&& ! IsLongClip(rp->al2.CigarData, 5)){
      return;
    }
    if((rp->al1.IsReverseStrand() && rp->al2.IsReverseStrand())
       || (! rp->al1.IsReverseStrand() && ! rp->al2.IsReverseStrand()) ){
      sameStrand = true;
    }
    if( abs(rp->al1.InsertSize) > *high){
      //      cerr << rp->al1.Name << " H " << rp->al1.InsertSize << endl;
      if(sameStrand){
	deviantInsertSize(rp, 'M', SM); 
      }
      else{
	deviantInsertSize(rp, 'H', SM); 
      }
    }    
    if( abs(rp->al1.InsertSize) < *low ){
      //      cerr << rp->al1.Name << " L " << rp->al1.InsertSize << endl;
      if(sameStrand){
	deviantInsertSize(rp, 'R', SM);
      }
      else{
	deviantInsertSize(rp, 'L', SM);
      }
    }
  }
  // one is not mapped
  else{
    //   mateNotMapped(rp, 'K');
  }
          
  if(rp->al1.GetTag("SA", sa1)){
    vector<saTag> parsedSa1;
    parseSA(parsedSa1, sa1, il);
    //    cerr << sa1 << endl;
    splitToGraph(rp->al1, parsedSa1, SM);
    //    cerr << "s1 processed " << endl;
  }
  if(rp->al2.GetTag("SA", sa2)){
    vector<saTag> parsedSa2;
    parseSA(parsedSa2, sa2, il);
    //    cerr << sa2 << endl;
    splitToGraph(rp->al2, parsedSa2, SM);
    //    cerr << "s2 processed " << endl;
  }
  
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : filename, seqid, start, end, refSeqs (bamtools object),
                   paired end store. 

 Function does   : Matches paired ends so they can be processed together.
                   The global pair store is loaded for the missing mates

 Function returns: NA

*/

bool runRegion(string filename,
	       int seqidIndex,
               int start,
               int end,
	       vector< RefData > seqNames,
	       map<string, readPair *> & globalPairStore,
	       map<string, int> & seqInverseLookup, string & SM){



#ifdef DEBUG
  cerr << "running region: " << seqidIndex << ":" << start << "-" << end << endl;
#endif 

  map<string, int> localInverseLookUp;

  omp_set_lock(&lock);
  for(map<string, int>::iterator itt = seqInverseLookup.begin();
      itt !=  seqInverseLookup.end(); itt++){    
    localInverseLookUp[itt->first] = itt->second;
  }
  insertDat localDists  = insertDists;
  omp_unset_lock(&lock);

  // local graph;
  graph localGraph;

  // local read pair store
  
  double high = localDists.mus[filename] + (2.5 * localDists.sds[filename]);
  double low  = localDists.mus[filename] - (2.5 * localDists.sds[filename]);

  if(low < 0){
    low = 100;
  }

  map<string, readPair *>pairStore;

  BamReader br;
  br.Open(filename);

  if(! br.LocateIndex()){
    vector<string> fileName = split(filename, ".");
    fileName.back() = "bai";
    string indexName = join(fileName, ".");
    if(! br.OpenIndex(indexName) ){
      cerr << "FATAL: cannot find bam index." << endl;
    }
  }
  if(!br.SetRegion(seqidIndex, start, seqidIndex, end)){
    return false;
  }

  BamAlignment al;

  while(br.GetNextAlignmentCore(al)){
    if((al.AlignmentFlag & 0x0800) != 0 ){
      continue;
    }
    if(! al.IsPaired()){
      continue;
    }
    
    if(! al.IsMapped() && ! al.IsMateMapped()){
      continue;
    }
    if(al.IsDuplicate()){
      continue;
    }
    if(! al.IsPrimaryAlignment()){
      continue;
    }

    // dna and read name
    al.BuildCharData();

    if(pairStore.find(al.Name) != pairStore.end()){
      pairStore[al.Name]->count += 1;
      if(pairStore[al.Name]->flag == 2){
	pairStore[al.Name]->al1 = al;
      }
      else{
	pairStore[al.Name]->al2 = al;
      }
      processPair(pairStore[al.Name], localInverseLookUp, &low, &high, SM);
      delete pairStore[al.Name];
      pairStore.erase(al.Name);	 
    }
    else{
      readPair * rp;
      rp = new readPair;
      rp->count = 1;
      if(al.IsFirstMate()){
	rp->flag = 1 ;
	rp->al1  = al;
      }
      else{
	rp->flag = 2 ;
	rp->al2  = al;
      } 
      pairStore[al.Name] = rp;
    }
  }

  // close the bam
  br.Close();

  // load lonely reads into the global struct;
  // if it finds a mate in the global store it processes and deletes
  omp_set_lock(&lock);

  for(map<string, readPair *>::iterator rps = pairStore.begin();
      rps != pairStore.end(); rps++){

    if(globalPairStore.find(rps->first) != globalPairStore.end()){
      globalPairStore[rps->first]->count += 1;
      
      if(globalPairStore[rps->first]->flag == 1 && (*rps->second).flag == 1){
	continue;
      }
      if(globalPairStore[rps->first]->flag == 2 && (*rps->second).flag == 2){
	continue;
      }

      if(globalPairStore[rps->first]->flag == 1){
//	cerr << "c: "     << globalPairStore[rps->first]->count << endl;
//	cerr << "f: "     << globalPairStore[rps->first]->flag << endl;
//	cerr << "test1A " << (*rps->second).al2.Name << endl;
//	cerr << "test1B " << (*rps->second).al1.Name << endl;
	(*globalPairStore[rps->first]).al2 = (*rps->second).al2;
      }
      else{
//	cerr << "c: "     << globalPairStore[rps->first]->count << endl;
//	cerr << "f: "     << globalPairStore[rps->first]->flag << endl;
//	cerr << "test2A " << (*rps->second).al1.Name << endl;
//	cerr << "test2B " << (*rps->second).al2.Name << endl;
      	(*globalPairStore[rps->first]).al1 = (*rps->second).al1;
      }
//      cerr << "INFO: about to process read pairs in different regions: " << globalPairStore[rps->first]->flag << " " 
//	   << globalPairStore[rps->first]->al1.Name      << " "
//	   << globalPairStore[rps->first]->al1.AlignmentFlag      << " "
//	   << globalPairStore[rps->first]->al1.RefID     << " " 
//	   << globalPairStore[rps->first]->al1.Position  << " " 
//	   << globalPairStore[rps->first]->al2.Name      << " " 
//	   << globalPairStore[rps->first]->al2.AlignmentFlag      << " "
//	   << globalPairStore[rps->first]->al2.RefID     << " " 
//	   << globalPairStore[rps->first]->al2.Position  << endl;
      processPair(globalPairStore[rps->first], localInverseLookUp, &low, &high, SM);
      delete globalPairStore[rps->first];
      delete pairStore[rps->first];
      globalPairStore.erase(rps->first);
    }
    else{
      globalPairStore[rps->first] = pairStore[rps->first];
    }
  }

  omp_unset_lock(&lock);
  return true;
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : string

 Function does   : reads bam into graph

 Function returns: NA

*/

void loadBam(string & bamFile){

  bool region    = false;
  int start      = 0;
  int end        = 0;
  string regionSID;

  omp_set_lock(&lock);
  
  if(!globalOpts.seqid.empty()){
    region = true;
    regionSID = globalOpts.seqid;
    start = globalOpts.region.front();
    end   = globalOpts.region.back();
  }

  omp_unset_lock(&lock);

  cerr << "INFO: reading bam file: " << bamFile << endl;

  BamReader br;
 
  if(! br.Open(bamFile)){
    cerr << "\n" << "FATAL: could not open: " << bamFile << endl;
    exit(1);
  }
  if(! br.LocateIndex()){    
    vector<string> fileName = split(bamFile, ".");
    fileName.back() = "bai";
    string indexName = join(fileName, ".");
    
    cerr << "INFO: Did not find *bam.bai" << endl;
    cerr << "INFO: Trying: " << indexName << endl;
    
    if(! br.OpenIndex(indexName) ){
      cerr << "FATAL: cannot find bam index." << endl;
    }
  }
  
  // grabbing the header
  SamHeader SH = br.GetHeader();
  
  if(!SH.HasReadGroups()){
    cerr << endl;
    cerr << "FATAL: No @RG detected in header.  WHAM uses \"SM:sample\"." << endl;
    cerr << endl;
  }
  
  SamReadGroupDictionary RG = SH.ReadGroups;
  
  if(RG.Size() > 1){
    cerr << endl;
    cerr << "WARNING: Multiple libraries (@RG). Assuming same library prep." << endl;
    cerr << "WARNING: Multiple libraries (@RG). Assuming same sample (SM)." << endl;
    cerr << endl;
  }

  string SM;
  
  if(!RG.Begin()->HasSample()){
    cerr << endl;
    cerr << "FATAL: No SM tag in bam file." << endl;
    exit(1);
    cerr << endl;
  }
  
  SM = RG.Begin()->Sample;
  
  // if the bam is not sorted die
  if(!SH.HasSortOrder()){
    cerr << "FATAL: sorted bams must have the @HD SO: tag in each SAM header: " << bamFile  << endl;
    exit(1);
  }

  RefVector sequences = br.GetReferenceData();

  //chunking up the genome
  
  vector< regionDat* > regions;

  // inverse lookup for split reads 
  map<string,int> seqIndexLookup;
  
  int seqidIndex = 0;
  for(vector< RefData >::iterator sit = sequences.begin(); sit != sequences.end(); sit++){
      seqIndexLookup[sequences[seqidIndex].RefName] = seqidIndex;
      seqidIndex += 1;
  }
  if(region){
    
    int p = start;
    int e = 0;
    for(; (p+1000000) <= end; p += 1000000){
      regionDat * regionInfo = new regionDat;
      regionInfo->seqidIndex = seqIndexLookup[regionSID];
      regionInfo->start      = p                      ;
      regionInfo->end        = 1000000 + p            ; 
      regions.push_back(regionInfo);   
      e = p + 1000000; 
    }
    if(e < end){
      regionDat * regionInfo = new regionDat;
      regionInfo->seqidIndex = seqIndexLookup[regionSID];
      regionInfo->start      = p                        ;
      regionInfo->end        = end                      ;
      regions.push_back(regionInfo);
    }

  }
  else{
    
    seqidIndex = 0;
    
    for(vector< RefData >::iterator sit = sequences.begin(); sit != sequences.end(); sit++){
      int start = 0;
      
      if(globalOpts.toSkip.find( (*sit).RefName ) == globalOpts.toSkip.end() && ((*sit).RefLength > 1000)){
	for(;start < (*sit).RefLength ; start += 1000000){
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
      else{
	cerr << "INFO: skipping: " << (*sit).RefName << endl;
	seqidIndex += 1;
      }
    }
  }
  // closing the bam reader before running regions
  br.Close();
  
  // global read pair store

  map<string, readPair*> pairStore;



  int Mb = 0;

  // running the regions with openMP
#pragma omp parallel for schedule(dynamic, 3)
 
  for(unsigned int re = 0; re < regions.size(); re++){
    if(! runRegion(bamFile, 
		   regions[re]->seqidIndex, 
		   regions[re]->start, 
		   regions[re]->end, 
		   sequences,
		   pairStore,
		   seqIndexLookup, SM)){
      omp_set_lock(&lock);
      cerr << "WARNING: region failed to run properly: "
           << sequences[regions[re]->seqidIndex].RefName
           << ":"  << regions[re]->start << "-"
           << regions[re]->end
           <<  endl;
      omp_unset_lock(&lock);
    }
    else{
      delete regions[re];
      omp_set_lock(&lock);
      Mb += 1;
      if((Mb % 10) == 0 ){
	cerr << "INFO: " << SM
	     << ": processed "
	     << Mb << "Mb of the genome." << endl;	
      }
      omp_unset_lock(&lock);
    }
  }
  cerr << "INFO: " << bamFile << " had " 
       << pairStore.size() 
       << " reads that were not processed" 
       << endl; 

  // freeing memory
  for(map<string, readPair *>::iterator rp = pairStore.begin();
      rp != pairStore.end(); rp++){
    delete rp->second;
  }

}
//-------------------------------   OPTIONS   --------------------------------
int parseOpts(int argc, char** argv)
    {
    int opt = 0;
    opt = getopt(argc, argv, optString);
    while(opt != -1){
      switch(opt){
      case 'k':
	{
	  globalOpts.skipGeno = true;
	  break;
	}
      case 'b':
	{
	  globalOpts.svs = optarg;
	  cerr << "INFO: WHAM-GRAPHENING will only genotype input: " << globalOpts.svs << endl;
	  break;
	}
      case 's':
	{
	  globalOpts.statsOnly = true;
	  break;
	}
      case 'g':
	{
	  globalOpts.graphOut = optarg;
	  cerr << "INFO: graphs will be written to: " <<  globalOpts.graphOut 
	       << endl;
	  break;
	}
      case 'a':
	{
	  globalOpts.fasta = optarg;
	  cerr << "INFO: fasta file: " << globalOpts.fasta << endl;
	  break;
	}
      case 'e':
	{
	  vector<string> seqidsToSkip = split(optarg, ",");
	  for(unsigned int i = 0; i < seqidsToSkip.size(); i++){
	    globalOpts.toSkip[seqidsToSkip[i]] = 1;
	    cerr << "INFO: WHAM will skip seqid: " << seqidsToSkip[i] << endl;
	  }
	  break;
	}
      case 'f':
	{
	  globalOpts.targetBams     = split(optarg, ",");
	  cerr << "INFO: target bams:\n" << joinReturn(globalOpts.targetBams) ;
	  break;
	}
      case 'h':
	{
	  printHelp();
	  exit(1);
	  break;
	}
      case '?':
	{
	  break;
	}
      case 'm':
	{
	  globalOpts.MQ = atoi(((string)optarg).c_str());
	  cerr << "INFO: Reads with mapping quality below " << globalOpts.MQ << " will be filtered. " << endl;
	  break;
	}
      case 'x':
	{
	  globalOpts.nthreads = atoi(((string)optarg).c_str());
	  cerr << "INFO: OpenMP will roughly use " << globalOpts.nthreads 
	       << " threads" << endl;
	  break;
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

      }
      opt = getopt( argc, argv, optString ); 
    }
    return 1;
    }

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector<nodes *>

 Function does   : dumps a graph in dot format

 Function returns: string

*/

string dotviz(vector<node *> & ns){

  stringstream ss;

  ss << "graph {\n";

  for(vector<node *>::iterator it = ns.begin(); 
      it != ns.end(); it++){
    for(vector<edge *>:: iterator iz = (*it)->eds.begin(); 
	iz != (*it)->eds.end(); iz++){



      if((*iz)->support['X'] > 0){
        if((*it)->pos != (*iz)->L->pos){
          ss << "     " << (*it)->seqid << "." << (*it)->pos << " -- " << (*iz)->L->seqid << "." << (*iz)->L->pos << " [style=dashed,penwidth=" << (*iz)->support['X'] << "];\n";
        }
        if((*it)->pos != (*iz)->R->pos){
          ss << "     " << (*it)->seqid << "." << (*it)->pos << " -- " << (*iz)->R->seqid << "." << (*iz)->R->pos << " [style=dashed,penwidth=" << (*iz)->support['X'] << "];\n";
        }
      }




      
      if((*iz)->support['R'] > 0){
	if((*it)->pos != (*iz)->L->pos){
	  ss << "     " << (*it)->seqid << "." << (*it)->pos << " -- " << (*iz)->L->seqid << "." << (*iz)->L->pos << " [color=yellow,penwidth=" << (*iz)->support['R'] << "];\n";
	}
	if((*it)->pos != (*iz)->R->pos){
	  ss << "     " << (*it)->seqid << "." << (*it)->pos << " -- " << (*iz)->R->seqid << "." << (*iz)->R->pos << " [color=yellow,penwidth=" << (*iz)->support['R'] << "];\n";
	}
      }      
      
      if((*iz)->support['M'] > 0){
        if((*it)->pos != (*iz)->L->pos){
          ss << "     " << (*it)->seqid << "." << (*it)->pos << " -- " << (*iz)->L->seqid << "." << (*iz)->L->pos << " [color=magenta,penwidth=" << (*iz)->support['M'] << "];\n";
        }
        if((*it)->pos != (*iz)->R->pos){
          ss << "     " << (*it)->seqid << "." << (*it)->pos << " -- " << (*iz)->R->seqid << "." << (*iz)->R->pos << " [color=magenta,penwidth=" << (*iz)->support['M'] << "];\n";
        }
      }



      if((*iz)->support['V'] > 0){
        if((*it)->pos != (*iz)->L->pos){
          ss << "     " << (*it)->seqid << "." << (*it)->pos << " -- " << (*iz)->L->seqid << "." << (*iz)->L->pos << " [colo\
r=green,penwidth=" << (*iz)->support['V'] << "];\n";
        }
        if((*it)->pos != (*iz)->R->pos){
          ss << "     " << (*it)->seqid << "." << (*it)->pos << " -- " << (*iz)->R->seqid << "." << (*iz)->R->pos << " [colo\
r=green,penwidth=" << (*iz)->support['V'] << "];\n";
        }
      }      

      if((*iz)->support['L'] > 0){
	if((*it)->pos != (*iz)->L->pos){
          ss << "     " << (*it)->seqid << "." << (*it)->pos << " -- " << (*iz)->L->seqid << "." << (*iz)->L->pos << " [color=brown,penwidth=" << (*iz)->support['L'] << "];\n";
	}
        if((*it)->pos != (*iz)->R->pos){
          ss << "     " << (*it)->seqid << "." << (*it)->pos << " -- " << (*iz)->R->seqid << "." << (*iz)->R->pos << " [color=brown,penwidth=" << (*iz)->support['L'] << "];\n";
	}
      }
      if((*iz)->support['H'] > 0){
        if((*it)->pos != (*iz)->L->pos){
          ss << "     " << (*it)->seqid << "." << (*it)->pos << " -- " << (*iz)->L->seqid << "." << (*iz)->L->pos << " [color=purple,penwidth=" << (*iz)->support['H'] << "];\n";
        }
        if((*it)->pos != (*iz)->R->pos){
          ss << "     " << (*it)->seqid << "." << (*it)->pos << " -- " << (*iz)->R->seqid << "." << (*iz)->R->pos << " [color=purple,penwidth=" << (*iz)->support['H'] << "];\n";
	}
      }

      if((*iz)->support['S'] > 0){
	if((*it)->pos != (*iz)->L->pos){
	  ss << "     " << (*it)->seqid << "." << (*it)->pos << " -- " << (*iz)->L->seqid << "." << (*iz)->L->pos << " [style=dotted,penwidth=" << (*iz)->support['S'] << "];\n";
	}
	if((*it)->pos != (*iz)->R->pos){
	  ss << "     " << (*it)->seqid << "." << (*it)->pos << " -- " << (*iz)->R->seqid << "." << (*iz)->R->pos << " [style=dotted,penwidth=" << (*iz)->support['S'] << "];\n";
	}
      }
      if((*iz)->support['I'] > 0){
	if((*it)->pos != (*iz)->L->pos){
	  ss << "     " << (*it)->seqid << "." << (*it)->pos <<  " -- " << (*iz)->L->seqid << "." << (*iz)->L->pos << " [color=red,penwidth=" << (*iz)->support['I'] << "];\n";
	}
	if((*it)->pos != (*iz)->R->pos){
	  ss << "     " << (*it)->seqid << "." << (*it)->pos << " -- " << (*iz)->R->seqid << "." << (*iz)->R->pos << " [color=red,penwidth=" << (*iz)->support['I'] << "];\n";
	}
      }
      if((*iz)->support['D'] > 0){
	if((*it)->pos != (*iz)->L->pos){
	  ss << "     " << (*it)->seqid << "." << (*it)->pos << " -- " << (*iz)->L->seqid << "." << (*iz)->L->pos << " [color=blue,penwidth=" << (*iz)->support['D'] << "];\n";
	}
	if((*it)->pos != (*iz)->R->pos){
	  ss << "     " << (*it)->seqid << "." << (*it)->pos << " -- " << (*iz)->R->seqid << "." << (*iz)->R->pos << " [color=blue,penwidth=" << (*iz)->support['D'] << "];\n";
	}
      }
    }
  }

  ss << "}";


  return ss.str();
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : purges poorly supported trees in forest

 Function does   : dumps and shrinks graph

 Function returns: NA

*/

    void thin(){
  
  map<int, map<int, int> > lookup;
  
  map<int, map<int, int> > toDelete; 
  
  
  for(map<int, map<int, node* > >::iterator it = globalGraph.nodes.begin();it != globalGraph.nodes.end(); it++){
    for(map<int, node*>::iterator itt = it->second.begin(); itt != it->second.end(); itt++){
      
      if(lookup[it->first].find(itt->first) != lookup[it->first].end() ){
      }
      else{
	lookup[it->first][itt->first] = 1;
	
	vector<node *> tree;
	
	getTree(globalGraph.nodes[it->first][itt->first], tree);
	
	int flag = 0;
	
	for(vector<node *>::iterator ir = tree.begin(); ir != tree.end(); ir++){
	  lookup[(*ir)->seqid][(*ir)->pos] = 1;
	  for(vector<edge *>::iterator iz = (*ir)->eds.begin(); 
	      iz != (*ir)->eds.end(); iz++){
	    if((*iz)->support['I'] > 2 || (*iz)->support['D'] > 2 || (*iz)->support['S'] > 2 || (*iz)->support['L'] > 2 || (*iz)->support['R'] > 2 ){
	      flag = 1;
	    }
	  }
	}
	
	if(flag == 0){
	  for(vector<node *>::iterator ir = tree.begin(); ir != tree.end(); ir++){
	    toDelete[(*ir)->seqid][(*ir)->pos] = 1;	    
	  }
	}
      }
    }
  }
  
  for(map< int, map<int, int> >::iterator td = toDelete.begin(); td != toDelete.end(); td++){
    for(map<int, int>::iterator tdz = toDelete[td->first].begin(); tdz != toDelete[td->first].end(); tdz++){
      
      for(vector<edge *>::iterator etd = globalGraph.nodes[td->first][tdz->first]->eds.begin(); 
	  etd != globalGraph.nodes[td->first][tdz->first]->eds.end(); etd++){
	//	delete (*etd);
      }

      delete globalGraph.nodes[td->first][tdz->first];
      globalGraph.nodes[td->first].erase(tdz->first);
    }
  }
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of node pointers

 Function does   : tries to resolve a deletion

 Function returns: NA

*/

bool detectInsertion(vector<node *> tree, breakpoints * bp){

  vector <node *> putative;
  
  for(vector<node * >::iterator t = tree.begin(); t != tree.end(); t++){
    
    if((*t)->eds.size() > 1 ){
      
      int tooClose    = 0;
      int splitR      = 0;
      int insertion   = 0;

      for(vector<edge *>::iterator es = (*t)->eds.begin(); es != (*t)->eds.end(); es++){
        tooClose  += (*es)->support['H'];
	splitR    += (*es)->support['S'];
	insertion += (*es)->support['I'];
      }
      if( tooClose > 0 && insertion  > 0 && splitR == 0){
	putative.push_back((*t));
      }
    } 
  }
  if(putative.size() == 2){

    sort(putative.begin(), putative.end(), sortNodesByPos);

    int lPos = putative.front()->pos;
    int rPos = putative.back()->pos ;

    int lhit = 0 ; int rhit = 0;

    for(vector<edge *>::iterator ed = putative.front()->eds.begin() ;
        ed != putative.front()->eds.end(); ed++){
      if(((*ed)->L->pos == rPos) || ((*ed)->R->pos == rPos)){
        lhit = 1;
        break;
      }
    }

    for(vector<edge *>::iterator ed = putative.back()->eds.begin() ;
        ed != putative.back()->eds.end(); ed++){
      if(((*ed)->L->pos == lPos) || ((*ed)->R->pos == lPos)){
	rhit = 1;
        break;
      }
    }
    
    if(lhit == 1 && rhit == 1){
      //      cerr << "insertion pair: " << putative.front()->seqid  << " " <<  lPos << "\t" << rPos << endl;
      return true;
    }
    else{
      cerr << "no linked putative breakpoints" << endl;
    }
    
  }

  return false;
}



//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : two node pointer

 Function does   : finds if nodes are directly connected

 Function returns: bool

*/

bool connectedNode(node * left, node * right){

  for(vector<edge *>::iterator l = left->eds.begin(); l != left->eds.end(); l++){
    
    if((*l)->L->pos == right->pos || (*l)->R->pos == right->pos){
      return true;
    }
  }
  return false;
}


/*
  Function input  : a vector of edge pointers 

  Function does   : does a tree terversal and joins nodes

  Function returns: NA

*/


bool findEdge(vector<edge *> & eds, edge ** e, int pos){
  
  //  cerr << "finding edge" << endl;

  for(vector<edge *>::iterator it = eds.begin(); it != eds.end(); it++){
    if((*it)->L->pos == pos || (*it)->R->pos == pos ){
      (*e) = (*it);
      // cerr << "found edge: " << (*e)->L->pos  << endl;
      return true;
    }
  }

  return false;
}



//------------------------------- SUBROUTINE --------------------------------
/*
  Function input  : a vector of node pointers, and a postion

  Function does   : removes edges that contain a position

  Function returns: NA

*/


void removeEdges(vector<node *> & tree, int pos){

  //  cerr << "removing edges" << endl;

  for(vector<node *>::iterator rm = tree.begin(); rm != tree.end(); rm++){

    vector<edge *> tmp;

    for(vector<edge *>::iterator e = (*rm)->eds.begin(); e != (*rm)->eds.end(); e++){
      if( (*e)->L->pos != pos && (*e)->R->pos != pos  ){
	tmp.push_back((*e));
      }
      else{
	// delete the edge

      }
    }
    (*rm)->eds.clear();
    (*rm)->eds.insert((*rm)->eds.end(), tmp.begin(), tmp.end());
  }
}


//------------------------------- SUBROUTINE --------------------------------
/*
  Function input  : a vector of node pointers

  Function does   : does a tree terversal and joins nodes

  Function returns: NA

*/

void joinNodes(node * L, node * R, vector<node *> & tree){


  // quantifying which node has more support

  int lSupport = 0;
  int rSupport = 0;
  
  for(vector<edge *>::iterator le = L->eds.begin(); le != L->eds.end(); le++){
    lSupport += (*le)->support['D'] + (*le)->support['S'] + (*le)->support['I']
      +  (*le)->support['H'] +  (*le)->support['L'] + (*le)->support['X'] 
      + (*le)->support['M'] + (*le)->support['R'] + (*le)->support['V'] ;
  }
  for(vector<edge *>::iterator re = R->eds.begin(); re != R->eds.end(); re++){
    lSupport += (*re)->support['D'] + (*re)->support['S'] + (*re)->support['I']
      +  (*re)->support['H'] +  (*re)->support['L'] +  (*re)->support['X'] 
      +  (*re)->support['M'] +  (*re)->support['R'] +  (*re)->support['V'] ;
  }
  
  if(lSupport <= rSupport){

    L->collapsed = true;



    for(vector<edge *>::iterator lc =  L->eds.begin(); lc != L->eds.end(); lc++){


      edge * e; 

      int otherP = (*lc)->L->pos;

      if(L->pos  == otherP){
	otherP = (*lc)->R->pos;
      }
      if(findEdge(R->eds, &e,  otherP)){

	e->support['I'] += (*lc)->support['I'];
	e->support['D'] += (*lc)->support['D'];
	e->support['S'] += (*lc)->support['S'];
	e->support['H'] += (*lc)->support['H'];
	e->support['L'] += (*lc)->support['L'];
	e->support['R'] += (*lc)->support['R'];
	e->support['M'] += (*lc)->support['M'];
	e->support['V'] += (*lc)->support['V'];
	e->support['X'] += (*lc)->support['X'];
      }
      else{
	if((*lc)->L->pos == otherP){
	  (*lc)->R = R;
	}
	else{
	  (*lc)->L = R;
	}

	R->eds.push_back(*lc);
      }

    }
    removeEdges(tree, L->pos);
 
  }
  else{

    //    cerr << "Joining right" << endl;
    R->collapsed = true;

    for(vector<edge *>::iterator lc =  R->eds.begin(); lc != R->eds.end(); lc++){

      edge * e;

      int otherP = (*lc)->L->pos;

      if(R->pos  == otherP){
	otherP = (*lc)->R->pos;
      }
      if(findEdge(L->eds, &e,  otherP)){
        e->support['I'] += (*lc)->support['I'];
        e->support['D'] += (*lc)->support['D'];
        e->support['S'] += (*lc)->support['S'];
	e->support['H'] += (*lc)->support['H'];
	e->support['L'] += (*lc)->support['L'];
	e->support['M'] += (*lc)->support['M'];
	e->support['X'] += (*lc)->support['X'];
	e->support['V'] += (*lc)->support['V'];
      }
      else{
        if((*lc)->L->pos == otherP){
          (*lc)->R = L;
	}
        else{
          (*lc)->L = L;
        }
        L->eds.push_back(*lc);
	
      }
    }
    removeEdges(tree, R->pos);
  }

}



//------------------------------- SUBROUTINE --------------------------------
/*
  Function input  : a vector of node pointers
  
  Function does   : does a tree terversal and joins nodes
  
  Function returns: NA
  
*/

void collapseTree(vector<node *> & tree){

  //  cerr << "Collapsing" << endl;

  vector<node *> tmp;

  for(vector<node *>::iterator tr = tree.begin(); tr != tree.end(); tr++){
    if((*tr)->collapsed){
//      cerr << "ci" << endl;
      continue;
    }
    for(vector<node *>::iterator tt = tree.begin(); tt != tree.end(); tt++){
      if((*tt)->collapsed){
	//      cerr << "cii" << endl;
        continue;
      }
      if( (*tr)->pos == (*tt)->pos ){
	//      cerr << "ciii" << endl;
        continue;
      }
      //      cerr << "N1: " << (*tr)->pos << " N2: " << (*tt)->pos << endl;

      if(abs( (*tr)->pos - (*tt)->pos ) < 20 && ! connectedNode((*tr), (*tt))){
	
	joinNodes((*tr), (*tt), tree);
	      
      }
    }
  }

  for(vector<node *>::iterator tr = tree.begin(); tr != tree.end(); tr++){
    if((*tr)->collapsed){
    }
    else{
      tmp.push_back((*tr));
    }
  }

  tree.clear();
  tree.insert(tree.end(), tmp.begin(), tmp.end());

}



//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of node pointers

 Function does   : tries to resolve an inversion

 Function returns: NA

*/

bool detectInversion(vector<node *> & tree, breakpoints * bp){

  vector <node *> putative;

  for(vector<node * >::iterator t = tree.begin(); t != tree.end(); t++){


    int tooFar           = 0;
    int splitR           = 0;
    int del              = 0;



    for(vector<edge *>::iterator es = (*t)->eds.begin(); es != (*t)->eds.end(); es++){

      tooFar += (*es)->support['M'] + (*es)->support['H'];
      splitR += (*es)->support['V'];
      del    += (*es)->support['D'];

    }

    if( (tooFar > 0 && splitR > 1) || (del > 1 && splitR > 0) || splitR > 1){
      putative.push_back((*t));
    }
  }
  
  if(putative.size() == 2){

    sort(putative.begin(), putative.end(), sortNodesByPos);
    
    int lPos = putative.front()->pos;
    int rPos = putative.back()->pos ;
    
    int lhit = 0 ; 
    int rhit = 0;
    
    int totalS = 0;

    for(vector<edge *>::iterator ed = putative.front()->eds.begin() ;
        ed != putative.front()->eds.end(); ed++){
      if(((*ed)->L->pos == rPos) || ((*ed)->R->pos == rPos)){
        lhit = 1;
	totalS += (*ed)->support['L'];
	totalS += (*ed)->support['H'];
	totalS += (*ed)->support['S'];
	totalS += (*ed)->support['I'];
	totalS += (*ed)->support['D'];
	totalS += (*ed)->support['V'];
	totalS += (*ed)->support['M'];
	totalS += (*ed)->support['R'];
	totalS += (*ed)->support['X'];
        break;
      }
    }

    for(vector<edge *>::iterator ed = putative.back()->eds.begin() ;
        ed != putative.back()->eds.end(); ed++){
      if(((*ed)->L->pos == lPos) || ((*ed)->R->pos == lPos)){
        rhit = 1;
        break;
      }
    }

    if((rPos - lPos) > 1500000 ){
      if(totalS < 5){
	return  false;
      }
    }

    if(lhit == 1 && rhit == 1){
      bp->two           = true                   ;
      bp->type          = 'V'                    ;
      bp->merged        =  0                     ; 
      bp->seqidIndexL   = putative.front()->seqid;
      bp->seqidIndexR   = putative.front()->seqid;
      bp->five          = lPos                   ;
      bp->three         = rPos                   ;
      bp->svlen         = rPos - lPos            ;
      bp->totalSupport  = totalS                 ;
      for(map<string, int>::iterator iz = putative.front()->sm.begin()
            ; iz != putative.front()->sm.end(); iz++){
        bp->sml.push_back(iz->first);
      }
      for(map<string, int>::iterator iz = putative.back()->sm.begin()
            ; iz != putative.back()->sm.end(); iz++){
        bp->smr.push_back(iz->first);
      }

      bp->supports.push_back(getSupport(putative.front()));
      bp->supports.push_back(getSupport(putative.back()));
      return true;
    }
    else{
      cerr << "no linked putative breakpoints" << endl;
    }
  }
  else if(putative.size() > 2){

    vector <node *> putativeTwo;
    cerr << "greater than two breakpoints: " << putative.size() << " "  << putative.front()->pos  <<\
      endl;

    sort(putative.begin(), putative.end(), sortNodesBySupport);

    while(putative.size() > 2){
      putative.pop_back();
    }

    sort(putative.begin(), putative.end(), sortNodesByPos);

    int lPos = putative.front()->pos;
    int rPos = putative.front()->pos;

    // cerr << lPos << " " << rPos << endl;

    int nhit   = 0;
    int totalS = 0;

    for(vector<edge *>::iterator ed = putative[0]->eds.begin() ;
        ed != putative[0]->eds.end(); ed++){
      if(((*ed)->L->pos == lPos) || ((*ed)->R->pos == lPos)){

        totalS += (*ed)->support['L'];
        totalS += (*ed)->support['H'];
        totalS += (*ed)->support['S'];
        totalS += (*ed)->support['I'];
        totalS += (*ed)->support['D'];
        totalS += (*ed)->support['V'];
        totalS += (*ed)->support['M'];
        totalS += (*ed)->support['R'];
        totalS += (*ed)->support['X'];

        nhit += 1;
        break;
      }
    }

    for(vector<edge *>::iterator ed = putative[1]->eds.begin() ;
	ed != putative[1]->eds.end(); ed++){
      if(((*ed)->L->pos == rPos) || ((*ed)->R->pos == rPos)){
	nhit += 1;
	break;
      }
    }

    if(nhit == 2){

      if((rPos - lPos) > 1500000 ){
        if(totalS < 5){
          return  false;
        }
      }
      bp->two           = true                   ;
      bp->type          = 'V'                    ;
      bp->merged        =  0                     ;
      bp->seqidIndexL   = putative.front()->seqid;
      bp->seqidIndexR   = putative.front()->seqid;
      bp->five          = lPos                   ;
      bp->three         = rPos                   ;
      bp->svlen         = rPos - lPos            ;
      bp->totalSupport  = totalS                 ;

      for(map<string, int>::iterator iz = putative.front()->sm.begin()
	    ; iz != putative.front()->sm.end(); iz++){
        bp->sml.push_back(iz->first);
      }
      for(map<string, int>::iterator iz = putative.back()->sm.begin()
            ; iz != putative.back()->sm.end(); iz++){
        bp->smr.push_back(iz->first);
      }

      bp->supports.push_back(getSupport(putative.front()));
      bp->supports.push_back(getSupport(putative.back()));
      return true;
    }  
  }
  return false;
}




//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of node pointers

 Function does   : tries to resolve a duplication

 Function returns: NA

*/

bool detectDuplication(vector<node *> tree, breakpoints * bp){

  vector <node *> putative;

  for(vector<node * >::iterator t = tree.begin(); t != tree.end(); t++){

    int tooFar   = 0;
    int tooClose = 0;
    int FlippedS = 0;

    for(vector<edge *>::iterator es = (*t)->eds.begin(); es != (*t)->eds.end(); es++){

      tooFar   += (*es)->support['H'];
      tooClose += (*es)->support['L'];
      FlippedS += (*es)->support['X'];

    }
    if( (tooFar > 0 && FlippedS  > 1) || FlippedS  > 1 || (tooClose  > 0 && FlippedS  > 1) ){
      putative.push_back((*t));
    }
  }
  if(putative.size() == 2){

    sort(putative.begin(), putative.end(), sortNodesByPos);

    int lPos = putative.front()->pos;
    int rPos = putative.back()->pos ;

    int lhit = 0 ; int rhit = 0;

    int totalS = 0;

    for(vector<edge *>::iterator ed = putative.front()->eds.begin() ;
        ed != putative.front()->eds.end(); ed++){
      if(((*ed)->L->pos == rPos) || ((*ed)->R->pos == rPos)){
        lhit = 1;
	totalS += (*ed)->support['L'];
        totalS += (*ed)->support['H'];
        totalS += (*ed)->support['S'];
        totalS += (*ed)->support['I'];
        totalS += (*ed)->support['D'];
        totalS += (*ed)->support['V'];
        totalS += (*ed)->support['M'];
        totalS += (*ed)->support['R'];
        totalS += (*ed)->support['X'];

        break;
      }
    }

    for(vector<edge *>::iterator ed = putative.back()->eds.begin() ;
        ed != putative.back()->eds.end(); ed++){
      if(((*ed)->L->pos == lPos) || ((*ed)->R->pos == lPos)){
        rhit = 1;
        break;
      }
    }
    
    if((rPos - lPos) > 1500000 ){
      if(totalS < 5){
	return  false;
      }
    }

    if(lhit == 1 && rhit == 1){
      bp->two         = true                   ;
      bp->type        = 'U'                    ;
      bp->merged      = 0                      ;
      bp->seqidIndexL = putative.front()->seqid;
      bp->seqidIndexR = putative.front()->seqid;
      bp->five        = lPos                   ;
      bp->three       = rPos                   ;
      bp->svlen       = rPos - lPos            ;
      bp->totalSupport  = totalS               ;

      for(map<string, int>::iterator iz = putative.front()->sm.begin()
	    ; iz != putative.front()->sm.end(); iz++){
	bp->sml.push_back(iz->first);
      }
      for(map<string, int>::iterator iz = putative.back()->sm.begin()
            ; iz != putative.back()->sm.end(); iz++){
	bp->smr.push_back(iz->first);
      }

      bp->supports.push_back(getSupport(putative.front()));
      bp->supports.push_back(getSupport(putative.back()));
      
      return true;
    }
    else{
      cerr << "no linked putative breakpoints" << endl;
    }
    
  }
  else if(putative.size() > 2){

    vector <node *> putativeTwo;
    cerr << "greater than two breakpoints: " << putative.size() << " "  << putative.front()->pos  << endl;
  
    sort(putative.begin(), putative.end(), sortNodesBySupport);

    while(putative.size() > 2){
      putative.pop_back();
    }

    sort(putative.begin(), putative.end(), sortNodesByPos);


    int lPos = putative.front()->pos;
    int rPos = putative.back()->pos;


    int nhit   = 0;
    int totalS = 0;

    for(vector<edge *>::iterator ed = putative[0]->eds.begin() ;
	ed != putative[0]->eds.end(); ed++){
      if(((*ed)->L->pos == lPos) || ((*ed)->R->pos == lPos)){
	
	totalS += (*ed)->support['L'];
	totalS += (*ed)->support['H'];
	totalS += (*ed)->support['S'];
	totalS += (*ed)->support['I'];
	totalS += (*ed)->support['D'];
	totalS += (*ed)->support['V'];
	totalS += (*ed)->support['M'];
	totalS += (*ed)->support['R'];
	totalS += (*ed)->support['X'];


	nhit += 1;
	break;
      }
    }
    
    if(nhit == 2){
      if((rPos - lPos) > 1500000 ){
	if(totalS < 5){
	  return  false;
	}
      }
      
      bp->two         = true                   ;
      bp->type        = 'U'                    ;
      bp->merged      = 0                      ;
      bp->seqidIndexL = putative.front()->seqid;
      bp->seqidIndexR = putative.front()->seqid;
      bp->five        = lPos                   ;
      bp->three       = rPos                   ;
      bp->svlen       = rPos - lPos            ;
      bp->totalSupport  = totalS               ;

      for(map<string, int>::iterator iz = putative.front()->sm.begin()
            ; iz != putative.front()->sm.end(); iz++){
        bp->sml.push_back(iz->first);
      }
      for(map<string, int>::iterator iz = putative.back()->sm.begin()
            ; iz != putative.back()->sm.end(); iz++){
        bp->smr.push_back(iz->first);
      }

      bp->supports.push_back(getSupport(putative.front()));
      bp->supports.push_back(getSupport(putative.back()));
      return true;
    }

  }
  else{
    // leaf node
  }
  
  return false;
}



//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of node pointers

 Function does   : tries to resolve a deletion with no links

 Function returns: NA

*/

bool detectHalfDeletion(vector<node *> & tree, breakpoints * bp, node ** n){
  
  vector <node *> putative;
  vector <int>    support ;

  for(vector<node * >::iterator t = tree.begin(); t != tree.end(); t++){

    int tooFar  = 0;
    int splitR  = 0;
    int del     = 0;

    for(vector<edge *>::iterator es = (*t)->eds.begin(); es != (*t)->eds.end(); es++){
      tooFar += (*es)->support['H'];
      splitR += (*es)->support['S'];
      del    += (*es)->support['D'];
    }


    if( tooFar > 2 || (splitR > 0 && tooFar > 0)){
      putative.push_back((*t));
      support.push_back(tooFar);
    }
  }
  


  if(putative.size() >= 1){

    int max   = support.front();
    int index = 0;
    int maxi  = 0;

    for(vector<int>::iterator it = support.begin(); 
	it != support.end();  it++){
      
      if(*it > max){
	max = *it;
	maxi = index;
      }
      index+=1;
    }
    
    bp->seqidIndexL = putative[maxi]->seqid;
    bp->five        = putative[maxi]->pos  ;
    *n = putative[maxi];

    return true;
  }

  return false;
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of node pointers

 Function does   : tries to resolve a deletion

 Function returns: NA

*/


bool detectDeletion(vector<node *> tree, breakpoints * bp){
  
  vector <node *> putative;

  for(vector<node * >::iterator t = tree.begin(); t != tree.end(); t++){
   
      int tooFar  = 0; 
      int splitR  = 0;
      int del     = 0;

      for(vector<edge *>::iterator es = (*t)->eds.begin(); es != (*t)->eds.end(); es++){
	tooFar += (*es)->support['H'];
	splitR += (*es)->support['S'];
	del    += (*es)->support['D'];
      }
      if( (tooFar > 0 && splitR > 1) || (del > 1 && splitR > 0) || splitR > 1){
	putative.push_back((*t));
      } 
  }
  
  if(putative.size() == 2){

    sort(putative.begin(), putative.end(), sortNodesByPos);
    
    int lPos = putative.front()->pos;
    int rPos = putative.back()->pos ;
    
    int lhit = 0 ; int rhit = 0;
    
    int totalS = 0;

    for(vector<edge *>::iterator ed = putative.front()->eds.begin() ;
	ed != putative.front()->eds.end(); ed++){
      if(((*ed)->L->pos == rPos) || ((*ed)->R->pos == rPos)){
	lhit = 1; 
	totalS += (*ed)->support['L'];
        totalS += (*ed)->support['H'];
        totalS += (*ed)->support['S'];
        totalS += (*ed)->support['I'];
        totalS += (*ed)->support['D'];
        totalS += (*ed)->support['V'];
        totalS += (*ed)->support['M'];
        totalS += (*ed)->support['R'];
        totalS += (*ed)->support['X'];

	break;
      }
    }
    
    for(vector<edge *>::iterator ed = putative.back()->eds.begin() ;
        ed != putative.back()->eds.end(); ed++){
      if(((*ed)->L->pos == lPos) || ((*ed)->R->pos == lPos)){
        rhit = 1;
	break;
      }
    }  

    if((rPos - lPos) > 1500000 ){
      if(totalS < 5){
	return  false;
      }
    }


    if(lhit == 1 && rhit == 1){
      bp->two         = true                   ;
      bp->type        = 'D'                    ;
      bp->seqidIndexL = putative.front()->seqid;
      bp->seqidIndexR = putative.front()->seqid;
      bp->merged      = 0                      ;
      bp->five        = lPos + 1               ;  // always starts after last node
      bp->three       = rPos - 1               ;
      bp->svlen       = rPos - lPos            ;
      bp->totalSupport = totalS                ;
      for(map<string, int>::iterator iz = putative.front()->sm.begin()
            ; iz != putative.front()->sm.end(); iz++){
        bp->sml.push_back(iz->first);
      }
      for(map<string, int>::iterator iz = putative.back()->sm.begin()
            ; iz != putative.back()->sm.end(); iz++){
        bp->smr.push_back(iz->first);
      }

      bp->supports.push_back(getSupport(putative.front()));
      bp->supports.push_back(getSupport(putative.back()));
      return true;
    }
    else{
      cerr << "no linked putative breakpoints" << endl;
    }

  }
  else if(putative.size() > 2){

    vector <node *> putativeTwo;
       cerr << "greater than two breakpoints: " << putative.size() << " "  << putative.front()->pos  << endl;
  
       sort(putative.begin(), putative.end(), sortNodesBySupport);

       while(putative.size() > 2){
	 putative.pop_back();
       }

       sort(putative.begin(), putative.end(), sortNodesByPos);

       int lPos = putative.front()->pos;
       int rPos = putative.back()->pos;
       
       int nhit   = 0;
       int totalS = 0;

       for(vector<edge *>::iterator ed = putative[0]->eds.begin() ;
	   ed != putative[0]->eds.end(); ed++){
	 if(((*ed)->L->pos == lPos) || ((*ed)->R->pos == lPos)){

	   totalS += (*ed)->support['L'];
	   totalS += (*ed)->support['H'];
	   totalS += (*ed)->support['S'];
	   totalS += (*ed)->support['I'];
	   totalS += (*ed)->support['D'];
	   totalS += (*ed)->support['V'];
	   totalS += (*ed)->support['M'];
	   totalS += (*ed)->support['R'];
	   totalS += (*ed)->support['X'];


	   nhit += 1;
	   break;
	 }
       }

       for(vector<edge *>::iterator ed = putative[1]->eds.begin() ;
           ed != putative[1]->eds.end(); ed++){
         if(((*ed)->L->pos == rPos) || ((*ed)->R->pos == rPos)){
           nhit += 1;
           break;
         }
       }

       if(nhit == 2){

	 if((rPos - lPos) > 1500000 ){
	   if(totalS < 5){
	     return  false;
	   }
	 }
	 bp->two         = true                   ;
	 bp->type        = 'D'                    ;
	 bp->seqidIndexL = putative.front()->seqid;
	 bp->seqidIndexR = putative.front()->seqid;
	 bp->merged      = 0                      ;
	 bp->five        = lPos + 1               ;  // always starts after last node
	 bp->three       = rPos - 1               ;
	 bp->svlen       = rPos - lPos            ;
	 bp->totalSupport = totalS                ;
	 for(map<string, int>::iterator iz = putative.front()->sm.begin()
	       ; iz != putative.front()->sm.end(); iz++){
	   bp->sml.push_back(iz->first);
	 }
	 for(map<string, int>::iterator iz = putative.back()->sm.begin()
	       ; iz != putative.back()->sm.end(); iz++){
	   bp->smr.push_back(iz->first);
	 }
	 bp->supports.push_back(getSupport(putative.front()));
	 bp->supports.push_back(getSupport(putative.back()));
	 return true;


       }

       
  }
  else{
    // leaf node
  }
  
  return false;
}

void callBreaks(vector<node *> & tree, 
		vector<breakpoints *> & allBreakpoints, 
		map < int , map <int, node *> > & hb){

  collapseTree(tree);

  breakpoints * bp;

  node * nr; 

  bp = new breakpoints;
  bp->fail    = false;
  bp->two     = false;
  bp->refined = 0;
  bp->lalt = 0;
  bp->lref = 0;
  bp->collapsed = 0;

  char hex[8 + 1];
  for(int i = 0; i < 8; i++) {
    sprintf(hex + i, "%x", rand() % 16);
  }

  stringstream xx ;
  xx << hex;
  bp->id = xx.str();

  if(detectDeletion(tree, bp)){
    omp_set_lock(&lock);
    allBreakpoints.push_back(bp);
    omp_unset_lock(&lock);
  }
  else if(detectDuplication(tree, bp)){
    omp_set_lock(&lock);
    allBreakpoints.push_back(bp);
    omp_unset_lock(&lock);
  }
  else if(detectInversion(tree, bp)){
    omp_set_lock(&lock);
    allBreakpoints.push_back(bp);
    omp_unset_lock(&lock);
  }
  else if(detectHalfDeletion(tree, bp, &nr)){
    omp_set_lock(&lock);
    hb[bp->seqidIndexL][bp->five] = nr;
    delete bp;
    omp_unset_lock(&lock);
  }
  else{
    delete bp;

  }
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : processes trees

 Function does   : tries to define type and breakpoints

 Function returns: NA

*/

void gatherTrees(vector<vector<node *> > & globalTrees){
  
  map<int, map<int, int> > lookup;
  
  for(map<int, map<int, node* > >::iterator it = globalGraph.nodes.begin();it != globalGraph.nodes.end(); it++){
    for(map<int, node*>::iterator itt = it->second.begin(); itt != it->second.end(); itt++){

      if(lookup[it->first].find(itt->first) != lookup[it->first].end() ){
      }
      else{
	lookup[it->first][itt->first] = 1;
	vector<node *> tree;
	getTree(globalGraph.nodes[it->first][itt->first], tree);
	for(vector<node *>::iterator ir = tree.begin(); ir != tree.end(); ir++){
          lookup[(*ir)->seqid][(*ir)->pos] = 1;
        }
	globalTrees.push_back(tree);
      }
    }
  }
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : nothing

 Function does   : dumps and shrinks graph

 Function returns: NA

*/


void dump(vector< vector< node *> > & allTrees){
  
  ofstream graphOutFile;

  graphOutFile.open(globalOpts.graphOut);

  for(vector< vector<node *> >::iterator it = allTrees.begin(); 
      it != allTrees.end(); it++){
    graphOutFile << dotviz(*it) << endl << endl;
  }

  graphOutFile.close();
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : bam file name and breakpoint

 Function does   : aligns reads to breakpoints

 Function returns: NA
*/

double pBases(vector<unsigned int> & cigar, string & baseQs){
  
  double phredSum   = 0;
  int runningPos = 0; 

  for(int i = 0; i < cigar.size(); i++){
    
    int op = cigar[i]&BAM_CIGAR_MASK  ;
    int l  = cigar[i]>>BAM_CIGAR_SHIFT  ;

    for(int j = 0; j < l; j++){      

      // match M
      if(op == 0 ){
	phredSum += log( 1 - unPhred(SangerLookup[int(baseQs[runningPos])])) ;
	runningPos += 1;
	continue;
      }
      // insertion I
      if(op == 1 ){
	phredSum += log(     unPhred(SangerLookup[int(baseQs[runningPos])])) ;
	runningPos += 1;
	continue;
      }
      // deletion D
      if(op == 2 ){
	phredSum += -3;
	//phredSum += -4.60517;
	//phredSum += -10;

	continue;
      }
      // skip N
      if(op == 3 ){
	runningPos += 1;
	continue;
      }
      // softclip S
      if(op == 4 ){
	phredSum += log(unPhred(SangerLookup[int(baseQs[runningPos])])) ;
	runningPos += 1;
	continue;
      }
      // hardclip H
      if(op == 5 ){
	cerr << "FATAL: hard clip made it into genotyping." << endl;
	cerr << "INFO : report bug to author.             " << endl;
	exit(1);
      }    
      // padding P
      if(op == 6 ){
        cerr << "FATAL: padded seq made it into genotyping." << endl;
        cerr << "INFO : report bug to author.             " << endl;
        exit(1);
      }
      // match =
      if(op == 7 ){

	cerr << baseQs[runningPos] << " " << SangerLookup[int(baseQs[runningPos])] << endl;

	phredSum += log( 1 - unPhred(SangerLookup[int(baseQs[runningPos])])) ;

	runningPos += 1;
	continue;
      }
      // mismatch X
      if(op == 8 ){
	phredSum += log(unPhred(SangerLookup[int(baseQs[runningPos])])) ;
	runningPos += 1;
	continue;
      }
    }
  }
  return phredSum;
}

void genotype(string & bamF, breakpoints * br){

  vector<string> bamFiles;

  bamFiles.push_back(bamF);

  vector<BamAlignment> reads;

  int buffer = 5;
    
  while(reads.size() < 2){
    getPopAlignments(bamFiles, br, reads, buffer);
    buffer = buffer + 1;     
  }

  alignHMM refHMM(int(reads.front().Length) +1,int(br->alleles.front().size()) +1);
  alignHMM altHMM(int(reads.front().Length) +1,int(br->alleles.back().size()) +1);

  double aal = 0;
  double abl = 0;
  double bbl = 0;

  int nref   = 0;
  int nalt   = 0;

  int nReads = 0;

  phredUtils pu;

  for(vector<BamAlignment>::iterator it = reads.begin(); it != reads.end(); it++){

    if(endsBefore(*it, br->five,20) || startsAfter(*it, br->three,20)){
      continue;
    }
    
    nReads += 1;
    

    refHMM.initPriors(br->alleles.front(), it->QueryBases, it->Qualities);
    refHMM.initTransProbs();
    refHMM.initializeDelMat();
    refHMM.updatecells();
    
    double pR = refHMM.finalLikelihoodCalculation();


    double div2 = log10(2);

    altHMM.initPriors(br->alleles.back(), it->QueryBases, it->Qualities);
    altHMM.initTransProbs();
    altHMM.initializeDelMat();
    altHMM.updatecells();
    double pA = altHMM.finalLikelihoodCalculation();

#ifdef DEBUG
    
    cerr << "RI: " << (*it).Position << " "
	 << (*it).GetEndPosition() << " "
	 << joinCig((*it).CigarData) << " "
	 << (*it).MapQuality         << endl;
    cerr << "pref palt:" << pR << " "<< pA << endl;

    cerr << br->alleles.front() << " REF " << it->QueryBases << " " << it->Qualities << " " << pR <<  " " << pR << endl; 
    cerr << br->alleles.back() << " ALT " << it->QueryBases << " " << it->Qualities << " " << pA <<  " " << pR << endl; 
    
    cerr << endl;

#endif

	      

    if(br->type == 'V'){

      string revcomp = string((*it).QueryBases.rbegin(), (*it).QueryBases.rend());
      Comp(revcomp);
      string revqual = string((*it).Qualities.rbegin(), (*it).Qualities.rend()); 

      altHMM.initPriors(br->alleles.back(), revcomp, revqual);
      altHMM.initTransProbs();
      altHMM.initializeDelMat();
      altHMM.updatecells();
      
      double pA1 =  altHMM.finalLikelihoodCalculation();
      
      if(pA1 > pA){
	pA = pA1;

	refHMM.initPriors(br->alleles.front(), revcomp, revqual);
	refHMM.initTransProbs();
	refHMM.initializeDelMat();
	refHMM.updatecells();

	pR = refHMM.finalLikelihoodCalculation();

      }
    }

    if(pR > pA){
      nref += 1;
    }
    else{
      nalt += 1;
    }

#ifdef DEBUG
    cerr << "ll: " << aal << " " << abl << " " << bbl << endl;
#endif 

    aal  += pu.log10Add((pR - div2), (pR - div2));
    abl  += pu.log10Add((pR - div2), (pA - div2));
    bbl  += pu.log10Add((pA - div2), (pA - div2));
    


  }



  int nnref = nref;
  int nnalt = nalt;

  while(nnref > 100 || nnalt > 100){
    nnref -= 1;
    nnalt -= 1;
    if(nnref < 0){
      nnref = 0;
    }
    if(nnalt < 0){
      nnalt = 0;
    }
  }
  


  vector<double> gl;
  gl.push_back(aal);
  gl.push_back(abl);
  gl.push_back(bbl);
 
  int index = -1;

  if(nref > 0 || nalt > 0){
    index = 0;
  }

  if(abl > aal && abl > bbl){
    index = 1;
  }
  if(bbl > aal && bbl > abl){
    index = 2;
  }

  br->genotypeLikelhoods.push_back(gl);
  br->genotypeIndex.push_back(index);
  br->nref.push_back(nref);
  br->nalt.push_back(nalt);


  br->lref += aal + (abl/2);
  br->lalt += bbl + (abl/2);

}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : bam file

 Function does   : dumps and shrinks graph

 Function returns: NA
*/
  
void gatherBamStats(string & targetfile){

  omp_set_lock(&lock);
  int quals [126];
  memcpy(quals, SangerLookup, 126*sizeof(int));

  FastaReference RefSeq;
  RefSeq.open(globalOpts.fasta);

  cerr << "INFO: gathering stats (may take some time) for bam: " << targetfile << endl;
  
  omp_unset_lock(&lock);

  vector<double> alIns        ;
  vector<double> nReads       ;
  vector<double> randomSWScore;

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

  // Declares a default Aligner
  StripedSmithWaterman::Aligner aligner;
  // Declares a default filter
  StripedSmithWaterman::Filter filter;
  // Declares an alignment that stores the result
  StripedSmithWaterman::Alignment alignment;
  // Aligns the query to the ref

  RefVector sequences = bamR.GetReferenceData();

  int i = 0; // index for while loop
  int n = 0; // number of reads

  BamAlignment al;

  int qsum = 0;
  int qnum = 0;

  int fail = 0;

  while(i < 8 || n < 100000){

    if((n % 100) == 0){
      omp_set_lock(&lock);
      cerr << "INFO: processed " << n << " reads for: " << targetfile << endl;
      omp_unset_lock(&lock);
    }

    fail += 1;
    if(fail > 1000000){
      cerr << "FATAL: was not able to gather stats on bamfile: " << targetfile << endl;
      exit(1);
    }

    unsigned int max = 20;

    if(sequences.size() < max){
      max = sequences.size() ;
    }

    int randomChr = 0;
    if(sequences.size() > 1){
      randomChr = rand() % (max -1);
    }
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
      
//      int randomChr2 = rand() % (max -1);
//      int randomPos2 = rand() % (sequences[randomChr].RefLength -1);
//      int randomEnd2 = randomPos + 10000;

//      int randomL = ( (rand() % 200) - 200 ) + al.Position;

//      while(randomEnd2 > sequences[randomChr].RefLength){
//	randomChr2 = rand() % (max -1);
//	randomPos2 = rand() % (sequences[randomChr].RefLength -1);
//	randomEnd2 = randomPos + 10000;
//      }

//      string RefChunk = RefSeq.getSubSequence(sequences[randomChr].RefName, randomL, 200);
					      
//      aligner.Align(al.QueryBases.c_str(), RefChunk.c_str(), RefChunk.size(), filter, &alignment);      

//      if(alignment.sw_score > 0){
//	randomSWScore.push_back(double(alignment.sw_score));
//      }      

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
  
  sort(alIns.begin(), alIns.end()     );

  int index = 0;

  if((alIns.size() % 2) != 0 ){
    index = ((alIns.size()+1)/2)-1;
  }
  else{
    index = (alIns.size()/2)-1;
  }

  double median   = alIns[index];
  double mu       = mean(alIns        );
  double mud      = mean(nReads       );
  double variance = var(alIns, mu     );
  double sd       = sqrt(variance     );
  double sdd      = sqrt(var(nReads, mud ));
  
  //  double sw_mu    = mean(randomSWScore);
  //  double sw_sd    = sqrt(var(randomSWScore, sw_mu));
 
 omp_set_lock(&lock);

 insertDists.mus[  targetfile ] = mu;
 insertDists.sds[  targetfile ] = sd;
 insertDists.avgD[ targetfile ] = mud;
 // insertDists.swm[ targetfile ] = sw_mu;
 // insertDists.sws[ targetfile ] = sw_sd;
 
 // cerr << "dist: " << join(randomSWScore, ",") << endl;

 cerr << "INFO: for file:" << targetfile << endl
      << "      " << targetfile << ": mean depth: ......... " << mud << endl
      << "      " << targetfile << ": sd depth: ........... " << sdd << endl
   //      << "      " << targetfile << ": mean SW alignments .. " << sw_mu << endl
   //      << "      " << targetfile << ": sd SW alignments .... " << sw_sd << endl
      << "      " << targetfile << ": mean insert length: . " << insertDists.mus[targetfile] << endl
      << "      " << targetfile << ": median insert length. " << median                      << endl
      << "      " << targetfile << ": sd insert length .... " << insertDists.sds[targetfile] << endl
      << "      " << targetfile << ": lower insert length . " << insertDists.mus[targetfile] - (2.5*insertDists.sds[targetfile])   << endl
      << "      " << targetfile << ": upper insert length . " << insertDists.mus[targetfile] + (2.5*insertDists.sds[targetfile])   << endl
      << "      " << targetfile << ": average base quality: " << double(qsum)/double(qnum) << " " << qsum << " " << qnum << endl
      << "      " << targetfile << ": number of reads used: " << n  << endl << endl;

  omp_unset_lock(&lock);
}


int avgP(node * n){
  
  double pi  = 0;
  double ni  = 0;

  for(vector<edge *>::iterator it = n->eds.begin(); it != n->eds.end(); it++){
    
    if((*it)->L->pos == n->pos){  
      pi += double((*it)->R->pos) * double((*it)->support['H']);
      ni += double((*it)->support['H'])                        ;  
    }
    else{
      pi += double((*it)->L->pos) * double((*it)->support['H']);
      ni += double((*it)->support['H'])                        ;  
    }
    
  }

  return int(double(pi) / double(ni));

}

void mergeDels(map <int, map <int, node * > > & hf, vector< breakpoints *> & br){

  map<int, map<int, int> > seen;

  for(map <int, map <int, node * > >::iterator hfs = hf.begin(); 
      hfs != hf.end(); hfs++){

    for(map<int, node*>::iterator hpos = hfs->second.begin(); 
	hpos != hfs->second.end(); hpos++){
    
      if(seen[hpos->second->seqid].find(hpos->second->pos) != seen[hpos->second->seqid].end()){
	continue;
      }

      int otherPos = avgP(hpos->second);
    
      for(map<int, node*>::iterator spos = hfs->second.begin(); spos != hfs->second.end(); spos++){

	if(hpos->second->pos == spos->second->pos ){
	  continue;
	}
	
	if(seen[spos->second->seqid].find(spos->second->pos) != seen[spos->second->seqid].end()){
	  continue;
	}

	int otherPosSecond = avgP(spos->second);
		
	if(abs(otherPos - spos->second->pos) < 500 && abs(hpos->second->pos - otherPosSecond) < 500 ){

	  seen[spos->second->seqid][spos->second->pos] = 1; 
	  seen[spos->second->seqid][hpos->second->pos] = 1; 

	  vector<node *> putative;
	  putative.push_back(hpos->second);
	  putative.push_back(spos->second);

	  sort(putative.begin(), putative.end(), sortNodesByPos);

	  int lPos = putative.front()->pos;
	  int rPos = putative.back()->pos;
	  

	  cerr << "INFO: joining deletion breakpoints: " <<  lPos << " " << rPos << endl;
	  
	  breakpoints * bp = new breakpoints;

	  char hex[8 + 1];
	  for(int i = 0; i < 8; i++) {
	    sprintf(hex + i, "%x", rand() % 16);
	  }

	  stringstream xx ;
	  xx << hex;
	  bp->id = xx.str();

	  
	  bp->fail         = false                  ;
	  bp->two          = true                   ;
	  bp->type         = 'D'                    ;
	  bp->seqidIndexL  = hpos->second->seqid    ;
	  bp->seqidIndexR  = hpos->second->seqid    ;
	  bp->merged       = 1                      ;
	  bp->five         = lPos                   ;
	  bp->three        = rPos                   ;
	  bp->svlen        = rPos - lPos            ;
	  bp->totalSupport = 0                      ;
	  bp->collapsed    = 0                      ; 
	  for(map<string, int>::iterator iz = putative.front()->sm.begin()
		; iz != putative.front()->sm.end(); iz++){
	    bp->sml.push_back(iz->first);
	  }
	  for(map<string, int>::iterator iz = putative.back()->sm.begin()
		; iz != putative.back()->sm.end(); iz++){
	    bp->smr.push_back(iz->first);
	  }
	  bp->supports.push_back(getSupport(putative.front())) ;
	  bp->supports.push_back(getSupport(putative.back()))  ;

	  br.push_back(bp);
	}		
      }
    }
  }
}


//-------------------------------    MAIN     --------------------------------
/*
 Comments:
*/

int main( int argc, char** argv)
{
  globalOpts.nthreads = -1;
  globalOpts.statsOnly = false;
  globalOpts.skipGeno  = false;
  globalOpts.MQ        = 20   ;

  int parse = parseOpts(argc, argv);
  if(parse != 1){
    cerr << "FATAL: unable to parse command line correctly. Double check commands." << endl;
    cerr << endl;
    printHelp();
    exit(1);
  }


  if(globalOpts.nthreads == -1){
  }
  else{
    omp_set_num_threads(globalOpts.nthreads);
  }

  if(globalOpts.fasta.empty()){
    cerr << "FATAL: no reference fasta provided." << endl << endl;
    printHelp();
    exit(1);
  }

  // gather the insert length and other stats

#pragma omp parallel for schedule(dynamic, 3)
  for(unsigned int i = 0; i < globalOpts.targetBams.size(); i++){     
    gatherBamStats(globalOpts.targetBams[i]);    
  }

  if(globalOpts.statsOnly){
    cerr << "INFO: Exiting as -s flag is set." << endl;
    cerr << "INFO: WHAM finished normally, goodbye! " << endl;
    return 0;
  }

 RefVector sequences;

 BamMultiReader mr;
 if(! mr.Open(globalOpts.targetBams)){
   cerr << "FATAL: issue opening all bams to extract header" << endl;
   exit(1);
 }
 else{
   sequences = mr.GetReferenceData();
   mr.Close();
 }


 map<string, int> inverse_lookup;
 int s = 0;
 
 for(vector<RefData>::iterator it = sequences.begin(); 
     it != sequences.end(); it++){
   
   inverse_lookup[(*it).RefName] = s;
   
   s+= 1;
 }
 
 
 // load bam has openMP inside for running regions quickly
 
 if(globalOpts.svs.empty()){
   
   cerr << "INFO: Loading discordant reads into graph." << endl;
   
   for(vector<string>::iterator bam = globalOpts.targetBams.begin();
       bam != globalOpts.targetBams.end(); bam++){
     
     cerr << "INFO: Reading: " << *bam << endl;
     
     loadBam(*bam);
     
     for(map<int, map<int, node * > >::iterator seqid = globalGraph.nodes.begin();
	 seqid != globalGraph.nodes.end(); seqid++){
       
       cerr << "INFO: Number of putative breakpoints for: " << sequences[seqid->first].RefName << ": " 
	    << globalGraph.nodes[seqid->first].size() << endl;
     }
   }
   
   cerr << "INFO: Finished loading reads." << endl;
 } 
 
 vector<breakpoints*> allBreakpoints ;
 vector<vector<node*> > globalTrees  ;
 
 if(globalOpts.svs.empty()){
   cerr << "INFO: Finding trees within forest." << endl;
   
   gatherTrees(globalTrees);
   
#ifdef DEBUG
   dump(globalTrees);
#endif
   
   map<int, map <int, node*> > delBreak;
   
   cerr << "INFO: Finding breakpoints in trees." << endl; 
#pragma omp parallel for schedule(dynamic, 3)
   for(unsigned int i = 0 ; i < globalTrees.size(); i++){
     
     if((i % 1000) == 0){
       omp_set_lock(&glock);     
       cerr << "INFO: Processed " << i << "/" << globalTrees.size() << " trees" << endl;
       omp_unset_lock(&glock);
     }
     
     if(globalTrees[i].size() > 200){
       omp_set_lock(&glock);
       cerr << "WARNING: Skipping tree, too many putative breaks." << endl;
       omp_unset_lock(&glock);
       continue;
     }
     callBreaks(globalTrees[i], allBreakpoints, delBreak);   
   }
   
 cerr << "INFO: Trying to merge deletion breakpoints: " << delBreak.size() << endl;
 
 mergeDels(delBreak, allBreakpoints);
 }
 
 if(! globalOpts.svs.empty()){
   cerr << "INFO: loading external SV calls" << endl;
   loadExternal(allBreakpoints, inverse_lookup);
 }
  
 cerr << "INFO: Gathering alleles." << endl;
 
 int nAlleles = 0;

#pragma omp parallel for
 for(unsigned  int z = 0; z < allBreakpoints.size(); z++){
   
   if(allBreakpoints[z]->fail){
     continue;
   }
   
   genAlleles(allBreakpoints[z], globalOpts.fasta, sequences);
   omp_set_lock(&glock);
   nAlleles += 1;
   if((nAlleles % 100) == 0){
     cerr << "INFO: generated " << nAlleles  << " alleles / " << allBreakpoints.size()  << endl;
   }
   omp_unset_lock(&glock);
 }
 
 if(globalOpts.skipGeno){
   printBEDPE(allBreakpoints, sequences);
   cerr << "INFO: Skipping genotyping: -k " << endl;
   cerr << "INFO: WHAM finished normally, goodbye! " << endl;
   return 0;
 }
 
 nAlleles = 0;
 cerr << "INFO: Refining breakpoints using SW alignments" << endl;
#pragma omp parallel for
 for(unsigned int z = 0; z < allBreakpoints.size(); z++){
   
   if(allBreakpoints[z]->fail){
     continue;
   }
   
   vector<BamAlignment> reads;
   int buffer = 1;
   while(reads.size() < 2){
     getPopAlignments(globalOpts.targetBams, allBreakpoints[z], reads, buffer);
     buffer +=1;
   }   
   double startingScore = totalAlignmentScore(reads, allBreakpoints[z]);
   double  oldScore = startingScore;
   int oldStart     = allBreakpoints[z]->five ;
   int oldEnd       = allBreakpoints[z]->three;
   int flag         = 0;
   
   breakpoints * secondary = new breakpoints;
   secondary->fail = false;
   
   for(int f = -2; f <= 2; f++){
     *secondary = *allBreakpoints[z];
     secondary->five = oldStart;
     secondary->five += f;
     if(secondary->five >= secondary->three){
       continue;
     }
     for(int t = -2; t <= 2; t++){
       secondary->three  = oldEnd;
       secondary->three  += t;
       if(secondary->three  <= secondary->five){
	 continue;
       }
       
       secondary->svlen = (secondary->three - secondary->five); 
       genAlleles(secondary, globalOpts.fasta, sequences);
       double newScore = totalAlignmentScore(reads, secondary);
       
       if(newScore > startingScore){
	 startingScore = newScore;
	 allBreakpoints[z]->five = secondary->five;
	 allBreakpoints[z]->three = secondary->three;
	 allBreakpoints[z]->svlen = (allBreakpoints[z]->three - allBreakpoints[z]->five);
	 flag = 1;
	 allBreakpoints[z]->refined = 1;
       }
     }
   }
   
   delete secondary;
  

   if(flag == 1){
     allBreakpoints[z]->svlen =   allBreakpoints[z]->three - allBreakpoints[z]->five;
//     omp_set_lock(&lock);  
//     cerr << "INFO: refined breakpoint pair: " << endl
//	  << "      "<< oldStart << " -> "  << allBreakpoints[z]->five << endl
//	  << "      "<< oldEnd   << " -> "  << allBreakpoints[z]->three << endl;
//     cerr << "      SW score: " << oldScore << " -> " << startingScore << endl;
//     omp_unset_lock(&lock);  
   }

   omp_set_lock(&glock);
   nAlleles += 1;
   if((nAlleles % 10) == 0){
     cerr << "INFO: refined " << nAlleles  << " breakpoint pairs / " << allBreakpoints.size()  << endl;
   }
   omp_unset_lock(&glock);  
 }

 cerr << "INFO: Finished breakpoints using SW alignments" << endl;

 int NGeno = 0;
 cerr << "INFO: Genotyping SVs." << endl;
#pragma omp parallel for
 for(unsigned int z = 0; z < allBreakpoints.size(); z++){

   if(allBreakpoints[z]->fail){
     continue;
   }

   //#pragma omp parallel for schedule(dynamic, 3)   
   for(unsigned int i = 0 ; i < globalOpts.targetBams.size(); i++){
     genotype(globalOpts.targetBams[i], allBreakpoints[z]        );
   }
   omp_set_lock(&glock);
   NGeno += 1;
   if((NGeno % 100) == 0 && z != 0){
     cerr << "Genotyped: " << NGeno  << "/" << allBreakpoints.size() << " SVs." << endl;
   }
   omp_unset_lock(&glock);
 }

 cerr << "INFO: Sorting "  << allBreakpoints.size() << " putative SVs." << endl;

 printBEDPE(allBreakpoints, sequences);
  
 if(!globalOpts.graphOut.empty()){
   dump(globalTrees);
 }

 cerr << "INFO: WHAM finished normally, goodbye! " << endl;
 return 0;
}
