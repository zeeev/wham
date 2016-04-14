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
#include "join.hpp"
#include "cigar.hpp"
#include "stats.hpp"
#include "dataStructures.hpp"

#include "fastahack/Fasta.h"
#include "ssw_cpp.h"

// openMP - swing that hammer
#include <omp.h>

// bamtools and my headers
#include "api/BamMultiReader.h"
#include "readPileUp.h"

// paired-like hmm
#include "phredUtils.h"
#include "alignHMM.h"

using namespace std;
using namespace BamTools;

struct options{
  std::vector<string> targetBams;
  bool statsOnly                ;
  bool skipGeno                 ;
  bool keepTrying               ;

  bool vcf                      ;
  int MQ                        ;
  int nthreads                  ;
  int lastSeqid                 ;
  string fasta                  ;
  string graphOut               ;
  map<string, int> toSkip       ;
  map<string, int> toInclude    ;
  string seqid                  ;
  vector<int> region            ;
  string svs                    ;
  map<string,string> SMTAGS     ;
  string saT                    ;
}globalOpts;


//GLOBAL STRUCTS
// the forest of graphs
graph        globalGraph;
// library stats (insertsize ...)
libraryStats insertDists;
// global read pair store
map<string, readPair*> globalPairStore;
// seqid->int index
map<string, int> inverse_lookup;



// options

static const char *optString = "c:i:u:b:m:r:a:g:x:f:e:hskz";

// omp lock

omp_lock_t lock;
omp_lock_t glock;

// read depth cuttoff

uint MAXREADDEPTH = 0;

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : a vector of breakpoint pointers and inverse chr->index
 Function does   : loads an external vcf file for genotyping into the
		   breakpoint vector
 Function returns: bool
*/

bool loadExternal(vector<breakpoints *> & br ){

  ifstream featureFile (globalOpts.svs.c_str());

  string line;

  if(featureFile.is_open()){

    while(getline(featureFile, line)){

      vector<string> SV   = split(line, "\t");

      if(SV.front()[0] == '#'){
	continue;
      }

      vector<string> info = split(SV[7], ";");

      map<string, string> infoDat;

      for(vector<string>::iterator iz = info.begin();
	  iz != info.end(); iz++){

	vector<string> kv = split(*iz, "=");

	infoDat[kv.front()] = kv.back();

      }


      if( (infoDat["SVTYPE"].compare("DEL") == 0)
	  || (infoDat["SVTYPE"].compare("DUP") == 0)
	  || (infoDat["SVTYPE"].compare("INV") == 0) ){

	breakpoints * bk = new breakpoints;

	bk->fail = false;
	bk->two  = true ;
	bk->type = 'D';
	bk->refBase = SV[3];

	if((infoDat["SVTYPE"].compare("DUP") == 0)){
	  bk->type = 'U';
	}
	if((infoDat["SVTYPE"].compare("INV") == 0)){
	  bk->type = 'V';
	}
	if(inverse_lookup.find(SV[0]) == inverse_lookup.end() ){
	  cerr << "FATAL: could not find seqid in inverse lookup: "
	       << SV[0] << endl;
	  exit(1);
	}

	vector<string> POS   = split(infoDat["POS"], ",");
	vector<string> SUP   = split(infoDat["SUPPORT"], ",");

	bk->fives            = infoDat["FIVE"];
	bk->threes           = infoDat["THREE"];
	bk->seqidIndexL = inverse_lookup[SV[0]];
	bk->seqidIndexR = inverse_lookup[SV[0]];
	bk->seqid       = SV[0];
	bk->merged      = false;
	bk->refined     = false;

	bk->five        = atoi(POS[0].c_str()) - 1;
	bk->three       = atoi(POS[1].c_str()) - 1;

	// use ID field if present
	bk->id         = (SV[2] == "" || SV[2] == ".") ? infoDat["ID"] : SV[2];
	bk->collapsed  = atoi(infoDat["COLLAPSED"].c_str());
	bk->lalt       = 0;
	bk->lref       = 0;

	bk->sml = split(infoDat["LID"], ",");
	bk->smr = split(infoDat["RID"], ",");

	vector<string> cipos = split(infoDat["CIPOS"], ",");
	vector<string> ciend = split(infoDat["CIEND"], ",");

	bk->posCIL = atoi(cipos.front().c_str());
	bk->posCIH = atoi(cipos.back().c_str());

	bk->endCIL = atoi(ciend.front().c_str());
	bk->endCIH = atoi(ciend.back().c_str());

	bk->supports.push_back(atoi(SUP[0].c_str()));
	bk->supports.push_back(atoi(SUP[1].c_str()));

	if(bk->five > bk->three){
	  cerr << "FATAL: SV starts before it ends: " << line << endl;
	  exit(1);
	}
	bk->svlen = (bk->three - bk->five);
	br.push_back(bk);
      }
      else{
	cerr << "FATAL: loading external breakpoints: unknown type: "
	     << infoDat["SVTYPE"] << endl;
	exit(1);
      }
    }
  }
  featureFile.close();
  return true;
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : bam alignment, pos
 Function does   : adds soft clip length to end of read to see if it overlaps
		   end pos
 Function returns: bool
*/

inline bool endsBefore(BamAlignment & r, int & pos, int b){
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
 Function returns: bool=
*/

inline bool startsAfter(BamAlignment & r, int & pos, int b){
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
 Function input  : reads and a breakpoint
 Function does   : aligns the reads to the breakpoint to get the sum of SW
 Function returns: double
*/

double totalAlignmentScore(map< string, vector<BamAlignment> > & reads,
			   breakpoints * br){

  int sum = 0;
  int n   = 0;

  // Declares a default Aligner
  StripedSmithWaterman::Aligner aligner;
  // Declares a default filter
  StripedSmithWaterman::Filter filter;
  // Declares an alignment that stores the result
  StripedSmithWaterman::Alignment alignment;
  // Aligns the query to the ref

  for(map<string, vector<BamAlignment> >::iterator it = reads.begin();
      it != reads.end(); it++){
    for(vector<BamAlignment>::iterator r = it->second.begin();
	r != it->second.end(); r++){

      if(endsBefore((*r), br->five,10) || startsAfter((*r), br->three,10)){
	continue;
      }

      if((*r).IsMapped()){
	if(((*r).CigarData.front().Type != 'S'
	    && (*r).CigarData.front().Length < 10) &&
	   ((*r).CigarData.back().Type  != 'S' &&
	    (*r).CigarData.back().Length  < 10)){
	  continue;
	}
      }
      n += 1;

      aligner.Align((*r).QueryBases.c_str(), br->alleles.back().c_str(),
		    br->alleles.back().size(),  filter, &alignment);
      sum +=  alignment.sw_score;
    }
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
 Function does   : loads up the reads
 Function returns: total number of reads
*/

int getPopAlignments(vector<string> & bamFiles,
		      breakpoints * br,
		      map< string, vector<BamAlignment> > & reads,
		      int buffer){

  int nreads = 0;


  for(vector<string>::iterator fs = bamFiles.begin();
      fs != bamFiles.end(); fs++){

    reads[*fs];

    BamReader bamR;
    if(! bamR.Open(*fs)){
      cerr << "FATAL: unable to open bamfile: " << *fs << endl;
      cerr << "INFO: if you have a large number of files, check ulimit: NCPU*NBAM. " << endl;
      exit(1);
    }

    if(!bamR.LocateIndex()){
      vector<string> fileName = split(*fs, ".");
      fileName.back() = "bai";
      string indexName = join(fileName, ".");
      if(! bamR.OpenIndex(indexName) ){
	cerr << "FATAL: cannot find bam index." << endl;
	cerr << "INFO: If you have a large number of files, check ulimit: NCPU*NBAM "
	     << endl;
      }
    }

    if(!bamR.SetRegion(br->seqidIndexL, br->five -buffer,
		       br->seqidIndexL, br->five +buffer)){
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

      reads[*fs].push_back(al);
      nreads += 1;
    }
    if(br->two){
      if(!bamR.SetRegion(br->seqidIndexL,
			 br->three-buffer,
			 br->seqidIndexL,
			 br->three+buffer)){
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
	reads[*fs].push_back(al);
	nreads += 1;
      }
    }

    bamR.Close();
  }

  return true;
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : string reference
 Function does   : revcomp
 Function returns: string

*/
inline void Comp(string & seq){

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
  cerr << "       WHAM-GRAPHENING -k -f my.bam -a my.fasta \\ " << endl;
  cerr << "       -g my.graph.out.txt -e M,GL000207.1 2> wham.err > wham.vcf"
       << endl;
  cerr << endl;
  cerr << " Required:  " << endl;
//------------------------------- XXXXXXXXXX --------------------------------

  cerr << "          -f  <STRING>  A sorted and indexed bam file or a list            " << endl;
  cerr << "                          of bams: a.bam,b.bam,...                         " << endl;
  cerr << "          -a  <STRING>  The reference genome (indexed fasta).              " << endl;
  cerr << endl;
  cerr << " Optional:  Recommended flags are noted with : *                           " << endl;
  cerr << "          -s  <FLAG>    Exits the program after the stats are              " << endl;
  cerr << "                        gathered. [false]                                  " << endl;
  cerr << "  *       -k  <FLAG>    Skip genotyping (much faster). [false]             " << endl;
  cerr << "          -g  <STRING>  File to write graph to (very large output). [false]" << endl;
  cerr << "  *|-c    -e  <STRING>  Comma sep. list of seqids to skip [false].         " << endl;
  cerr << "  *|-e    -c  <STRING>  Comma sep. list of seqids to keep [false].         " << endl;
  cerr << "          -r  <STRING>  Region in format: seqid:start-end [whole genome]   " << endl;
  cerr << "  *       -x  <INT>     Number of CPUs to use [all cores].                 " << endl;
  cerr << "          -m  <INT>     Mapping quality filter [20].                       " << endl;
  cerr << "          -b  <STRING>  External file to genotype [false].                 " << endl;
  cerr << "          -i  <STRING>  non standard split read tag [SA]                   " << endl;
  cerr << "          -z  <FLAG>    Sample reads until success. [false]                " << endl;

  cerr << endl;
  cerr << " Output:  " << endl;
  cerr << "        STDERR: Run statistics and bam stats                        " << endl;
  cerr << "        STOUT : SV calls in VCF or BEDPE format                     " << endl;
  cerr << endl;
  cerr << " Details:  " << endl;
  cerr << "        -z  <FLAG>    WHAM-GRAPHENING can fail if does not sample        " << endl;
  cerr << "                      enough reads. This flag prevents WHAM-GRAPHENING   " << endl;
  cerr << "                      from exiting. If your bam header has seqids not in " << endl;
  cerr << "                      the bam (e.g. split by region) use -z.             " << endl;
  cerr << "        -k  <FLAG>    The WHAM-GRAPHENING pipeline can genotype after    " << endl;
  cerr << "                      samples are merged (-b).  This will save time for  " << endl;
  cerr << "                      population level calling.                          " << endl;
  cerr << "        -b  <STRING>  The VCF output of WHAM-GRAPHENING for genotyping.  " << endl;
  cerr << "                      WHAM-GRAPHENING will genotype any BAM file at the  " << endl;
  cerr << "                      positions in the -b file.                          " << endl;
  cerr << "        -i  <STRING>  WHAM-GRAPHENING uses the optional bwa-mem SA tag.  " << endl;
  cerr << "                      Older version of bwa-mem used XP.                  " << endl;
  cerr << "     -e|-c  <STRING>  A list of seqids to include or exclude while       " << endl;
  cerr << "                      sampling insert and depth.  For humans you should  " << endl;
  cerr << "                      use the standard chromosomes 1,2,3...X,Y.          " << endl;
  cerr << endl;

  printVersion();
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of pointers to breakpoints and a bamtools RefVector
 Function does   : prints a vcf format
 Function returns: nada
*/

void printVCF(vector<breakpoints *> & calls, RefVector & seqs){

  int index = 0;

  sort(calls.begin(), calls.end(), sortBreak);

  stringstream header;

  header << "##fileformat=VCFv4.2" << endl;
  header << "##source=WHAM-GRAPHENING:" << VERSION << endl;
  header << "##reference=" << globalOpts.fasta << endl;
  header << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">" << endl;
  header << "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">" << endl;
  header << "##INFO=<ID=ID,Number=1,Type=String,Description=\"Unique hexadecimal identifier\">" << endl;
  header << "##INFO=<ID=SUPPORT,Number=2,Type=Integer,Description=\"Number of reads supporting POS and END breakpoints\">" << endl;
  header << "##INFO=<ID=MERGED,Number=1,Type=Integer,Description=\"SV breakpoints were joined without split read support 0=false 1=true\">" << endl;
  header << "##INFO=<ID=REFINED,Number=1,Type=Integer,Description=\"SV breakpoints were refined based on SW alignment 0=false 1=true\">" << endl;

  header << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">" << endl;
  header << "##INFO=<ID=POS,Number=2,Type=String,Description=\"POS and END\">" << endl;
  header << "##INFO=<ID=FIVE,Number=.,Type=Integer,Description=\"collapsed POS\">" << endl;
  header << "##INFO=<ID=THREE,Number=.,Type=Integer,Description=\"collapsed END\">" << endl;
  header << "##INFO=<ID=LID,Number=.,Type=String,Description=\"POS breakpoint support came from SM, independent of genotype\">" << endl;
  header << "##INFO=<ID=RID,Number=.,Type=String,Description=\"END breakpoint support came from SM, independent of genotype\">" << endl;
  header << "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">" << endl;
  header << "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">" << endl;
  header << "##INFO=<ID=COLLAPSED,Number=1,Type=Integer,Description=\"Number of SV calls merged into record\">" << endl;
  header << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
  header << "##FORMAT=<ID=GL,Number=.,Type=Float,Description=\"Genotype likelihoods comprised of comma separated floating point log10-scaled likelihoods for all possible genotypes given the set of alleles defined in the REF and ALT fields.\">" << endl;
  header << "##FORMAT=<ID=AS,Number=1,Type=Integer,Description=\"Number of reads that align better to ALT allele\">" << endl;
  header << "##FORMAT=<ID=RS,Number=1,Type=Integer,Description=\"Number of reads that align better to REF allele\">" << endl;
  header << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" ;

  for(vector<string>::iterator iz = globalOpts.targetBams.begin();
      iz !=  globalOpts.targetBams.end(); iz++){

    if(globalOpts.SMTAGS.find(*iz) == globalOpts.SMTAGS.end()){
      cerr << "FATAL: could not find SM tag for: " << *iz << endl;
      exit(1);
    }
    header << "\t" << globalOpts.SMTAGS[*iz];
  }

  cout << header.str() << endl;

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

    string type = "BND";
    switch((*c)->type){
    case 'D':
      type = "DEL";
      svlen = -svlen;
      break;
    case 'U':
      type = "DUP";
      break;
    case 'I':
      type = "INS";
      break;
    case 'V':
      type = "INV";
      break;
    default:
      break;
    }

    ss << seqs[(*c)->seqidIndexL].RefName
       << "\t"
       << ((*c)->five + 1)
       << "\t"
       << "WG:" << type << ":" << (*c)->id
       << "\t"
       << (*c)->refBase
       << "\t"
       << "<" << type << ">"
       << "\t"
       << "."
       << "\t"
       << ".";

    ss << "\tSVTYPE=" << type << ";SVLEN=" << svlen << ";ID=" << (*c)->id << ";"
       << "SUPPORT=" << (*c)->supports[0] << "," << (*c)->supports[1] << ";"
       << "MERGED=" << (*c)->merged << ";"
       << "REFINED=" << (*c)->refined << ";"
       << "END=" << ((*c)->three +1) << ";"
       << "POS=" << ((*c)->five + 1) << "," << ((*c)->three +1) << ";" ;

    stringstream tmp1;
    stringstream tmp2;

    if((*c)->fives.empty()){
      tmp1 << (*c)->five  + 1;
      tmp2 << (*c)->three + 1;

      (*c)->fives = tmp1.str();
      (*c)->threes = tmp2.str();
    }
    ss << "FIVE=" << (*c)->fives << ";";
    ss << "THREE=" << (*c)->threes << ";";

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
    ss << "RID=" << SMR << ";" ;
    ss << "CIPOS=" << (*c)->posCIL << "," << (*c)->posCIH << ";";
    ss << "CIEND=" << (*c)->endCIL << "," << (*c)->endCIH << ";";
    ss << "COLLAPSED=" << (*c)->collapsed ;
    ss << "\tGT:GL:AS:RS";

    if((*c)->genotypeIndex.size() != globalOpts.SMTAGS.size()){

      for(int gi = 0; gi < globalOpts.SMTAGS.size(); gi++){
	ss << "\t" << ".:.:.:." ;
      }
      ss << endl;
      cout << ss.str();
    }
    else{
      for(unsigned int i = 0; i < (*c)->genotypeIndex.size(); i++){
	if((*c)->genotypeIndex[i] == -1){
	  ss << "\t" << "./.:" << "."
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
	cerr << "FATAL: printVCF: unknown genotype." << endl;
	exit(1);
	}
      }
      ss << endl;
      cout << ss.str();
    }
  }
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

  /* if something is pushed to the back of the vector it changes the
     positions ! be warned. */

  while(!edges.empty()){

     uint hit = 0;

     if(seen.find(edges.back()->L->pos) != seen.end()
	&& seen.find(edges.back()->R->pos) != seen.end() ){
       hit = 1;
     }
     else if(seen.find(edges.back()->L->pos) == seen.end()
	     && seen.find(edges.back()->R->pos) == seen.end()){
       seen[edges.back()->L->pos] = 1;
       seen[edges.back()->R->pos] = 1;
       ns.push_back(edges.back()->L);
       ns.push_back(edges.back()->R);
       edges.insert(edges.end(),
		    edges.back()->L->eds.begin(),
		    edges.back()->L->eds.end());
       edges.insert(edges.end(),
		    edges.back()->R->eds.begin(),
		    edges.back()->R->eds.end());
     }
     else if(seen.find(edges.back()->L->pos) == seen.end()){

       seen[edges.back()->L->pos] = 1;

       ns.push_back(edges.back()->L);

       edges.insert(edges.end(),
		    edges.back()->L->eds.begin(),
		    edges.back()->L->eds.end() );
     }
     else{
       seen[edges.back()->R->pos] = 1;
       ns.push_back(edges.back()->R);

       edges.insert(edges.end(),
		    edges.back()->R->eds.begin(),
		    edges.back()->R->eds.end() );
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

    ref = rs.getSubSequence(rv[bp->seqidIndexL].RefName, bp->five - 200,
			    bp->svlen + 400 );
    alt = rs.getSubSequence(rv[bp->seqidIndexL].RefName, bp->five - 200,
			    200) +
      // bp->five = first base of deletion -1 last ref base + 1 for fasta
      rs.getSubSequence(rv[bp->seqidIndexL].RefName, bp->three+1, 200);
    // start one after deletion ends

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
    ref = rs.getSubSequence(rv[bp->seqidIndexL].RefName, bp->five -200,
			    bp->svlen + 400);
    alt = rs.getSubSequence(rv[bp->seqidIndexL].RefName, bp->five - 200,
			    bp->svlen + 200)
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
    string inv = rs.getSubSequence(rv[bp->seqidIndexL].RefName, bp->five,
				   (bp->svlen) );
    inv = string(inv.rbegin(), inv.rend());
    Comp(inv);

    ref = rs.getSubSequence(rv[bp->seqidIndexL].RefName, bp->five -200,
			    (bp->svlen + 400)) ;
    alt = rs.getSubSequence(rv[bp->seqidIndexL].RefName, bp->five-200, 200)
      + inv + rs.getSubSequence(rv[bp->seqidIndexL].RefName, bp->three, 200);


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

    alt = alt.substr(0,400) + alt.substr(bp->svlen, 400)
      + alt.substr(alt.size() -400, 400);

  }

  bp->alleles.clear();

  bp->refBase = rs.getSubSequence(rv[bp->seqidIndexL].RefName, (bp->five), 1);

  bp->alleles.push_back(ref) ;
  bp->alleles.push_back(alt) ;

  // upper case alleles
  std::transform(bp->alleles.front().begin(), bp->alleles.front().end(),
		 bp->alleles.front().begin(), ::toupper);
  std::transform(bp->alleles.back().begin(), bp->alleles.back().end(),
		 bp->alleles.back().begin(), ::toupper);

  return true;
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : read pair pointer

 Function does   : -->s s<-- tests if pairs suggests insertion

 Function returns: bool

*/

inline bool isPointIn(readPair * rp){

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
 Function input  : left position, right position, local graph

 Function does   : adds nodes or builds connections

 Function returns: void

*/

void addIndelToGraph(int refIDL, 
		     int refIDR,
		     int l, 
		     int r, 
		     char s, 
		     string & SM){

  omp_set_lock(&glock);

  if( ! isInGraph(refIDL, l, globalGraph)
      &&  ! isInGraph(refIDR, r, globalGraph) ){

    node * nodeL;
    node * nodeR;
    edge * ed   ;

    nodeL = new node;
    nodeR = new node;
    ed    = new edge;

    nodeL->sm[SM] = 1;
    nodeR->sm[SM] = 1;

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
    nodeL->seqid = refIDL;
    nodeR->pos = r;
    nodeR->seqid = refIDR;

    globalGraph.nodes[refIDL][l] = nodeL;
    globalGraph.nodes[refIDR][r] = nodeR;

  }
 else if(isInGraph(refIDL, l, globalGraph)
	 &&  ! isInGraph(refIDR, r, globalGraph)){

   node * nodeR;
   edge * ed;

   nodeR = new node;
   ed    = new edge;

   nodeR->sm[SM] = 1;
  
   nodeR->collapsed = false;

   nodeR->beginSupport = 0;
   nodeR->endSupport   = 0;

   initEdge(ed);
   ed->support[s] += 1;

   nodeR->pos      = r;
   nodeR->seqid    = refIDR;
   ed->L           = globalGraph.nodes[refIDL][l];
   ed->R           = nodeR;

   nodeR->eds.push_back(ed);

   globalGraph.nodes[refIDL][l]->eds.push_back(ed);
   globalGraph.nodes[refIDR][r] = nodeR;

   if(globalGraph.nodes[refIDL][l]->sm.find(SM) == 
      globalGraph.nodes[refIDL][l]->sm.end()){

     globalGraph.nodes[refIDL][l]->sm[SM] =  1;

   }
   else{
     globalGraph.nodes[refIDL][l]->sm[SM] += 1;
   }

 }
 else if(! isInGraph(refIDL, l, globalGraph)
	 &&  isInGraph(refIDR, r, globalGraph)){

   node * nodeL;
   edge * ed;

   nodeL = new node;
   ed    = new edge;

   nodeL->collapsed    = false;
   nodeL->beginSupport = 0;
   nodeL->endSupport   = 0;

   initEdge(ed);
   ed->support[s] += 1;
   nodeL->pos      = l;
   nodeL->seqid    = refIDL;
   ed->R = globalGraph.nodes[refIDR][r];
   ed->L = nodeL;

   nodeL->eds.push_back(ed);

   globalGraph.nodes[refIDR][r]->eds.push_back(ed);
   globalGraph.nodes[refIDL][l] = nodeL;

   if(globalGraph.nodes[refIDR][r]->sm.find(SM) ==
      globalGraph.nodes[refIDR][r]->sm.end()){

     globalGraph.nodes[refIDR][r]->sm[SM] = 1;

   }
   else{
     globalGraph.nodes[refIDR][r]->sm[SM] += 1;
   }
 }
 else{
   uint hit = 0;

   if(globalGraph.nodes[refIDL][l]->sm.find(SM)
      != globalGraph.nodes[refIDL][l]->sm.end()){
     globalGraph.nodes[refIDL][l]->sm[SM] = 1;
   }
   else{
     globalGraph.nodes[refIDL][l]->sm[SM] += 1;
   }
   if(globalGraph.nodes[refIDR][r]->sm.find(SM)
      != globalGraph.nodes[refIDR][r]->sm.end()){
     globalGraph.nodes[refIDR][r]->sm[SM] = 1;
   }
   else{
     globalGraph.nodes[refIDR][r]->sm[SM] += 1;
   }

   for(vector<edge *>::iterator ite
	 = globalGraph.nodes[refIDL][l]->eds.begin();
       ite != globalGraph.nodes[refIDL][l]->eds.end(); ite++){
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
     ne->L =      globalGraph.nodes[refIDL][l];
     ne->R =      globalGraph.nodes[refIDR][r];
     globalGraph.nodes[refIDL][l]->eds.push_back(ne);
     globalGraph.nodes[refIDR][r]->eds.push_back(ne);
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
	addIndelToGraph(ba.RefID, ba.RefID, p, (p + ci->Length), 'I', SM);
	break;
      }
    case 'D':
      {
	hit = true;
	addIndelToGraph(ba.RefID, ba.RefID, p, (p + ci->Length), 'D', SM);
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
 Function input  : read pair to filter
 Function does   : fails poor quality or non-sv read pairs
 Function returns: true = bad, false = good;
*/

inline bool pairFailed(readPair * rp){

  if(rp->al1.IsMapped() && rp->al2.IsMapped()){
    if(rp->al1.Length == rp->al1.CigarData[0].Length
       && rp->al1.CigarData[0].Type == 'M' &&
       rp->al2.Length == rp->al2.CigarData[0].Length
       && rp->al2.CigarData[0].Type == 'M' ){
      return true;
    }
    if(rp->al1.MapQuality < globalOpts.MQ
       && rp->al2.MapQuality < globalOpts.MQ){
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
 Function input  : bam alignment, vector<saTag>
 Function does   : finds links for splitter to put in graph
 Function returns: NA
*/

void splitToGraph(BamAlignment  & al,
		  vector<saTag> & sa,
		  string & SM       ){

  if(!al.IsMapped()){
    return;
  }

  // too many chimeras ... for me
  if(sa.size() > 1){
    return;
  }

  // split read is trimmed on both side skip
  if(sa.front().cig.front().Type == 'S'
     && sa.front().cig.back().Type == 'S'){
    return;
  }

  // alignment is trimmed on both sides skip
  if(al.CigarData.front().Type   == 'S'
     && al.CigarData.back().Type == 'S'){
    return;
  }

  char support = 'S';

  // split reads are on different strands
  if((sa[0].strand && ! al.IsReverseStrand())
     || (! sa[0].strand && al.IsReverseStrand() )){
    support = 'V';
  }


  int start = al.Position    ;
  int end   = sa.front().pos ;

  /* since both sides are not clipped if the back is clipped we know
     the front is not */

  if(al.CigarData.back().Type == 'S'){
    start = al.GetEndPosition(false,true);
  }

  /* we also know that both sides of the split read are not trimmed */

  if(sa.front().cig.back().Type == 'S'){
    endPos(sa[0].cig, &end)  ;
  }

  addIndelToGraph(al.RefID, sa[0].seqid, start, end, support, SM);
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : pointer to readPair;

 Function does   : adds high and low insert sizes to graph. both pairs are
		   treated as mapped

 Function returns: NA

*/

bool deviantInsertSize(readPair * rp, char supportType, string & SM){

  if(! rp->al1.IsMapped() || ! rp->al2.IsMapped() ){
    return false;
  }
  
  //
  if(rp->al1.RefID != rp->al2.RefID){
    return false;
  }

  // both reads are clipped, skip
  if(! IsLongClip(rp->al1.CigarData, 1)
     && ! IsLongClip(rp->al2.CigarData, 1)){
    return false;
  }
  
  int start = rp->al1.Position;
  int end   = rp->al2.Position;

  if(rp->al1.CigarData.front().Type == 'S'
     || rp->al1.CigarData.back().Type == 'S'){

    if(rp->al2.CigarData.back().Type == 'S'){
      end = rp->al2.GetEndPosition();
    }
    if(rp->al1.CigarData.back().Type == 'S'){
      start = rp->al1.GetEndPosition(false,true);
    }
  }
  else{
    if(rp->al1.CigarData.back().Type == 'S'){
      end = rp->al1.GetEndPosition();
    }
    if(rp->al2.CigarData.back().Type == 'S'){
      start = rp->al2.GetEndPosition(false,true);
    }
  }
  
  addIndelToGraph(rp->al2.RefID, rp->al2.RefID, start, end, supportType, SM);

  return true;
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : pointer to readPair ; seqid to index
 Function does   : processes Pair
 Function returns: NA

*/

void processPair(readPair * rp, 
		 double   * low, 
		 double * high, 
		 string & SM){

  string sa1;
  string sa2;

  bool sameStrand = false;

  if(pairFailed(rp)){
    return;
  }

  indelToGraph(rp->al1, SM);
  indelToGraph(rp->al2, SM);

  if( ! IsLongClip(rp->al1.CigarData, 5)
      && ! IsLongClip(rp->al2.CigarData, 5)){
    return;
  }

  if((rp->al1.IsReverseStrand() && rp->al2.IsReverseStrand())
     || (! rp->al1.IsReverseStrand() && ! rp->al2.IsReverseStrand()) ){
    sameStrand = true;
  }
  
  if( abs(rp->al1.InsertSize) > *high){
    if(sameStrand){
      deviantInsertSize(rp, 'M', SM);
    }
    else{
      deviantInsertSize(rp, 'H', SM);
    }
  }

  if(abs(rp->al1.InsertSize) < *low ){
    if(sameStrand){
      deviantInsertSize(rp, 'R', SM);
    }
    else{
      deviantInsertSize(rp, 'L', SM);
    }
  }
  // put the split reads in the graph
  if(rp->al1.GetTag(  globalOpts.saT, sa1)){
    vector<saTag> parsedSa1;
    parseSA(parsedSa1, sa1, globalOpts.saT, inverse_lookup);
    splitToGraph(rp->al1, parsedSa1, SM);
  }
  if(rp->al2.GetTag(  globalOpts.saT, sa2)){
    vector<saTag> parsedSa2;
    parseSA(parsedSa2, sa2, globalOpts.saT, inverse_lookup);
    splitToGraph(rp->al2, parsedSa2, SM);
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
	       string & SM){

  if(seqidIndex > 4){
    return true;
  }

#ifdef DEBUG
  cerr << "running region: " << seqidIndex 
       << ":" << start << "-" << end << endl;
#endif
  
  // local graph;
  graph localGraph;
  // local read pair store
  map<string, readPair *>pairStoreLocal;

  double high = insertDists.mus[filename] + (2.5 * insertDists.sds[filename]);
  double low  = insertDists.mus[filename] - (2.5 * insertDists.sds[filename]);

  if(low < 0){
    low = 100;
  }

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

    if(pairStoreLocal.find(al.Name) != pairStoreLocal.end()){
      pairStoreLocal[al.Name]->count += 1;
      if(pairStoreLocal[al.Name]->flag == 2){
	pairStoreLocal[al.Name]->al1 = al;
      }
      else{
	pairStoreLocal[al.Name]->al2 = al;
      }
      processPair(pairStoreLocal[al.Name], &low, &high, SM);
      delete pairStoreLocal[al.Name];
      pairStoreLocal.erase(al.Name);
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
      pairStoreLocal[al.Name] = rp;
    }
  }

  // close the bam
  br.Close();

  // load lonely reads into the global struct;
  // if it finds a mate in the global store it processes and deletes
  omp_set_lock(&lock);

  for(map<string, readPair *>::iterator rps = pairStoreLocal.begin();
      rps != pairStoreLocal.end(); rps++){

    if(globalPairStore.find(rps->first) != globalPairStore.end()){
      globalPairStore[rps->first]->count += 1;

      if(globalPairStore[rps->first]->flag == 1 && (*rps->second).flag == 1){
	continue;
      }
      if(globalPairStore[rps->first]->flag == 2 && (*rps->second).flag == 2){
	continue;
      }
      if(globalPairStore[rps->first]->flag == 1){
	(*globalPairStore[rps->first]).al2 = (*rps->second).al2;
      }
      else{
      	(*globalPairStore[rps->first]).al1 = (*rps->second).al1;
      }
      processPair(globalPairStore[rps->first], &low, &high, SM);
      delete globalPairStore[rps->first];
      delete pairStoreLocal[rps->first];
      globalPairStore.erase(rps->first);
    }
    else{
      globalPairStore[rps->first] = pairStoreLocal[rps->first];
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
    region    = true;
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
  int seqidIndex = 0;

  if(region){

    int p = start;
    int e = 0;
    for(; (p+1000000) <= end; p += 1000000){
      regionDat * regionInfo = new regionDat;
      regionInfo->seqidIndex = inverse_lookup[regionSID];
      regionInfo->start      = p                      ;
      regionInfo->end        = 1000000 + p            ;
      regions.push_back(regionInfo);
      e = p + 1000000;
    }
    if(e < end){
      regionDat * regionInfo = new regionDat;
      regionInfo->seqidIndex = inverse_lookup[regionSID];
      regionInfo->start      = p                        ;
      regionInfo->end        = end                      ;
      regions.push_back(regionInfo);
    }
  }
  else{
    seqidIndex = 0;

    for(vector< RefData >::iterator sit = sequences.begin(); 
	sit != sequences.end(); sit++){
      int start = 0;

      if(globalOpts.toSkip.find( (*sit).RefName )
	 == globalOpts.toSkip.end() && ((*sit).RefLength > 1000)){
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


  int Mb = 0;

  // running the regions with openMP
#pragma omp parallel for schedule(dynamic, 3)

  for(unsigned int re = 0; re < regions.size(); re++){
    if(! runRegion(bamFile                ,
		   regions[re]->seqidIndex,
		   regions[re]->start     ,
		   regions[re]->end       ,
		   sequences              ,
		   SM                      )){
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
       << globalPairStore.size()
       << " reads that were not processed"
       << endl;

  // cleanup remaining reads
  for(map<string, readPair*>::iterator rps = globalPairStore.begin();
      rps != globalPairStore.end(); rps++){
    delete rps->second;
  }
}

//-------------------------------   OPTIONS   --------------------------------
int parseOpts(int argc, char** argv)
{
  int opt = 0;
  opt = getopt(argc, argv, optString);
  while(opt != -1){
    switch(opt){
    case 'i':
      {
	globalOpts.saT = optarg;

	if(globalOpts.saT != "XP" && globalOpts.saT != "XP"){
	  cerr << "FATAL: only SA and XP optional tags are supported for split reads" << endl;
	  exit(1);
	}

	cerr << "INFO: You are using a non standard split-read tag: " << globalOpts.saT << endl;


	break;
      }
    case 'u':
      {
      cerr << "INFO: You are using a hidden flag." << endl;
      globalOpts.lastSeqid =  atoi(((string)optarg).c_str());
      cerr << "INFO: Random sampling will only go up to: " << globalOpts.lastSeqid << endl;
      break;
      }
    case 'z':
      {
	globalOpts.keepTrying = true;
	cerr << "INFO: WHAM-GRAPHENING will not give up sampling reads: -z set" << globalOpts.svs << endl;
	break;
      }
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
    case 'c':
      {
	vector<string> seqidsToInclude = split(optarg, ",");
	for(unsigned int i = 0; i < seqidsToInclude.size(); i++){
	  globalOpts.toInclude[seqidsToInclude[i]] = 1;
	  cerr << "INFO: WHAM will only sample seqid: " << seqidsToInclude[i] << endl;
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


    for(map<string,int>::iterator iz = L->sm.begin(); iz != L->sm.end(); iz++){

      if(R->sm.find(iz->first) != R->sm.end()){
	R->sm[iz->first] += iz->second;
      }
      else{
	R->sm[iz->first] = iz->second;
      }
    }

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

    for(map<string,int>::iterator iz = R->sm.begin(); iz != R->sm.end(); iz++){

      if(L->sm.find(iz->first) != L->sm.end()){
	L->sm[iz->first] += iz->second;
      }
      else{
	L->sm[iz->first] = iz->second;
      }
    }



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
      continue;
    }
    for(vector<node *>::iterator tt = tree.begin(); tt != tree.end(); tt++){
      if((*tt)->collapsed){
	continue;
      }
      if( (*tr)->pos == (*tt)->pos ){
	continue;
      }

      if(abs( (*tr)->pos - (*tt)->pos ) < 20 
	 && ! connectedNode((*tr), (*tt))
	 && ( (*tr)->seqid == (*tt)->seqid ) ){
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


    sort(putative.begin(), putative.end(), sortNodesBySupport);

    while(putative.size() > 2){
      putative.pop_back();
    }

    sort(putative.begin(), putative.end(), sortNodesByPos);

    int lPos = putative.front()->pos;
    int rPos = putative.back()->pos;

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
  bp->fail      = false;
  bp->two       = false;
  bp->refined   = 0;
  bp->lalt      = 0;
  bp->lref      = 0;
  bp->collapsed = 0;

  bp->posCIL = -10;
  bp->posCIH =  10;
  bp->endCIL = -10;
  bp->endCIH =  10;

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
 Function input  : vector of reads for an individual

 Function does   : genotype the individual with pairHMM

 Function returns: NA

*/

void genotype(vector<BamAlignment> & reads,
	      string & bamF, breakpoints * br){


  bool toohigh = false;

  int max = insertDists.avgD[bamF] + (4 * sqrt(insertDists.avgD[bamF]));

  if(reads.size() > max ){
    toohigh = true;
  }

  double aal = 0;
  double abl = 0;
  double bbl = 0;

  int nref   = 0;
  int nalt   = 0;

  int nReads = 0;

  phredUtils pu;

  // currently no illumina reads greater than 300
  // single allocation for alignment object

  alignHMM refHMM(300,int(br->alleles.front().size()) +1);
  alignHMM altHMM(300,int(br->alleles.back().size())  +1);

  for(vector<BamAlignment>::iterator it = reads.begin();
      it != reads.end(); it++){

    if(endsBefore(*it, br->five,20) || startsAfter(*it, br->three,20)
       || toohigh || (*it).MapQuality < 10){
      continue;
    }

    refHMM.clear(int((*it).Length) +1, int(br->alleles.front().size()) +1);
    altHMM.clear(int((*it).Length) +1, int(br->alleles.back().size())  +1);

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

    if(br->type == 'V'){

      string revcomp = string((*it).QueryBases.rbegin(),
			      (*it).QueryBases.rend());
      Comp(revcomp);
      string revqual = string((*it).Qualities.rbegin(),
			      (*it).Qualities.rend());

      altHMM.initPriors(br->alleles.back(), revcomp, revqual);
      altHMM.initTransProbs();
      altHMM.initializeDelMat();
      altHMM.updatecells();

      double pA1 =  altHMM.finalLikelihoodCalculation();

      refHMM.initPriors(br->alleles.front(), revcomp, revqual);
      refHMM.initTransProbs();
      refHMM.initializeDelMat();
      refHMM.updatecells();

      double pR1 = refHMM.finalLikelihoodCalculation();

      if(pA1 > pA){
	pA = pA1;
      }
      if(pR1 > pR){
	pR = pR1;
      }
    }

    if(pR > pA){
      nref += 1;
    }
    else{
      nalt += 1;
    }

    aal  += pu.log10Add((pR - div2), (pR - div2));
    abl  += pu.log10Add((pR - div2), (pA - div2));
    bbl  += pu.log10Add((pA - div2), (pA - div2));

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
    vector<string> fileName = split(targetfile, ".");
    fileName.back() = "bai";
    string indexName = join(fileName, ".");
    if(! bamR.OpenIndex(indexName) ){
      cerr << "FATAL: cannot find bam index." << endl;
    }
  }

  SamHeader SH = bamR.GetHeader();
  if(!SH.HasSortOrder()){
  cerr << "FATAL: sorted bams must have the @HD SO: tag in each SAM header." << endl;
    exit(1);
  }

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

  omp_set_lock(&lock);

  globalOpts.SMTAGS[targetfile] = SM;

  omp_unset_lock(&lock);

  RefVector sequences = bamR.GetReferenceData();

  int i = 0; // index for while loop
  int n = 0; // number of reads

  BamAlignment al;

  int qsum = 0;
  int qnum = 0;

  int fail = 0;

  while(i < 8 || n < 100000){

    if((n % 10000) == 0 && fail < 10){
      omp_set_lock(&lock);
      cerr << "INFO: processed " << n << " reads for: " << targetfile << endl;
      omp_unset_lock(&lock);
    }

    fail += 1;
    if(fail > 1000000 && (! globalOpts.keepTrying) ){
      cerr << "FATAL: Unable to gather stats on bamfile: " << targetfile << endl;
      cerr << "INFO:  Consider using -z if bamfile was split by region." << endl;
      exit(1);
    }

    uint max = sequences.size() ;

    if(globalOpts.lastSeqid > 0){
      max = globalOpts.lastSeqid;
    }

    int randomChr = 0;
    bool exclude = true;

    while(exclude){
      if(sequences.size() > 1){
	int prand = rand() % (max -1);
	if(globalOpts.toSkip.find(sequences[prand].RefName) == globalOpts.toSkip.end()){
	  randomChr = prand;
	  exclude = false;
	}
	if(globalOpts.toInclude.size() > 0
	   && globalOpts.toInclude.find(sequences[prand].RefName) == globalOpts.toInclude.end()){
	  exclude = true;
	}
      }
      else{
	exclude = false;
      }
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
      if(al.GetTag(  globalOpts.saT, any)){
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

  omp_set_lock(&lock);

 insertDists.mus[  targetfile ] = mu;
 insertDists.sds[  targetfile ] = sd;
 insertDists.avgD[ targetfile ] = mud;

 stringstream whereTo;

 whereTo << "INFO: for file:" << targetfile << endl
      << "      " << targetfile << ": mean depth: ......... " << mud << endl
      << "      " << targetfile << ": sd depth: ........... " << sdd << endl
      << "      " << targetfile << ": mean insert length: . " << insertDists.mus[targetfile] << endl
      << "      " << targetfile << ": median insert length. " << median                      << endl
      << "      " << targetfile << ": sd insert length .... " << insertDists.sds[targetfile] << endl
      << "      " << targetfile << ": lower insert length . " << insertDists.mus[targetfile] - (2.5*insertDists.sds[targetfile])   << endl
      << "      " << targetfile << ": upper insert length . " << insertDists.mus[targetfile] + (2.5*insertDists.sds[targetfile])   << endl
      << "      " << targetfile << ": average base quality: " << double(qsum)/double(qnum) << endl
      << "      " << targetfile << ": number of reads used: " << n  << endl << endl;



 if(globalOpts.statsOnly){
   cout << whereTo.str();
 }
 else{
   cerr << whereTo.str();
 }


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
	  bp->posCIL = -10;
	  bp->posCIH =  10;
	  bp->endCIL = -10;
	  bp->endCIH =  10;

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

void allStats(void){

#pragma omp parallel for schedule(dynamic, 1)
  for(unsigned int i = 0; i < globalOpts.targetBams.size(); i++){
    gatherBamStats(globalOpts.targetBams[i]);
  }
  if(globalOpts.statsOnly){
    cerr << "INFO: Exiting as -s flag is set." << endl;
    cerr << "INFO: WHAM finished normally, goodbye! " << endl;
    exit(0);
  }
}

void loadReads(std::vector<RefData> & sequences){
  // load bam has openMP inside for running regions quickly
  if(globalOpts.svs.empty()){
    cerr << "INFO: Loading discordant reads into forest." << endl;

    for(vector<string>::iterator bam = globalOpts.targetBams.begin();
	bam != globalOpts.targetBams.end(); bam++){

      cerr << "INFO: Reading: " << *bam << endl;

      loadBam(*bam);

      for(map<int, map<int, node * > >::iterator seqid
	    = globalGraph.nodes.begin();
	  seqid != globalGraph.nodes.end(); seqid++){

	cerr << "INFO: Number of putative breakpoints for: "
	     << sequences[seqid->first].RefName << ": "
	     << globalGraph.nodes[seqid->first].size() << endl;
      }
    }
    cerr << "INFO: Finished loading reads." << endl;
  }
}

void processAlleles(vector<breakpoints*> & allBreakpoints,
		    vector<RefData> & sequences){

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
      cerr << "INFO: generated " << nAlleles
	   << " alleles / " << allBreakpoints.size()  << endl;
    }
    omp_unset_lock(&glock);
  }
}

void refineBreakpoint(int nReads,
		      breakpoints * br,
		      map<string, vector< BamAlignment > > & ReadsPerPerson,
		      vector<RefData> & sequences){


  if(nReads > (MAXREADDEPTH * 3)){
    return;
  }

  double startingScore = totalAlignmentScore(ReadsPerPerson, br);
  int oldStart     = br->five ;
  int oldEnd       = br->three;
  int flag         = 0;

  breakpoints * secondary = new breakpoints;
  secondary->fail = false;

  for(int f = -2; f <= 2; f++){
    *secondary = *br;
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
      double newScore = totalAlignmentScore(ReadsPerPerson, secondary);

      if(newScore > startingScore){
	startingScore = newScore;
	br->five = secondary->five;
	br->three = secondary->three;
	br->svlen = br->three - br->five;
	genAlleles(br, globalOpts.fasta, sequences);
	flag = 1;
	br->refined = 1;
      }
    }
  }

  delete secondary;

  if(flag == 1){
    br->svlen = br->three - br->five;
  }
}

//-------------------------------    MAIN     --------------------------------
/*
 Comments:
*/

int main( int argc, char** argv)
{
  globalOpts.nthreads   = 1    ;
  globalOpts.lastSeqid  = 0    ;
  globalOpts.MQ         = 20   ;
  globalOpts.saT        = "SA" ;
  globalOpts.keepTrying = false;
  globalOpts.statsOnly  = false;
  globalOpts.skipGeno   = false;
  globalOpts.vcf        = true ;

  int parse = parseOpts(argc, argv);
  if(parse != 1){
    cerr << "FATAL: unable to parse command line correctly." << endl;
     printHelp();
    exit(1);
  }

  omp_set_num_threads(globalOpts.nthreads);

  if(globalOpts.fasta.empty()){
    cerr << "FATAL: no reference fasta provided." << endl << endl;
    printHelp();
    exit(1);
  }
  // gather stats for each bam (depth, insert l, ...)
  allStats();

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

 int s = 0;

 for(vector<RefData>::iterator it = sequences.begin();
     it != sequences.end(); it++){
   inverse_lookup[(*it).RefName] = s;
   s+= 1;
 }

 loadReads(sequences);

 vector<breakpoints*> allBreakpoints ;
 vector<vector<node*> > globalTrees  ;

 if(globalOpts.svs.empty()){
   cerr << "INFO: Finding trees within forest." << endl;

   gatherTrees(globalTrees);

   map<int, map <int, node*> > delBreak;

   cerr << "INFO: Finding breakpoints in trees." << endl;
#pragma omp parallel for schedule(dynamic, 3)
   for(unsigned int i = 0 ; i < globalTrees.size(); i++){

     if((i % 100000) == 0){
       omp_set_lock(&glock);
       cerr << "INFO: Processed " << i
	    << "/" << globalTrees.size() << " trees" << endl;
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

 cerr << "INFO: Trying to merge deletion breakpoints: "
      << delBreak.size() << endl;

 mergeDels(delBreak, allBreakpoints);
 }

 if(! globalOpts.svs.empty()){
   cerr << "INFO: loading external SV calls" << endl;
   loadExternal(allBreakpoints);
 }

 processAlleles(allBreakpoints, sequences);

 if(globalOpts.skipGeno){

   printVCF(allBreakpoints, sequences);
   cerr << "INFO: Skipping genotyping: -k set" << endl;
   cerr << "INFO: WHAM finished normally, goodbye! " << endl;
   return 0;
 }

 int count = 0;

#pragma omp parallel for schedule(dynamic, 3)
 for(unsigned int z = 0; z < allBreakpoints.size(); z++){

   if(allBreakpoints[z]->fail) continue;

   if((count % 100) == 0){

     omp_set_lock(&glock);
     cerr << "INFO: Refined and genotyped " << count
	  << "/" << allBreakpoints.size() << " breakpoints" << endl;
     omp_unset_lock(&glock);

   }

   omp_set_lock(&lock);
   count += 1;
   omp_unset_lock(&lock);

   map<string, vector< BamAlignment > > ReadsPerPerson;

   int nReads = getPopAlignments(globalOpts.targetBams,
				 allBreakpoints[z],
				 ReadsPerPerson, 5);

   refineBreakpoint(nReads, allBreakpoints[z],
		    ReadsPerPerson, sequences);

   for(unsigned int p = 0; p < globalOpts.targetBams.size(); p++){
     genotype(ReadsPerPerson[ globalOpts.targetBams[p] ],
	      globalOpts.targetBams[p], allBreakpoints[z]);
   }
 }

 printVCF(allBreakpoints, sequences);

 if(!globalOpts.graphOut.empty()){
   dump(globalTrees);
 }

 for(vector<breakpoints*>::iterator bks = allBreakpoints.begin();
     bks != allBreakpoints.end(); bks++){
   delete (*bks);
 }

 cerr << "INFO: WHAM finished normally, goodbye! " << endl;
 return 0;
}
