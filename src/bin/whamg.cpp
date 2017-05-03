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

// openMP - swing that hammer
#include <omp.h>

// bamtools and my headers
#include "api/BamMultiReader.h"
#include "readPileUp.h"

// phred scaling
#include "phredUtils.h"

using namespace std;
using namespace BamTools;

typedef map< char, map< int, map< int, breakpoint * > > >  unpairedSVs;
typedef map< int, map< int, breakpoint * > > unpairedSVsChr;
typedef map< int, breakpoint * > unpairedSVsTypePos;


struct options{
    std::vector<string> targetBams;
    bool statsOnly                ;
    bool skipGeno                 ;
    bool keepTrying               ;
    bool noInterSeqid             ;
    int MQ                        ;
    int NM                        ;
    int nthreads                  ;
    int lastSeqid                 ;
    int minPairMatch              ;
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
//int->seqid index
map<int, string> forward_lookup;

// options

static const char *optString = "c:i:u:m:r:a:g:x:f:e:d:hsz";

// omp lock
omp_lock_t lock;
// omp lock for the graph
omp_lock_t glock;
// omp lock for the pairstor
omp_lock_t pslock;

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
  cerr << " Basic usage:  " << endl;
  cerr << "       whamg -f my.bam -a my.fasta \\ " << endl;
  cerr << "             -e M,GL000207.1 2> wham.err > wham.vcf"
       << endl;
  cerr << endl;
  cerr << " Required:  " << endl;
//------------------------------- XXXXXXXXXX --------------------------------

  cerr << "          -f  <STRING>  Comma separated list of bam files or a file with    " << endl;
  cerr << "                        one bam (full path) per line.                       " << endl;
  cerr << "          -a  <STRING>  The reference genome (indexed fasta).               " << endl;
  cerr << endl;
  cerr << " Optional:  Recommended flags are noted with : *                            " << endl;
  cerr << "          -s  <FLAG>    Exits the program after the stats are               " << endl;
  cerr << "                        gathered. [false]                                   " << endl;
  cerr << "          -g  <STRING>  File to write graph to (very large output). [false] " << endl;
  cerr << "  *|-c    -e  <STRING>  Comma sep. list of seqids to skip [false].          " << endl;
  cerr << "  *|-e    -c  <STRING>  Comma sep. list of seqids to keep [false].          " << endl;
  cerr << "          -r  <STRING>  Region in format: seqid:start-end [whole genome]    " << endl;
  cerr << "  *       -x  <INT>     Number of CPUs to use [1 CPU].                      " << endl;
  cerr << "          -m  <INT>     Mapping quality filter [20].                        " << endl;
  cerr << "          -i  <STRING>  non standard split read tag [SA]                    " << endl;
  cerr << "          -z  <FLAG>    Sample reads until success. [false]                 " << endl;
  cerr << "          -d  <INT>     Minimum number of matching bases (both reads).[100] " << endl;

  cerr << endl;
  cerr << " Output:  " << endl;
  cerr << "        STDERR: Run statistics and bam stats                                " << endl;
  cerr << "        STOUT : SV calls in VCF                                             " << endl;
  cerr << endl;
  cerr << " Details:  " << endl;
  cerr << "        -z  <FLAG>    WHAM-GRAPHENING can fail if does not sample           " << endl;
  cerr << "                      enough reads. This flag prevents whamg                " << endl;
  cerr << "                      from exiting. If your bam header has seqids not in    " << endl;
  cerr << "                      the bam (e.g. split by region) use -z.                " << endl;
  cerr << "        -i  <STRING>  WHAM-GRAPHENING uses the optional bwa-mem SA tag.     " << endl;
  cerr << "                      Older version of bwa-mem used XP.                     " << endl;
  cerr << "     -e|-c  <STRING>  A list of seqids to include or exclude while          " << endl;
  cerr << "                      sampling insert and depth.  For humans you should     " << endl;
  cerr << "                      use the standard chromosomes 1,2,3...X,Y.             " << endl;
  cerr << endl;

  printVersion();
}

//------------------------------- SUBROUTINE --------------------------------
bool breakSort(breakpoint * L, breakpoint * R){

    if(L->IsMasked() || R->IsMasked()){
        return false;
    }

    if(L->nodeL->seqid == R->nodeL->seqid){
        if(L->nodeL->pos <= R->nodeL->pos){
            return true;
        }
        else{
            return false;
        }
    }
    else{
        if(L->nodeL->seqid < R->nodeL->seqid){
            return true;
        }
        else{
            return false;
        }
    }
    return false;
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : NA
 Function does   : prints vcf header
 Function returns: NA
*/
void printVCF(std::vector<breakpoint*> & bp){

    std::vector<breakpoint*> scrubbed;

    for(std::vector<breakpoint*>::iterator it = bp.begin();
        it != bp.end(); it++){
        if((*it)->IsMasked()){
            continue;
        }
        if(!(*it)->IsPrint()){
            continue;
        }
        if((*it)->nodeL->seqid != (*it)->nodeR->seqid){
            continue;
        }

        if((*it)->getType() == 'T'){
            continue;
        }
        if((*it)->getTotalSupport() < 3){
            continue;
        }
        scrubbed.push_back(*it);
    }


    sort(scrubbed.begin(), scrubbed.end(), breakSort);

    stringstream header;

    header << "##fileformat=VCFv4.2"      << std::endl;
    header << "##source=WHAM-GRAPHENING:" << VERSION << std::endl;
    header << "##reference=" << globalOpts.fasta << std::endl;
    header << "##INFO=<ID=A,Number=1,Type=Integer,Description=\"Total pieces of evidence\">" << std::endl;
    header << "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">" << std::endl;
    header << "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">" << std::endl;
    header << "##INFO=<ID=CF,Number=1,Type=Float,Description=\"Fraction of reads in graph that cluster with SVTYPE pattern\">" << std::endl;
    header << "##INFO=<ID=CW,Number=5,Type=Float,Description=\"SVTYPE weight 0-1; DEL,DUP,INV,INS,BND\">" << std::endl;
    header << "##INFO=<ID=D,Number=1,Type=Integer,Description=\"Number of reads supporting a deletion\">" << std::endl;
    header << "##INFO=<ID=DI,Number=1,Type=Float,Description=\"Average distance of mates to breakpoint\">" << std::endl;
    header << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">" << std::endl;
    header << "##INFO=<ID=EV,Number=1,Type=Integer,Description=\"Number everted mate-pairs\">" << std::endl;
    header << "##INFO=<ID=I,Number=1,Type=Integer,Description=\"Number of reads supporting an insertion\">" << endl;
    header << "##INFO=<ID=SR,Number=1,Type=Integer,Description=\"Number of split-reads supporing SV\">" << std::endl;
    header << "##INFO=<ID=SS,Number=1,Type=Integer,Description=\"Number of split-reads supporing SV\">" << std::endl;
    header << "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">" << std::endl;
    header << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">" << std::endl;
    header << "##INFO=<ID=T,Number=1,Type=Integer,Description=\"Number of reads supporting a BND\">" << std::endl;
    header << "##INFO=<ID=TAGS,Number=.,Type=String,Description=\"SM tags with breakpoint support\">" << std::endl;
    header << "##INFO=<ID=TF,Number=1,Type=Integer,Description=\"Number of reads mapped too far\">" << std::endl;
    header << "##INFO=<ID=U,Number=1,Type=Integer,Description=\"Number of reads supporting a duplication\">" << std::endl;
    header << "##INFO=<ID=V,Number=1,Type=Integer,Description=\"Number of reads supporting an inversion\">" << std::endl;
    header << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
    header << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">" << std::endl;
    header << "##FORMAT=<ID=SP,Number=1,Type=Integer,Description=\"Per sample SV support\">" << std::endl;
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

    for(std::vector<breakpoint*>::iterator it = scrubbed.begin();
        it != scrubbed.end(); it++){
        if((*it)->IsMasked()){
            continue;
        }
        if((*it)->getType() == 'T'){
            continue;
        }
        else{
            (*it)->loadSMSupport();
            std::cout << **it << "\tGT:DP:SP";
            for(vector<string>::iterator iz = globalOpts.targetBams.begin();
                iz !=  globalOpts.targetBams.end(); iz++){
                std::cout << "\t.:.:" << (*it)->getSMSupport(globalOpts.SMTAGS[*iz]);
            }
            std::cout << std::endl;
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

  map<edge *, int>  seenEdges ;
  map<node *, int>  seenNodes ;
  vector<edge *>        edges ;

  edges.insert(edges.end(), n->eds.begin(), n->eds.end());

  /* if something is pushed to the back of the vector it changes the
     positions ! be warned. */

  while(!edges.empty()){

      edge * last = edges.back();
      edges.pop_back();

      for(std::vector<edge *>::iterator it = last->L->eds.begin();
          it != last->L->eds.end(); it++){

          if(seenEdges.find(*it) == seenEdges.end()){
              edges.push_back(*it);
              seenEdges[*it]      = 1;
              seenNodes[(*it)->L] = 1;
              seenNodes[(*it)->R] = 1;
          }
      }

      for(std::vector<edge *>::iterator it = last->R->eds.begin();
          it != last->R->eds.end(); it++){

          if(seenEdges.find(*it) == seenEdges.end()){
              edges.push_back(*it);
              seenEdges[*it]      = 1;
              seenNodes[(*it)->L] = 1;
              seenNodes[(*it)->R] = 1;
          }
      }

  }
  for(std::map<node *, int>::iterator it = seenNodes.begin();
      it != seenNodes.end(); it++){
      ns.push_back(it->first);
  }

}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : read pair pointer
 Function does   : returns true is reads are everted
 Function returns: bool
*/

inline bool isEverted(readPair * rp){

    if(!rp->al1.IsMapped() || !rp->al2.IsMapped()){
        return false;
    }
    if(rp->al1.RefID != rp->al2.RefID){
        return false;
    }

    if(rp->al1.Position <= rp->al2.Position){
        if(rp->al1.IsReverseStrand()
           &&  (!rp->al2.IsReverseStrand()) ){
            return true;
        }
    }
    else{
        if(rp->al2.IsReverseStrand()
           &&  (!rp->al1.IsReverseStrand()) ){
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
      == globalGraph.nodes[refIDL][l]->sm.end()){
     globalGraph.nodes[refIDL][l]->sm[SM] = 1;
   }
   else{
     globalGraph.nodes[refIDL][l]->sm[SM] += 1;
   }
   if(globalGraph.nodes[refIDR][r]->sm.find(SM)
      == globalGraph.nodes[refIDR][r]->sm.end()){
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
 Function input  : bam alignment, sm tag
 Function does   : finds the positions of indels
 Function returns: string
*/

bool indelToGraph(BamAlignment & ba, string & SM){

  if(ba.MapQuality < 30){
    return false;
  }

  bool hit = false;

  int p = ba.Position;

  for(vector<CigarOp>::iterator ci = ba.CigarData.begin();
      ci != ba.CigarData.end(); ci++){

      // 20bp seems reasonable
      if(ci->Length < 20){
          continue;
      }

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

    if(! rp->al1.IsMapped() || ! rp->al2.IsMapped() ){
        return true;
    }

    if(rp->al1.MapQuality < 5
       || rp->al2.MapQuality < 5 ){
        return true;
    }

    if(rp->al1.MapQuality < globalOpts.MQ
       && rp->al2.MapQuality < globalOpts.MQ){
        return true;
    }
    if(rp->al1.Length == rp->al1.CigarData[0].Length
       && rp->al1.CigarData[0].Type == 'M' &&
       rp->al2.Length == rp->al2.CigarData[0].Length
       && rp->al2.CigarData[0].Type == 'M' ){
        return true;
    }
    if((match(rp->al1.CigarData) + match(rp->al2.CigarData)) < globalOpts.minPairMatch){
      return true;
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

  // too many chimeras ... for me
  if(sa.size() > 1){
    return;
  }
  // too many mismatches
  if(sa.front().NM > globalOpts.NM){
      return;
  }
  // map quality of zero...
  if(sa.front().mapQ < 5){
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

  bool readFront = true;

  if(al.CigarData.back().Type == 'S'){
    start = al.GetEndPosition(false,true);
    readFront = false;
  }

  /* we also know that both sides of the split read are not trimmed */

  if(sa.front().cig.back().Type == 'S'){
    endPos(sa[0].cig, &end)  ;
  }

  if(readFront){
      if(end > start){
          support = 'Z';
      }
  }
  else{
      if(end < start){
          support = 'Z';
      }
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

  //can't be deviant on different seqids
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
        end = rp->al2.GetEndPosition(false,true);
    }
    if(rp->al1.CigarData.back().Type == 'S'){
      start = rp->al1.GetEndPosition(false,true);
    }
  }
  else{
    if(rp->al1.CigarData.back().Type == 'S'){
        end = rp->al1.GetEndPosition(false,true);
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
 Function input  : pointer to readPair
 Function does   : calulates if overlap is too much
 Function returns: bool
*/

inline bool tooMuchOverlap(readPair * rp){

    if(rp->al1.RefID != rp->al2.RefID){
        return false;
    }

    if(rp->al1.Position > rp->al2.GetEndPosition(false,true) ||
       rp->al2.Position > rp->al1.GetEndPosition(false,true)){
        return false;
    }

    long int maxStart = std::max(rp->al1.Position, rp->al2.Position);
    long int minEnd   = std::min(rp->al1.GetEndPosition(false,true),
                                 rp->al1.GetEndPosition(false,true));

    double perOver = double(minEnd - maxStart)
        / double(rp->al1.Length + rp->al2.Length);

    if(perOver > 0.2){
        return true;
    }


    return false;
}



//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : pointer to readPair ; seqid to index
 Function does   : processes Pair
 Function returns: NA
*/

void processPair(readPair * rp,
         double * low,
         double * high,
         string & SM){

  if(globalOpts.noInterSeqid && rp->al1.RefID != rp->al2.RefID){
      return;
  }

  if(pairFailed(rp)){
    return;
  }
  if(isPointIn(rp)){
      return;
  }
  if(tooMuchOverlap(rp)){
      return;
  }

  string sa1;
  string sa2;

  bool sameStrand = false;
  bool everted    = isEverted(rp);

  string xa;
  std::vector<std::string> xas;

  if(rp->al1.GetTag("XA", xa)){
      xas = split(xa, ";");
      if(xas.size() > 3){
          return;
      }
  }
  if(rp->al2.GetTag("XA", xa)){
      xas = split(xa, ";");
      if(xas.size() > 3){
          return;
      }
  }

  if(rp->al1.RefID == rp->al2.RefID){
      indelToGraph(rp->al1, SM);
      indelToGraph(rp->al2, SM);
  }

  // nm filter
  std::string nm;
  if(rp->al1.GetTag("NM", nm)){
      if(atoi(nm.c_str()) > globalOpts.NM){
          return;
      }
  }

  if(rp->al2.GetTag("NM", nm)){
      if(atoi(nm.c_str()) > globalOpts.NM){
          return;
      }
  }

  if( ! IsLongClip(rp->al1.CigarData, 5)
      && ! IsLongClip(rp->al2.CigarData, 5)){
      return;
  }

  if((rp->al1.IsReverseStrand() && rp->al2.IsReverseStrand())
     || (! rp->al1.IsReverseStrand() && ! rp->al2.IsReverseStrand()) ){
    sameStrand = true;
  }

  if(! everted){
      if( abs(rp->al1.InsertSize) > *high){
          if(sameStrand){
              deviantInsertSize(rp, 'M', SM);
          }
          else{
              deviantInsertSize(rp, 'H', SM);
          }
      }
      else if(abs(rp->al1.InsertSize) < *low ){
          if(sameStrand){
              deviantInsertSize(rp, 'R', SM);
          }
          else{
              deviantInsertSize(rp, 'L', SM);
          }
      }
      else{
          if(sameStrand){
              deviantInsertSize(rp, 'A', SM);
          }
      }
  }
  else{
      deviantInsertSize(rp, 'X', SM);
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

  // local read pair store
  map<string, readPair *>pairStoreLocal;

  double high = insertDists.upr[filename] ;
  double low  = insertDists.low[filename] ;

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
    if((al.AlignmentFlag & 0x0100) != 0 ){
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
  omp_set_lock(&pslock);

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
  omp_unset_lock(&pslock);
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
          if((Mb % 100) == 0 ){
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
  omp_set_lock(&pslock);
  for(map<string, readPair*>::iterator rps = globalPairStore.begin();
      rps != globalPairStore.end(); rps++){
      delete rps->second;
  }
  globalPairStore.clear();
  omp_unset_lock(&pslock);
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : string
 Function does   : loads filenames for bam files
 Function returns: bool
*/

bool loadTargetBamFileNames(std::string fns){

    std::string line;

    std::vector<std::string> ext = split(fns, ".");
    if(ext.back() == "bam"){
            globalOpts.targetBams  = split(fns, ",");
    }
    else{
        ifstream myfile(fns.c_str());

        if(myfile){
            while (getline( myfile, line )){
                globalOpts.targetBams.push_back(line);
            }
        }
    }
    return true;
}
//-------------------------------   OPTIONS   --------------------------------
int parseOpts(int argc, char** argv)
{
  int opt = 0;
  opt = getopt(argc, argv, optString);
  while(opt != -1){
    switch(opt){
    case 'd':
        {
            globalOpts.minPairMatch = atoi(((string)optarg).c_str());
            cerr << "INFO: whamg will keep read pairs with " << globalOpts.minPairMatch << " matching bases." << endl;
            break;
        }
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
            cerr << "INFO: WHAM will skip seqid: " << optarg << endl;
            for(unsigned int i = 0; i < seqidsToSkip.size(); i++){
                globalOpts.toSkip[seqidsToSkip[i]] = 1;
            }
            break;
        }
    case 'c':
        {
            vector<string> seqidsToInclude = split(optarg, ",");
            for(unsigned int i = 0; i < seqidsToInclude.size(); i++){
                globalOpts.toInclude[seqidsToInclude[i]] = 1;
                cerr << "INFO: WHAM will analyze seqid: " << seqidsToInclude[i] << endl;
            }
            break;
        }

    case 'f':
        {
            loadTargetBamFileNames(string(optarg));
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
    default:
        {
            cerr << "FATAL: unknown command line option" << endl;
            exit(1);
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
          if((*iz)->support['A'] > 0){
              if((*it)->pos != (*iz)->L->pos){
                  ss << "     " << (*it)->seqid << "." << (*it)->pos << " -- " << (*iz)->L->seqid << "." << (*iz)->L->pos << " [c\
olor=orange,penwidth=" << (*iz)->support['A'] << "];\n";
              }
              if((*it)->pos != (*iz)->R->pos){
                  ss << "     " << (*it)->seqid << "." << (*it)->pos << " -- " << (*iz)->R->seqid << "." << (*iz)->R->pos << " [c\
olor=orange,penwidth=" << (*iz)->support['A'] << "];\n";
              }
          }
      }
  }

  ss << "}";

  return ss.str();
}
//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector<nodes *>, and a breakpoint
 Function does   : checks if inversion meets basic requirements
 Function returns: NA
*/

void doubleCheckIns(std::vector<breakpoint *> & bks, unpairedSVs & up){

    for(std::vector<breakpoint *>::iterator it = bks.begin();
        it != bks.end(); it++){

        if((*it)->IsMasked()){
            continue;
        }
        if((*it)->getType() != 'I'){
            continue;
        }
        if(((*it)->getTooCloseCount() > 2 && (*it)->getInsCount() > 2)
           ||  (*it)->getInsCount()  > 4){

            // switching starts because one is a terminal node
            if((*it)->nodeR->eds.size() > (*it)->nodeL->eds.size()){
                node * tmp;
                tmp = (*it)->nodeR;
                (*it)->nodeR = (*it)->nodeL;
                (*it)->nodeL = tmp;
            }


        }
        else{
            (*it)->unSetPrint();

        }
    }
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector<nodes *>, and a breakpoint
 Function does   : checks if inversion meets basic requirements
 Function returns: NA
*/

void doubleCheckInv(std::vector<breakpoint *> & bks, unpairedSVs & up){

    for(std::vector<breakpoint *>::iterator it = bks.begin();
        it != bks.end(); it++){

        if((*it)->IsMasked()){
            continue;
        }
        if((*it)->getType() != 'V'){
            continue;
        }
        if((*it)->getClustFrac() > 0.1
           && (*it)->getSameStrandCount() > 2
           && (*it)->getInvCount() > 2 ){
        }
        else{
            if((*it)->getLength()  < 500
               && ((*it)->getSameStrandCount() < 3
                   || (*it)->getInvCount() < 3
                   || (*it)->getSplitReadCount() < 1)
               ){
                (*it)->unSetPrint();
            }
        }
    }
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector<nodes *>, and a breakpoint
 Function does   : finds the best within graph link
 Function returns: string
*/

void doubleCheckDup(std::vector<breakpoint *> & bks, unpairedSVs & up){

    for(std::vector<breakpoint *>::iterator it = bks.begin();
        it != bks.end(); it++){

        if((*it)->IsMasked()){
            continue;
        }
        if((*it)->getType() != 'U'){
            continue;
        }
        if(((*it)->getClustFrac()   > 0.1 &&
           (*it)->getEvertCount()  > 1
            && (*it)->getDupCount() > 1)
           || ((*it)->getSplitReadCount() > 3 && (*it)->getLength() < 500)
           ){
        }
        else{
            (*it)->unSetPrint();
        }
    }
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector<nodes *>, and a breakpoint
 Function does   : finds the best within graph link
 Function returns: string
*/

void doubleCheckDel(std::vector<breakpoint *> & bks, unpairedSVs & up){

    for(std::vector<breakpoint *>::iterator it = bks.begin();
        it != bks.end(); it++){

        if((*it)->IsMasked()){
            continue;
        }
        if((*it)->getType() != 'D'){
            continue;
        }
        if(((*it)->getClustFrac() > 0.1
           && (*it)->getTooFarCount() > 1
           && (*it)->getDelCount() > 1 )
           || ((*it)->getInternalDelCount() > 2 && (*it)->getSplitReadCount() > 2)
           ){
        }
        else{
            (*it)->unSetPrint();

            up[(*it)->getType()][(*it)->nodeL->seqid][(*it)->nodeL->pos] = *it;
            up[(*it)->getType()][(*it)->nodeR->seqid][(*it)->nodeR->pos] = *it;
        }
    }
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : unpairedSVs & up, type char
 Function does   : tries to cluster and link deletions
 Function returns: NA
*/
void clusterUnlinked(unpairedSVs & up, char type){

    for(unpairedSVsChr::iterator it = up[type].begin();
        it != up[type].end(); it++ ){

        for(unpairedSVsTypePos::iterator iz = it->second.begin();
            iz != it->second.end(); iz++){

            unpairedSVsTypePos::iterator forwardScan = iz;
            forwardScan++;
            for(; forwardScan != it->second.end(); forwardScan++){


                if(iz->second->delClusterCheck(iz->second->nodeL,
                                       forwardScan->second->nodeR)){

                    std::cerr << "R hit " << iz->second->nodeL->pos
                              << " " << forwardScan->second->nodeR->pos
                              << " " << iz->second->getClustFrac()
                              << " " << iz->second->getNClustered()
                              << " " << iz->second->getAvgDist() << std::endl;

                }
                if(iz->second->delClusterCheck(iz->second->nodeL,
                                               forwardScan->second->nodeL)){

                    std::cerr << "L hit " << iz->second->nodeL->pos
                              << " " << forwardScan->second->nodeL->pos
                              << " " << iz->second->getClustFrac()
                              << " " << iz->second->getNClustered()
                              << " " << iz->second->getAvgDist() << std::endl;


                }

            }
        }
    }
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector<nodes *>, and a breakpoint
 Function does   : finds the best within graph link
 Function returns: string
*/

void findPairs(vector<node*> & tree,
               breakpoint * bp,
               map<int, map<int, breakpoint* > > & lookup){

    if(tree.size() < 3 || tree.size() > 200){
        bp->setMasked();
        return;
    }

    int lmax = 0;
    int rmax = 0;
    int tmax = 0;
    int edsn = 0;

    node * finalL = NULL;
    node * finalR = NULL;

    for(vector<node *>::iterator lit = tree.begin();
        lit != tree.end(); lit++){

      for(vector<node *>::iterator rit = tree.begin();
          rit != tree.end(); rit++){

          if(*lit == *rit){
              continue;
          }

          if(!connectedNode(*lit, *rit)){
              continue;
          }
          int lmaxTmp = getSupport(*lit);
          int rmaxTmp = getSupport(*rit);
          int tmaxTmp = lmaxTmp + rmaxTmp;
          int edsnTmp = ((*lit)->eds.size() + (*rit)->eds.size());

          if(tmaxTmp > tmax){
              finalL  = *lit;
              finalR  = *rit;
              lmax    = lmaxTmp;
              rmax    = rmaxTmp;
              tmax    = tmaxTmp;
              edsn    = edsnTmp;
          }
      }
    }

    bp->add(finalL);
    bp->add(finalR);

    if(finalL->eds.size() == 1 || finalR->eds.size() == 1){

        if(finalL->eds.front()->support['D'] < 5 &&
           finalL->eds.front()->support['I'] < 5 &&
           finalL->eds.front()->support['S'] < 5 &&
           finalL->eds.front()->support['V'] < 5 &&
           finalR->eds.front()->support['D'] < 5 &&
           finalR->eds.front()->support['I'] < 5 &&
           finalR->eds.front()->support['S'] < 5 &&
           finalR->eds.front()->support['V'] < 5   ){
            bp->setBadPair();
        }
    }
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : processes trees
 Function does   : tries to define type and breakpoints
 Function returns: NA
*/

void gatherTrees(vector<vector<node *> > & trees){

  map<int, map<int, int> > lookup;

  for(map<int, map<int, node* > >::iterator it = globalGraph.nodes.begin();
      it != globalGraph.nodes.end(); it++){
    for(map<int, node*>::iterator itt = it->second.begin();
        itt != it->second.end(); itt++){

      if(lookup[it->first].find(itt->first) != lookup[it->first].end() ){
      }
      else{
          lookup[it->first][itt->first] = 1;
          vector<node *> tree;
          getTree(globalGraph.nodes[it->first][itt->first], tree);
          for(vector<node *>::iterator ir = tree.begin(); ir != tree.end(); ir++){
              lookup[(*ir)->seqid][(*ir)->pos] = 1;
          }
          trees.push_back(tree);
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
 Function input  : bam file
 Function does   : generates stats for bam file
 Function returns: NA
*/

void gatherBamStats(string & targetfile){

  omp_set_lock(&lock);
  int quals [126];
  memcpy(quals, SangerLookup, 126*sizeof(int));

  cerr << "INFO: gathering stats (may take some time) for bam: " << targetfile << endl;

  omp_unset_lock(&lock);

  vector<double> alIns        ;
  vector<double> nReads       ;
  vector<double> randomSWScore;
  vector<int>    mapQuality   ;

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
          mapQuality.push_back(al.MapQuality);
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
  double muQ      = mean(mapQuality);

  omp_set_lock(&lock);

 insertDists.mus[  targetfile ] = mu;
 insertDists.sds[  targetfile ] = sd;
 insertDists.avgD[ targetfile ] = mud;

 insertDists.low[ targetfile ] = insertDists.mus[targetfile]
   - (1.6*insertDists.sds[targetfile]);

 if(insertDists.low[ targetfile ] < 0 ){
   insertDists.low[ targetfile ] = 0;
 }

 insertDists.upr[ targetfile ] = insertDists.mus[targetfile]
   + (2.5*insertDists.sds[targetfile]);


 stringstream whereTo;

 whereTo << "INFO: for Sample:" << SM << endl
         << "  STATS:    " << SM << ": mean depth ..............: " << mud
         << std::endl
         << "  STATS:    " << SM << ": sd depth ................: " << sdd
         << std::endl
         << "  STATS:    " << SM << ": mean insert length: .....: " << insertDists.mus[targetfile]
         << std::endl
         << "  STATS:    " << SM << ": median insert length ....: " << median
         << std::endl
         << "  STATS:    " << SM << ": sd insert length ........: " << insertDists.sds[targetfile]
         << std::endl
         << "  STATS:    " << SM << ": lower insert cutoff .....: " << insertDists.low[targetfile]
         << std::endl
         << "  STATS:    " << SM << ": upper insert cutoff .....: " << insertDists.upr[targetfile]
         << std::endl
         << "  STATS:    " << SM << ": average base quality ....: " <<  double(qsum)/double(qnum)
         << std::endl
         << "  STATS:    " << SM << ": average mapping quality .: " << muQ
         << std::endl
         << "  STATS:    " << SM << ": number of reads used ....: " << n
         << std::endl << std::endl;

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

//-------------------------------    MAIN     --------------------------------
/*
 Comments:
*/

int main( int argc, char** argv)
{
  globalOpts.nthreads     = 1    ;
  globalOpts.lastSeqid    = 0    ;
  globalOpts.MQ           = 20   ;
  globalOpts.NM           = 10   ;
  globalOpts.minPairMatch = 100  ;
  globalOpts.saT          = "SA" ;
  globalOpts.keepTrying   = false;
  globalOpts.statsOnly    = false;
  globalOpts.skipGeno     = false;
  globalOpts.noInterSeqid = true ;

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
      forward_lookup[s] = (*it).RefName;
      s+= 1;
  }

  loadReads(sequences);

  vector<breakpoint*>                allBreakpoints;
  vector<vector<node*> >                globalTrees;
  map<int, map<int, breakpoint*> > breakpointLookup;

  if(globalOpts.svs.empty()){
      cerr << "INFO: Gathering graphs from forest." << endl;

      gatherTrees(globalTrees);

      for(unsigned int i = 0 ; i < globalTrees.size(); i++){
          breakpoint * tmp = new breakpoint;
          allBreakpoints.push_back(tmp);
      }

      cerr << "INFO: Matching breakpoints." << endl;
#pragma omp parallel for schedule(dynamic, 3)
      for(unsigned int i = 0 ; i < globalTrees.size(); i++){
          if((i % 10000) == 0){
              omp_set_lock(&lock);
              cerr << "INFO: Processed " << i
                   << "/" << globalTrees.size() << " graphs" << endl;
              omp_unset_lock(&lock);
          }
          findPairs(globalTrees[i], allBreakpoints[i], breakpointLookup);

          if(allBreakpoints[i]->IsMasked() ){
              continue;
          }
          allBreakpoints[i]->countSupportType();
          allBreakpoints[i]->calcType();
          allBreakpoints[i]->delClusterCheck();
          allBreakpoints[i]->invClusterCheck();
          allBreakpoints[i]->dupClusterCheck();
          allBreakpoints[i]->getRefBases(globalOpts.fasta, forward_lookup);
      }
      cerr << "INFO: Printing." << endl;

      unpairedSVs unMatched;

      doubleCheckDel(allBreakpoints, unMatched );
      doubleCheckDup(allBreakpoints, unMatched );
      doubleCheckInv(allBreakpoints, unMatched );
      doubleCheckIns(allBreakpoints, unMatched );


      //      std::cerr << "INFO: clustering unlinked deletions." << std::endl;
      //      clusterUnlinked(unMatched, 'D');


      printVCF(allBreakpoints);

      std::cerr << "INFO: done processing trees" << std::endl;

      if(!globalOpts.graphOut.empty()){
          dump(globalTrees);
      }
  }

  if(!globalOpts.graphOut.empty()){
      dump(globalTrees);
  }

  for(vector<breakpoint*>::iterator bks = allBreakpoints.begin();
     bks != allBreakpoints.end(); bks++){
      delete (*bks);
  }

 cerr << "INFO: WHAM finished normally, goodbye! " << endl;
 return 0;
}
