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


using namespace std;
using namespace BamTools;

struct options{
  std::vector<string> targetBams;
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
  string strand;
  vector<cigDat> cig;
};

struct node;
struct edge;

struct edge{
  node * L;
  node * R;
  int forwardSupport;
  int reverseSupport;
  map<char,int> support;

};

struct node{
  int   seqid;
  int    pos;
  vector <edge *> eds  ;
};


struct graph{
  map< int, map<int, node *> > nodes;
  vector<edge *>   edges;
}globalGraph;


// options

static const char *optString = "f:h";

// omp lock

omp_lock_t lock;
omp_lock_t glock;

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
  
  //  seen[n->pos] = 1;

  //  ns.push_back(n);

  // if something is pushed to the back of the stack it changes the positions ! be warned.

  while(!edges.empty()){
    
     uint hit = 0;
     
     if(seen.find(edges.back()->L->pos) != seen.end() && seen.find(edges.back()->R->pos) != seen.end() ){
       hit = 1;
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
 Function input  : edge pointer 

 Function does   : init

 Function returns: void

*/

void initEdge(edge * e){
  e->L = NULL;
  e->R = NULL;
  e->forwardSupport  = 0;
  e->reverseSupport  = 0;

  e->support['S'] = 0;
  e->support['I'] = 0;
  e->support['D'] = 0;

}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of strings and separator

 Function does   : joins vector with separator

 Function returns: string

*/

string join(vector<string> strings, string sep){

  string joined = "";

  for(vector<string>::iterator sit = strings.begin(); sit != strings.end(); sit++){
    joined = joined + sep + (*sit) ;
  }
  return joined;
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of strings

 Function does   : joins vector with returns;

 Function returns: string

*/

string joinReturn(vector<string> strings){

  string joined = "";

  for(vector<string>::iterator sit = strings.begin(); sit != strings.end(); sit++){
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
    
    if(rp->al1.CigarData.back().Type == 'S' && rp->al2.CigarData.front().Type == 'S'){
      return true;
    }
  }
  else{
    if(rp->al1.CigarData.front().Type == 'S' && rp->al2.CigarData.back().Type == 'S'){
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

void addIndelToGraph(int refID, int l, int r, char s){

  omp_set_lock(&glock);

  if( ! isInGraph(refID, l, globalGraph) &&  ! isInGraph(refID, r, globalGraph) ){
    
    node * nodeL;
    node * nodeR;
    edge * ed   ;

    nodeL = new node;
    nodeR = new node;
    ed    = new edge;

    initEdge(ed);

    ed->support[s] +=1;

    ed->forwardSupport += 1;
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
 else if(isInGraph(refID, l, globalGraph) &&  ! isInGraph(refID, r, globalGraph)){
   node * nodeR;
   edge * ed;

   nodeR = new node;
   ed    = new edge;

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
 else if(! isInGraph(refID, l, globalGraph) &&  isInGraph(refID, r, globalGraph)){
   node * nodeL;
   edge * ed;

   nodeL = new node;
   ed    = new edge;

   initEdge(ed);
   ed->support[s] +=1;
   nodeL->pos   = l;
   nodeL->seqid = refID;
   ed->R =  globalGraph.nodes[refID][r];
   ed->L = nodeL;
   
   nodeL->eds.push_back(ed);

   globalGraph.nodes[refID][r]->eds.push_back(ed);
   globalGraph.nodes[refID][l] = nodeL;  
}
 else{
   uint hit = 0;

   for(vector<edge *>::iterator ite = globalGraph.nodes[refID][l]->eds.begin();
       ite != globalGraph.nodes[refID][l]->eds.end(); ite++){
     if((*ite)->L->pos == l && (*ite)->R->pos == r){

       (*ite)->support[s] += 1;

       (*ite)->forwardSupport += 1;
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

bool indelToGraph(BamAlignment & ba){

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
    case 'X':
      {
	p += ci->Length;
	break;
      }
    case 'I':
      {
	hit = true;

	addIndelToGraph(ba.RefID, p, p + ci->Length, 'I');
	
	break;
      }
    case 'D':
      {
	hit = true;
	addIndelToGraph(ba.RefID, p, p + ci->Length, 'D');

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

  for(vector<CigarOp>::iterator sit = strings.begin(); sit != strings.end(); sit++){
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

bool IsLongClip(vector<CigarOp> & ci, int len){

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
    if(rp->al1.Length == rp->al1.CigarData[0].Length && rp->al1.CigarData[0].Type == 'M' &&
       rp->al2.Length == rp->al2.CigarData[0].Length && rp->al2.CigarData[0].Type == 'M' ){
      return true;
    }
    if(rp->al1.MapQuality < 20 && rp->al2.MapQuality < 20){
      return true;
    }
    if((match(rp->al1.CigarData) + match(rp->al2.CigarData)) < 100){
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

void parseCigar(vector<cigDat> parsedCigar, string cigar){

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

    sDat.seqid = il[sat[0]];
    sDat.pos   = atoi(sat[1].c_str()) - 1;
    sDat.strand = sat[2];
    parseCigar(sDat.cig,sat[3]);
    parsed.push_back(sDat);

  }
  
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : bam alignment, vector<saTag> 

 Function does   : finds links for splitter to put in graph

 Function returns: NA

*/


void splitToGraph(BamAlignment al, vector<saTag> & sa){

  if(!al.IsMapped()){
    cerr << "fail: not mapped" << endl;
    return;
  }

  if(sa.size() > 1){
    cerr << "fail: too many fragments in split" << endl;
    return;
  }

  if(sa[0].seqid != al.RefID){
    cerr << "fail: split to different chr" << endl;
    return;
  }

  if(al.CigarData[0].Type == 'S' && al.CigarData[al.CigarData.size() -1 ].Type == 'S'){
    return;
  }
  else if(al.CigarData[0].Type == 'S'){
    int start = al.Position; 
    int end   = sa[0].pos  ;
    
    if(start > end){
      start = sa[0].pos;
      end   = al.Position;
    }
    
    addIndelToGraph(al.RefID, start, end, 'S');
  }
  else{
    int start =  al.GetEndPosition(false,true);
    int end   = sa[0].pos                  ;

    if(start > end){
      start = sa[0].pos;
      end   = al.GetEndPosition(false,true);
    }
    addIndelToGraph(al.RefID, start, end, 'S');
  }
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : pointer to readPair ; seqid to index

 Function does   : processes Pair

 Function returns: NA

*/

void processPair(readPair * rp, map<string, int> & il){
  
  string sa1;
  string sa2;

  if(pairFailed(rp)){
    return;
  }

  if(rp->al1.RefID != rp->al2.RefID){
    return;
  }

  if(rp->al1.IsMapped() && rp->al2.IsMapped()){
    if(! IsLongClip(rp->al1.CigarData, 5) && ! IsLongClip(rp->al2.CigarData, 5)){
      return;
    }
  }

  if(rp->al1.GetTag("SA", sa1)){
    vector<saTag> parsedSa1;
    parseSA(parsedSa1, sa1, il);
    //    cerr << sa1 << endl;
    splitToGraph(rp->al1, parsedSa1);
    //    cerr << "s1 processed " << endl;
  }
  if(rp->al2.GetTag("SA", sa2)){
    vector<saTag> parsedSa2;
    parseSA(parsedSa2, sa2, il);
    //    cerr << sa2 << endl;
        splitToGraph(rp->al1, parsedSa2);
    //    cerr << "s2 processed " << endl;
  }

  indelToGraph(rp->al1);
  indelToGraph(rp->al2);

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
	       map<string, int> & seqInverseLookup){


  map<string, int> localInverseLookUp;



  omp_set_lock(&lock);
  for(map<string, int>::iterator itt = seqInverseLookup.begin();
      itt !=  seqInverseLookup.end(); itt++){    
    localInverseLookUp[itt->first] = itt->second;
  }
  omp_unset_lock(&lock);

  // local graph;
  graph localGraph;


  // local read pair store

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
      processPair(pairStore[al.Name], localInverseLookUp);
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
      processPair(globalPairStore[rps->first], localInverseLookUp);
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

  // if the bam is not sorted die
  if(!SH.HasSortOrder()){
    cerr << "FATAL: sorted bams must have the @HD SO: tag in each SAM header: " << bamFile  << endl;
    exit(1);
  }

  RefVector sequences = br.GetReferenceData();

  //chunking up the genome
  int seqidIndex = 0;
  vector< regionDat* > regions;

  // inverse lookup for split reads 
  map<string,int> seqIndexLookup;
  

  for(vector< RefData >::iterator sit = sequences.begin(); sit != sequences.end(); sit++){
    int start = 0;
    if(seqidIndex != 0){
      seqidIndex+=1;
      continue;
    }
    seqIndexLookup[sequences[seqidIndex].RefName] = seqidIndex;

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

  // closing the bam reader before running regions
  br.Close();

  // global read pair store

  map<string, readPair*> pairStore;


  // running the regions with openMP
#pragma omp parallel for schedule(dynamic, 3)
  
  for(unsigned int re = 0; re < regions.size(); re++){
    
    omp_set_lock(&lock);
    cerr << "INFO: running region: " 
	 << sequences[regions[re]->seqidIndex].RefName 
	 << ":" << regions[re]->start 
	 << "-" << regions[re]->end << endl;
    omp_unset_lock(&lock);
    
    if(! runRegion(bamFile, 
		   regions[re]->seqidIndex, 
		   regions[re]->start, 
		   regions[re]->end, 
		   sequences,
		   pairStore,
		   seqIndexLookup)){
      omp_set_lock(&lock);
      cerr << "WARNING: region failed to run properly: "
           << sequences[regions[re]->seqidIndex].RefName
           << ":"  << regions[re]->start << "-"
           << regions[re]->end
           <<  endl;
      omp_unset_lock(&lock);
    }
  }


  cerr << "INFO: " << bamFile << " had " << pairStore.size() << " reads that were not processed" << endl; 
}
//-------------------------------   OPTIONS   --------------------------------
int parseOpts(int argc, char** argv)
    {
    int opt = 0;
    opt = getopt(argc, argv, optString);
    while(opt != -1){
      switch(opt){
      case 'f':
	{
	  globalOpts.targetBams     = split(optarg, ",");
	  cerr << "INFO: target bams:\n" << joinReturn(globalOpts.targetBams) ;
	  break;
	}
      case 'h':
	{
	  cerr << "nada" << endl;
	  break;
	}
      case '?':
	{
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


string dotviz(vector<node *> ns){

  stringstream ss;

  ss << "graph {\n";

  for(vector<node *>::iterator it = ns.begin(); it != ns.end();
      it++){
    for(vector<edge *>:: iterator iz = (*it)->eds.begin(); 
	iz != (*it)->eds.end(); iz++){

      if((*iz)->support['S'] > 0){
	if((*it)->pos != (*iz)->L->pos){
	  ss << "     " << (*it)->seqid << "." << (*it)->pos << " -- " << (*iz)->L->seqid << "." << (*iz)->L->pos << " [style=dotted] ;\n" ;	
	}
	if((*it)->pos != (*iz)->R->pos){
	  ss << "     " << (*it)->seqid << "." << (*it)->pos << " -- " << (*iz)->R->seqid << "." << (*iz)->R->pos << " [style=dotted] ;\n" ;
	}
      }
      else if((*iz)->support['I'] > 0){
	if((*it)->pos != (*iz)->L->pos){
	  ss << "     " << (*it)->seqid << "." << (*it)->pos <<  " -- " << (*iz)->L->seqid << "." << (*iz)->L->pos << " [color=red] ;\n" ;
	}
	if((*it)->pos != (*iz)->R->pos){
	  ss << "     " << (*it)->seqid << "." << (*it)->pos << " -- " << (*iz)->R->seqid << "." << (*iz)->R->pos << " [color=red] ;\n" ;
	}
      }
      else if((*iz)->support['D'] > 0){
	if((*it)->pos != (*iz)->L->pos){
	  ss << "     " << (*it)->seqid << "." << (*it)->pos << " -- " << (*iz)->L->seqid << "." << (*iz)->L->pos << " [color=blue] ;\n" ;
	}
	if((*it)->pos != (*iz)->R->pos){
	  ss << "     " << (*it)->seqid << "." << (*it)->pos << " -- " << (*iz)->R->seqid << "." << (*iz)->R->pos << " [color=blue] ;\n" ;
	}
      }
      else{
	if((*it)->pos != (*iz)->L->pos){
	  ss << "     " << (*it)->seqid << "." << (*it)->pos << " -- " << (*iz)->L->seqid << "." << (*iz)->L->pos << ";\n" ;
	}
	if((*it)->pos != (*iz)->R->pos){
	  ss << "     " << (*it)->seqid << "." << (*it)->pos << " -- " << (*iz)->R->seqid << "." << (*iz)->R->pos << ";\n" ;
	}
      }
    }
  }

  ss << "}";


  return ss.str();
}




//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : nothing

 Function does   : dumps and shrinks graph

 Function returns: NA

*/


void dump(){
  
  map<int, map<int, int> > lookup;
  
  for(map<int, map<int, node* > >::iterator it = globalGraph.nodes.begin();it != globalGraph.nodes.end(); it++){
      for(map<int, node*>::iterator itt = it->second.begin(); itt != it->second.end(); itt++){
	
	if(lookup[it->first].find(itt->first) != lookup[it->first].end() ){
	}
	else{
	  lookup[it->first][itt->first] = 1;
	  
	  vector<node *> tree;
	  

	  getTree(globalGraph.nodes[it->first][itt->first], tree);
	  
	  if(tree.size() < 20 || tree.size() > 100){
	    continue;
	  }


	  string dotvizg = dotviz(tree);
	  cout << "Tree: " << endl;
	  for(vector<node *>::iterator ir = tree.begin(); ir != tree.end(); ir++){
	    cout << "NODE: " << (*ir)->pos << " " << (*ir)->seqid << endl;
	    lookup[(*ir)->seqid][(*ir)->pos] = 1;

	    for(vector<edge *>::iterator iz = (*ir)->eds.begin(); iz != (*ir)->eds.end(); iz++){
	      cout << " EDGE: L: " << (*iz)->L->pos << " R: " <<  (*iz)->R->pos << endl;
	      cout << " SUPPORT: " << (*iz)->forwardSupport << " I:" << (*iz)->support['I'] << " D:" << (*iz)->support['D'] << " S:" << (*iz)->support['S'] << endl; 

	    }
	  }
	  cout << endl;
	  cout << endl << dotvizg << endl;
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
int parse = parseOpts(argc, argv);


 for(vector<string>::iterator bam = globalOpts.targetBams.begin();
     bam != globalOpts.targetBams.end(); bam++){

   loadBam(*bam);

   cerr << "INFO: Graph now has " << globalGraph.nodes.size() << " nodes." << endl;
 }

 dump();

 return 0;

}
