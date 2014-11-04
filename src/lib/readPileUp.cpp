//
//  readPileUp.cpp
//  wham
//
//  Created by Zev Kronenberg on 1/29/13.
//  Copyright (c) 2013 Zev Kronenberg. All rights reserved.
//

#include "readPileUp.h"

using namespace std;
using namespace BamTools;

bool sameStrand(BamAlignment & al){

  if(( al.IsReverseStrand() && al.IsMateReverseStrand() ) 
     || ( ! al.IsReverseStrand() && ! al.IsMateReverseStrand() )){
    return true;
  }
  return false;
}

bool readPileUp::processDiscordant(BamAlignment & al, string & saTag){

  nDiscordant++;

  if(sameStrand(al)){
    nsameStrandDiscordant++;
  }

  if((al.RefID =! al.MateRefID)){
    ndiscordantCrossChr++;
  }

  odd[al.Name]++;
  
  clusterFrontOrBackPrimary(al, true, saTag);

  return true;

}

bool readPileUp::processSplitRead(BamAlignment & al, string & saTag){

  vector<string> sas = split(saTag, ";");

  if(sas.size() > 2){
    return true;
  }

  odd[al.Name]++;

  nsplitRead += 1;

  vector<string> saData = split(sas[0], ",");
  
  // checking if the two splitread fragments
  // are on the same strand

  if(saData[2].compare("+") == 0){
    if(!al.IsReverseStrand() ){
      nf1f2SameStrand += 1;
    }
  }
  else{
    if(al.IsReverseStrand() ){
      nf1f2SameStrand++;
    }
  }
  
  // checking fragment 1 and fragment 2 
  // vs the mate pair

  if(al.IsMateMapped()){
    
    // checking the first fragment 
    // against the mate pair

    if(sameStrand(al)){
      nf1SameStrand++;
    }

    // checking the second fragment
    // against the mate pair

    if(saData[2].compare("+") == 0){
      if(!al.IsReverseStrand()){
	nf2SameStrand++;
      }
    }
    else{
      nsplitMissingMates++;
    }
    
    // splitread translocation 

    if(al.RefID != al.MateRefID){
      nsplitReadCrossChr++;
    }
  }

  clusterFrontOrBackPrimary(al, false, saTag);
  
  return true;
  
}


bool readPileUp::processMissingMate(BamAlignment & al, string & saTag){

  nMatesMissing++;

  clusterFrontOrBackPrimary(al, true, saTag);

  odd[al.Name]++;

  return true;

}

bool readPileUp::processProperPair(BamAlignment & al, string & saTag){

  nPaired++ ;

  if(sameStrand(al)){
    nSameStrand += 1;
    odd[al.Name]++;
  }
  if(al.RefID != al.MateRefID){
    nCrossChr++;
    odd[al.Name]++;
  }

  clusterFrontOrBackPrimary(al, true, saTag);

  vector< CigarOp > cd = al.CigarData;
  for(vector< CigarOp >::iterator cig = cd.begin(); 
      cig != cd.end(); cig++){

    switch((*cig).Type){
    case 'I':
      {
	if((*cig).Length > 25){
	  internalInsertion += 1;
	  odd[al.Name]++;
	}
	break;
      }
    case 'D':
      {
	if((*cig).Length > 25){
	  internalDeletion += 1;
	  odd[al.Name]++;
	}
	break;
      }
    default:
      {
      }
    }
  }
  
  return true;
}


bool readPileUp::clusterFrontOrBackPrimary(BamAlignment & al, bool p, string & saTag){

  vector< CigarOp > cd = al.CigarData;  

  if((al.AlignmentFlag & 0x0800) != 0){
    if(cd.front().Type == 'H'){
      allCount[al.Position]++;
      supplementCount[al.Position]++;
      supplement[al.Position].push_back(al);
    }
    if(cd.back().Type == 'H'){
      allCount[al.Position]++;
      supplementCount[al.Position]++;
      supplement[al.GetEndPosition()].push_back(al);
    }
  }

  else{
    if(cd.front().Type == 'S'){
      nClippedFront++;
      allCount[al.Position]++;
      primaryCount[al.Position]++;
      primary[al.Position].push_back(al);
      if(! saTag.empty()){
	supplement[al.Position].push_back(al);
      }
    }
    if(cd.back().Type == 'S'){
      nClippedBack++;
      allCount[al.GetEndPosition()]++;
      primaryCount[al.GetEndPosition()]++;
      primary[al.GetEndPosition()].push_back(al);
      if(! saTag.empty()){
	supplement[al.GetEndPosition()].push_back(al);
      }
    }
  }
  return true;
}

void readPileUp::printPileUp(void){
  for(list<BamAlignment>::iterator r = currentData.begin();
      r != currentData.end(); r++){
    cerr << (*r).Name 
	 << "\t"
	 << (*r).Position
	 << "\t"
	 << (*r).QueryBases
	 << endl;
  }

}

void readPileUp::processPileup(long int * pos){

  clearClusters();
  clearStats();

  for(list<BamAlignment>::iterator r = currentData.begin(); 
      r != currentData.end(); r++){
  
    // trailing pileup data
    if((*r).Position > *pos){
      continue;
#ifdef DEBUG
      cerr << "Too far ahead: " << (*r).Name << endl;
#endif
    }

    numberOfReads += 1;

    string saTag;

    // split reads
    if( (*r).GetTag("SA", saTag) ){
      processSplitRead(*r, saTag);
#ifdef DEBUG
      cerr << "Split Read: " << (*r).Name << endl;
#endif
      continue;
    }
    // discordant reads
    if(!(*r).IsProperPair()){
      processDiscordant(*r, saTag);
#ifdef DEBUG
      cerr << "Discordant Read: " << (*r).Name << endl;
#endif
      continue;
    }
    // mates missing
    if(!(*r).IsMateMapped()){
      processMissingMate(*r, saTag);
#ifdef DEBUG
      cerr << "Mate Missing: " << (*r).Name << endl;
#endif
      continue;
    }
    // good data
    if((*r).IsMateMapped() && (*r).IsProperPair() ){
      
#ifdef DEBUG
      cerr << "before count: " << nPaired << endl;
#endif

      processProperPair(*r, saTag);
#ifdef DEBUG
      cerr << "Mate paired: " << (*r).Name << endl;
      cerr << "after count: " << nPaired << endl;
#endif
      continue;
    }    
  
#ifdef DEBUG
    cerr << "Bleed through: " << (*r).Name << endl;
#endif
    
  }
}

void readPileUp::clearClusters(void){
  odd.clear();
  primaryCount.clear();
  supplementCount.clear();
  primary.clear();
  supplement.clear();
  allCount.clear();
}

void readPileUp::clearStats(void){
  internalInsertion = 0;
  internalDeletion  = 0;
  numberOfReads = 0;
  nMatesMissing = 0;
  nPaired = 0;
  nCrossChr = 0;
  nsplitRead = 0;
  nsplitReadCrossChr = 0;
  nf1SameStrand = 0;
  nf2SameStrand = 0;
  nf1f2SameStrand = 0;
  nsplitMissingMates = 0;
  nDiscordant = 0;
  ndiscordantCrossChr = 0;
  nsameStrandDiscordant = 0;
  nSameStrand = 0;
  nDiscordant = 0;
  nSoftClipped = 0;
  nClippedFront = 0; 
  nClippedBack = 0;
}

readPileUp::readPileUp(){
  CurrentPos   = 0;
  CurrentStart = 0;
}

readPileUp::~readPileUp(){}

void readPileUp::processAlignment(BamTools::BamAlignment Current_alignment, long int pos){
  currentData.push_back(Current_alignment);
  CurrentStart    = Current_alignment.Position;
  CurrentPos      = pos;
}

void readPileUp::purgeAll(void){
  currentData.clear();
}


void readPileUp::purgePast(void){
  
  int nreads  = currentData.size();
  int counter = 0;
  
  while (!currentData.empty()) {
    
    BamAlignment read = currentData.front();
    
    if( (read.Position + read.Length) < CurrentStart){ 
      currentData.pop_front();
    } 
   
    counter += 1;
    
    if(counter == nreads){
      break;
    }
  }
}

std::list<BamAlignment> readPileUp::pileup(void){
  readPileUp::purgePast();
  return currentData;
}

int readPileUp::currentPos(void){
  return CurrentPos;
}

int readPileUp::currentStart(void){
  return CurrentStart;
}

int readPileUp::nReads(void){
  return currentData.size();
}
