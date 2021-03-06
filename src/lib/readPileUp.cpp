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

  if(!al.IsMateMapped()){
    nMatesMissing++;
    odd[al.Name] = 1;
    clusterFrontOrBackPrimary(al, true, saTag);
    return true;
  }

  if(al.RefID != al.MateRefID){
    nCrossChr++;
    ndiscordantCrossChr++; 
 }
  
  if(sameStrand(al)){
    nSameStrand += 1;
    odd[al.Name] = 1;
  }
  
  if(sameStrand(al)){
    nsameStrandDiscordant++;
  }

  odd[al.Name] = 1;
  
  clusterFrontOrBackPrimary(al, true, saTag);

  return true;

}

bool readPileUp::processSplitRead(BamAlignment & al, string & saTag){

  vector<string> sas = split(saTag, ";");


  odd[al.Name] = 1;
  nsplitRead  += 1;

  if(al.IsProperPair()){
    nPaired++;
  }
  else{
    nDiscordant++;
  }

  //everts
  if(al.IsReverseStrand()
     && ! al.IsMateReverseStrand()
     && al.Position < al.MatePosition){
    evert++;
  }
  else if(! al.IsReverseStrand()
	  && al.IsMateReverseStrand()
	  && al.Position > al.MatePosition){
    evert++;
  }


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

  odd[al.Name] = 1;

  return true;

}

bool readPileUp::processPair(BamAlignment & al, string & saTag){

  #ifdef DEBUG
  //  cerr << "processing pair: " << al.Name << " " << al.Position << endl;
  #endif

  nPaired++;

  if(al.IsMateMapped()){
    if(al.RefID != al.MateRefID){
      nCrossChr++;
      odd[al.Name] = 1;
    }
    else if(sameStrand(al)){
      nSameStrand += 1;
      odd[al.Name] = 1;
    }
    else if(al.IsReverseStrand() 
	    && ! al.IsMateReverseStrand() 
	    && al.Position < al.MatePosition){
      evert++;
    }
    else if(! al.IsReverseStrand()
	    && al.IsMateReverseStrand()
	    && al.Position > al.MatePosition){
      evert++;
    }
  }
  
  clusterFrontOrBackPrimary(al, true, saTag);

  vector< CigarOp > cd = al.CigarData;

  for(vector< CigarOp >::iterator cig = cd.begin(); 
      cig != cd.end(); cig++){

    switch((*cig).Type){
    case 'I':
      {
	if((*cig).Length > 3){
	  internalInsertion += 1;
	  odd[al.Name] = 1;
	}
	break;
      }
    case 'D':
      {
	if((*cig).Length > 3){
	  internalDeletion += 1;
	  odd[al.Name] = 1;
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


#ifdef DEBUG
  //  cerr << "clustering read: " << al.Name << " " << al.Position << endl;
#endif

  if((al.AlignmentFlag & 0x0800) != 0){
    if(cd.front().Type == 'H'){
      supplement[al.Position].push_back(al);
    }
    if(cd.back().Type == 'H' ){
      supplement[al.GetEndPosition(false,true)].push_back(al);
    }
  }
  else{
    if(cd.front().Type == 'S'){
   
      #ifdef DEBUG
      //      cerr << "front clip: " << al.Name << " " << al.Position << endl; 
      #endif
      
      nClippedFront++;
      primary[al.Position].push_back(al);
      odd[al.Name] = 1;
      if(! saTag.empty()){
	supplement[al.Position].push_back(al);
      }
    }
    if(cd.back().Type == 'S'){


     #ifdef DEBUG
      //      cerr << "back clip: " << al.Name << " " << al.GetEndPosition(false,true) << endl;
      #endif

      nClippedBack++;
      primary[al.GetEndPosition(false,true)].push_back(al);
      odd[al.Name] = 1;
      if(! saTag.empty()){
	supplement[al.GetEndPosition(false,true)].push_back(al);
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

#ifdef DEBUG
      //      cerr << "fail pos greater: " << *pos << " " << (*r).Name << " " << (*r).Position << endl;
#endif
      continue;
    }

    vector< CigarOp > cd = (*r).CigarData;
    if(cd.size() > 3){
      tooManyCigs += 1;
      
#ifdef DEBUG
      //      cerr << "fail too many cigars" << endl;
#endif 

    }

#ifdef DEBUG
    //    cerr << "currentDat: " << (*r).Name << " " << *pos << endl;
#endif

    mapQsum += (*r).MapQuality ;

    numberOfReads += 1;

    if((*r).MapQuality < 50){
      nLowMapQ += 1;
    }
    
    if(! (*r).IsMateMapped()){
      nMatesMissing += 1;
    }

    string saTag;   
    // split reads
    if( (*r).GetTag("SA", saTag ) ){
      processSplitRead(*r, saTag);
      continue;
    }
   
    // paired end data
    if((*r).IsPaired()){
      if(!(*r).IsProperPair()){
	nDiscordant++;
      }
      processPair(*r, saTag);
      continue;
    }
    
#ifdef DEBUG
    cerr << "Bleed through: " << (*r).Name << endl;
#endif
    
  }
}

void readPileUp::clearClusters(void){
  odd.clear();
  primary.clear();
  supplement.clear();
}

void readPileUp::clearStats(void){
  tooManyCigs           = 0;
  nLowMapQ              = 0;
  mapQsum               = 0;
  mateTooClose          = 0;
  mateTooFar            = 0;
  evert                 = 0;
  internalInsertion     = 0;
  internalDeletion      = 0;
  numberOfReads         = 0;
  nMatesMissing         = 0;
  nPaired               = 0;
  nCrossChr             = 0;
  nsplitRead            = 0;
  nsplitReadCrossChr    = 0;
  nf1SameStrand         = 0;
  nf2SameStrand         = 0;
  nf1f2SameStrand       = 0;
  nsplitMissingMates    = 0;
  nDiscordant           = 0;
  ndiscordantCrossChr   = 0;
  nsameStrandDiscordant = 0;
  nSameStrand           = 0;
  nDiscordant           = 0;
  nSoftClipped          = 0;
  nClippedFront         = 0; 
  nClippedBack          = 0;
}

readPileUp::readPileUp(){
  CurrentPos   = 0;
  CurrentStart = 0;
}

readPileUp::~readPileUp(){}

void readPileUp::processAlignment(BamTools::BamAlignment Current_alignment){
  currentData.push_back(Current_alignment);
  CurrentStart    = Current_alignment.Position;
}

void readPileUp::purgeAll(void){
  currentData.clear();
}

void readPileUp::purgePast(long int * delPos){
  
  int nreads  = currentData.size();
  int counter = 0;
  
  while (!currentData.empty()) {
    
    BamAlignment read = currentData.front();
    currentData.pop_front();    

    if( read.GetEndPosition(false,true) >= *delPos){ 
      currentData.push_back(read);
    } 
   
    counter += 1;
    
    if(counter == nreads){
      break;
    }
  }
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
