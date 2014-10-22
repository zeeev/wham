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

readPileUp::~readPileUp(){}

void readPileUp::processAlignment(BamTools::BamAlignment Current_alignment, long int pos){
  currentData.push_back(Current_alignment);
  CurrentStart    = Current_alignment.Position;
  CurrentPos      = pos;
}

void readPileUp::purgeAll(void){
  currentData.clear();
}

bool readPileUp::softClipAtEnds(void){
  
  bool clipped = false;
  
  for( list< BamAlignment >::iterator it = currentData.begin(); it != currentData.end(); it++ ){

    vector< CigarOp > cd = (*it).CigarData;
    if(cd.back().Type == 'S' || cd.front().Type == 'S' || cd.back().Type == 'H' ||  cd.front().Type == 'H'){
      clipped = true;
      break;
    }  
  }
  return clipped;
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
