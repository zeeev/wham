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

void readPileUp::purgePast(void){
  while (!currentData.empty()) {
    BamAlignment read = currentData.front();
    if (read.GetEndPosition() > CurrentPos) {
      break;
    }
    else {
      currentData.pop_front();
    }
  }
}

std::list<BamAlignment> readPileUp::pileup(void){
  readPileUp::purgePast();
  BamTools::BamAlignment last_read_in_pile =  currentData.back();
  return currentData;
}

int readPileUp::currentPos(void){
  return CurrentPos;
}

int readPileUp::currentStart(void){
  return CurrentStart;
}
