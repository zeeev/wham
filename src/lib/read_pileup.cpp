//
//  read_pileup.cpp
//  wham
//
//  Created by Zev Kronenberg on 1/29/13.
//  Copyright (c) 2013 Zev Kronenberg. All rights reserved.
//

#include "read_pileup.h"

using namespace std;
using namespace BamTools;

void read_pileup::proccess_alignment(BamTools::BamAlignment Current_alignment){
  if(Current_alignment.RefID > CurrentId){
    CurrentPos = 1;
    CurrentId  = Current_alignment.RefID;
    purge_all();
  }
  Current_data.push_back(Current_alignment);
  CurrentStart = Current_alignment.Position;
  // if(start > CurrentPos){
  //    pileup();
  // }
}
void read_pileup::purge_all(void){
  Current_data.clear();
}

void read_pileup::purge_past(void){
  bool inside = false;
  while (!inside) {
    BamAlignment read = Current_data.front();
    if (read.GetEndPosition() < CurrentPos) {
      Current_data.pop_front();
    }
    else {
      inside = true;
    }
  }
}

std::list<BamAlignment> read_pileup::pileup(void){
  read_pileup::purge_past();
  BamTools::BamAlignment last_read_in_pile =  Current_data.back();
  CurrentPos = last_read_in_pile.Position;
  return Current_data;
}

int read_pileup::currentPos(void){
  return CurrentPos;
}

int read_pileup::currentStart(void){
  return CurrentStart;
}
