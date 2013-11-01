//
//  flag.cpp
//  wham
//
//  Created by Zev Kronenberg on 2/19/13.
//  Copyright (c) 2013 Zev Kronenberg. All rights reserved.
//

#include "flag.h"

#define NONE                    0x0
#define MULTIPLE_SEGMENTS       0x1   // paired end ie two reads
#define PAIR_ALIGNER_BLESS      0x2   // both pairs are aligned - insert size ect via aligner
#define UNMAPPED                0x4   // read is not mapped
#define PARTNER_UNMAPPED        0x8   // pair is not mapped
#define REVERSE_STRAND          0x10
#define PARTNER_REVERSE_STRAND  0x20
#define FIRST                   0x40
#define LAST                    0x80
#define SECONDARY_ALIGNMENT     0x100
#define FAIL_QUALITY            0x200
#define DUPLICATE               0x400
#define BOTH_UNMAPPED           0xD
#define BOTH_REVERSE            0x30

void flag::addFlag(int flag){
  bamflag = flag;
}


bool flag::isPaired(void){
  if((bamflag & MULTIPLE_SEGMENTS) != 0){
    return true;
  }
  else{
    return false;
  }
}


bool flag::isPairAlignmentPass(void){
  if((bamflag & PAIR_ALIGNER_BLESS) != 0){
    return true;
  }
  else{
    return false;
  }
    
}

bool flag::isUnMapped(void){
  if((bamflag & UNMAPPED) != 0){
    return true;
  }
  else{
    return false;
  }
    
}

bool flag::isPairMapped(void){
  if((bamflag & PARTNER_UNMAPPED) == 0 ){
    return true;
  }
  else{
    return false;
  }
}


bool flag::bothUnmapped(void){
  if((bamflag & BOTH_UNMAPPED) != 0){
    return true;
  }
  else{
    return false;
  }
}

// RL
bool flag::sameStrand(void){
  if((bamflag & BOTH_REVERSE) !=0){
    return true;
  }
  else{
    return false;
  }
}

int flag::returnFlag(void){
  return bamflag;
}


