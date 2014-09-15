//
//  read_pileup.h
//  wham
//
//  Created by Zev Kronenberg on 1/29/13.
//  Copyright (c) 2013 Zev Kronenberg. All rights reserved.
//

#ifndef readPileUp_h
#define readPileUp_h

#include  "api/api_global.h"
#include  "api/BamReader.h"

#include <list>



class readPileUp {
    
 public:

  int  CurrentId;
  long int  CurrentPos;
  long int  CurrentStart;
  std::list <BamTools::BamAlignment> currentData;

  readPileUp() : CurrentId(0) , CurrentPos(0) {};
  void processAlignment( BamTools::BamAlignment Current_alignment, long int pos);
  void purgeAll(void);
  void purgePast(void);
  bool softClipAtEnds(void);
  std::list<BamTools::BamAlignment> pileup(void);
  int  currentPos(void);
  int  currentStart(void);
  ~readPileUp();
};


#endif

