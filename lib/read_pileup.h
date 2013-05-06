//
//  read_pileup.h
//  wham
//
//  Created by Zev Kronenberg on 1/29/13.
//  Copyright (c) 2013 Zev Kronenberg. All rights reserved.
//

#ifndef wham_read_pileup_h
#define wham_read_pileup_h

#include  "api/api_global.h"
#include  "api/BamReader.h"

#include <list>

using namespace std;
using namespace BamTools;

class read_pileup {
    
 protected:
    
  int  CurrentId;
  int  CurrentPos;
  int  CurrentStart;
  std::list<BamTools::BamAlignment> Current_data;

 public:
  
  read_pileup() : CurrentId(-1) , CurrentPos(1) {};
    void proccess_alignment(BamTools::BamAlignment Current_alignment);
    void purge_all(void);
    void purge_past(void);
    std::list<BamAlignment> pileup(void);
    int  currentPos(void);
    int  currentStart(void);
   
};


#endif

