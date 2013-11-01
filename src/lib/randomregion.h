//
//  randomregion.h
//  wham
//
//  Created by Zev Kronenberg on 2/21/13.
//  Copyright (c) 2013 Zev Kronenberg. All rights reserved.
//

#ifndef __wham__randomregion__
#define __wham__randomregion__

#include  "api/api_global.h"
#include  "api/BamReader.h"
#include  "api/SamHeader.h"
#include  "api/SamSequenceDictionary.h"

class randomRegion{
 private:
  BamTools::BamReader * reader;
  BamTools::RefVector ids;
  int regionSize;
  int numberRefs;
 public:
  void init(std::string filename, int region_size);
  BamTools::BamReader * getRandom(void);
   
};

#endif /* defined(__wham__randomregion__) */
