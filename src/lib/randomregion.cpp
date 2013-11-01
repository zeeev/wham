//
//  randomregion.cpp
//  wham
//
//  Created by Zev Kronenberg on 2/21/13.
//  Copyright (c) 2013 Zev Kronenberg. All rights reserved.
//

#include "randomregion.h"

void randomRegion::init(std::string filename, int regions_size){
  srand(time(NULL));
  reader = new BamTools::BamReader();
  reader->Open(filename);
  numberRefs = reader->GetReferenceCount();
  ids = reader->GetReferenceData();
  regionSize = regions_size;
}

BamTools::BamReader * randomRegion::getRandom(void){
  int random_ref = rand() % numberRefs;
  int seq_len    = ids.at(random_ref).RefLength;
  int random_pos = rand() % seq_len;
  while(1){
    if((random_pos - regionSize > 0) && (random_pos + regionSize < seq_len)){
      break;
    }
    else{
      random_ref = rand() % numberRefs;
      random_pos = rand() % seq_len;
    }
  }
    
  BamTools::BamRegion random_region;

  random_region.LeftRefID     = random_ref;
  random_region.RightRefID    = random_ref;
  random_region.LeftPosition  = random_pos;
  random_region.RightPosition = (random_pos + regionSize);
    
  //std::cout << random_ref << "\t" << random_pos << "\t" << regionSize << "\n";
  reader->SetRegion(random_region);
  return reader;

}
