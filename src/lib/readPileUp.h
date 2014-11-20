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
#include  "split.h"

#include <list>
#include <map>
#include <vector>

class readPileUp {
    
 public:

  int  CurrentId;
  long int  CurrentPos;
  long int  CurrentStart;
  std::list <BamTools::BamAlignment> currentData;

  std::map <std::string, int> odd;
  std::map <long int, int > primaryCount, supplementCount, allCount;
  std::map <long int, std::vector<BamTools::BamAlignment> > primary, supplement;

  int numberOfReads;

  // internal shit
  int internalInsertion;
  int internalDeletion ;
  // split read info
  int nsplitRead;
  int nsplitReadCrossChr;
  int nf1SameStrand;
  int nf2SameStrand;
  int nf1f2SameStrand;
  int nsplitMissingMates;

  // discordant reads
  int nDiscordant;
  int nsameStrandDiscordant;
  int ndiscordantCrossChr;

  // paired info
  int nPaired      ;
  int nMatesMissing;
  int nSameStrand  ;
  int nCrossChr    ;

  //clustered info
  int nSoftClipped ;
  int nClippedFront;
  int nClippedBack ;
  
  readPileUp() ;
  ~readPileUp();

  bool clusterFrontOrBackPrimary(BamTools::BamAlignment &, bool, std::string&);
  bool processSplitRead(BamTools::BamAlignment&, std::string&);
  bool processDiscordant(BamTools::BamAlignment &, std::string&);
  bool processSupplement(BamTools::BamAlignment &, std::string&);
  bool processMissingMate(BamTools::BamAlignment &, std::string&);
  bool processProperPair(BamTools::BamAlignment &, std::string&);

  void processAlignment( BamTools::BamAlignment);
  void processPileup(long int *);
  void printPileUp(void);
  void purgeAll(void);
  void purgePast(long int *);
  void clearStats(void);
  void clearClusters(void);

  int  currentPos(void);
  int  currentStart(void);
  int  nReads(void);
};

#endif

