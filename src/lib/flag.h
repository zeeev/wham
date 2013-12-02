//
//  flag.h
//  wham
//
//  Created by Zev Kronenberg on 2/19/13.
//  Copyright (c) 2013 Zev Kronenberg. All rights reserved.
//

#ifndef __wham__flag__
#define __wham__flag__

#endif /* defined(__wham__flag__) */



class flag{
 private:
  int bamflag;
 public:

  ///addFlag takes a SAM specified flag 
   
  void addFlag(int flag);

  ///isPaired determines if the added flag is paired

  bool isPaired(void);

  /// isPairAlignmentPass determines if the added flag passes the mapper's "proper pair" - very subjective

  bool isPairAlignmentPass(void);

  /// isUnMapped determines if the added flag is mapped 

  bool isUnMapped(void);

  /// isPairMapped determines if the added flag's mate is mapped

  bool isPairMapped(void);

  /// bothUnmapped determines if both pair1 and pair2 are unmapped

  bool bothUnmapped(void);

  /// sameStrand determines if both pair1 and pair2 are on the same strand

  bool sameStrand(void);

  /// bothForStrand determines if both pair1 and pair2 are on the reverse strand

  bool bothRevStrand(void);

  /// bothForStrand determines if both pair1 and pair2 are on the forward strand

  bool bothForStrand(void);

  /// returnFlag returns the added flag

  int  returnFlag(void);
};
