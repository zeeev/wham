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
   
  void addFlag(int flag);
  bool isPaired(void);
  bool isPairAlignmentPass(void);
  bool isUnMapped(void);
  bool isPairMapped(void);
  bool bothUnmapped(void);
  bool sameStrand(void);
  int  returnFlag(void);
};
