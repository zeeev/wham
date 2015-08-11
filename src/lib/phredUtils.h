//
//  phredUtils.h
//  alignHMM
//
//  Created by Zev Kronenberg on 7/29/15.
//  Copyright (c) 2015 Zev Kronenberg. All rights reserved.
//

#ifndef alignHMM_phredUtils_h
#define alignHMM_phredUtils_h

#include <cmath>
#include "string.h"
#include "float.h"

const int MAXQUAL = 60;
const int MAXQUAL_LOG10 = -6;

const int SangerLookup[127] = {-1,-1,-1,-1,-1, -1,-1,-1,-1,-1, // 0-9     1-10
			       -1,-1,-1,-1,-1, -1,-1,-1,-1,-1, // 10-19   11-20
			       -1,-1,-1,-1,-1, -1,-1,-1,-1,-1, // 20-29   21-30
			       -1,-1,-1, 0, 1,  2, 3, 4, 5, 6, // 30-39   31-40
			        7, 8, 9, 10,11, 12,13,14,15,16, // 40-49   41-50
			        17,18,19,20,21, 22,23,24,25,26, // 50-59   51-60
			        27,28,29,30,31, 32,33,34,35,36, // 60-69   61-70
			        37,38,39,40,41, 42,43,44,45,46, // 70-79   71-80
			        47,48,49,50,51, 52,53,54,55,56, // 80-89   81-90
			        57,58,59,60,61, 62,63,64,65,66, // 90-99   91-100
			        67,68,69,70,71, 72,73,74,75,76, // 100-109 101-110
			        77,78,79,80,81, 82,83,84,85,86, // 110-119 111-120
			        87,88,89,90,91, 92,93           }; // 120-119 121-130

const int IlluminaOneThree[127] = {-1,-1,-1,-1,-1, -1,-1,-1,-1,-1, // 0-9     1-10
			           -1,-1,-1,-1,-1, -1,-1,-1,-1,-1, // 10-19   11-20
			           -1,-1,-1,-1,-1, -1,-1,-1,-1,-1, // 20-29   21-30
			           -1,-1,-1,-1,-1, -1,-1,-1,-1,-1, // 30-39   31-40
			           -1,-1,-1,-1,-1  -1,-1,-1,-1,-1, // 40-49   41-50
			           -1,-1,-1,-1,-1, -1,-1,-1,-1,-1, // 50-59   51-60
			           -1,-1,-1,-1,-1,  0, 1, 2, 3, 4, // 60-69   61-70
                                    5, 6, 7, 8, 9, 10,11,12,13,14, // 70-79   71-80
			            15,16,17,18,19, 20,21,22,23,24, // 80-89   81-90
			            25,26,27,28,29, 30,31,32,33,34, // 90-99   91-100
			            35,36,37,38,39, 40,-1,-1,-1,-1, // 100-109 101-110
			            -1,-1,-1,-1,-1, -1,-1,-1,-1,-1, // 110-119 111-120
			            -1,-1,-1,-1,-1, -1,-1        }; // 120-119 121-130

int LookUp[126];

class phredUtils{

    
private:

    
    
    bool IsFiniteNumber(double x)
    {
        return (x <= DBL_MAX && x >= -DBL_MAX);
    }
    
    
public:
    phredUtils(void){
                memcpy(LookUp, SangerLookup, 127*sizeof(int));
    }
    void setIllumina(void){
              memcpy(LookUp, IlluminaOneThree, 127*sizeof(int));
    }
    int qualToPhred(char b){
        return LookUp[int(b)];
    }
    double qualToProb(char b){
        int  v = LookUp[int(b)];
        return pow(10.0,(-1*double(v)/10));
    }
    double qualToProbLog10(char b){
        return log10(1-qualToProb(b));
    }
    double qualToProbErrorLog10(char b){
         return log10(qualToProb(b));
    }
    double phredToLog10(int i){
        double p = pow(10.0,(-1*double(i)/10));
        return log10(p);
    }
    double phredToProb(int i){
        return pow(10.0,(-1*double(i)/10));
    }
    double log10Add(double a, double b){
        if(a > b){
            return log10Add(b, a);
        }
        if(! IsFiniteNumber(a)){
            return b;
        }
        double diff = b - a;

        if(diff >= 8){
            return b;
        }
	
        return  b + log10(1.0 + pow(10, -1*diff));
    }
 
};

#endif
