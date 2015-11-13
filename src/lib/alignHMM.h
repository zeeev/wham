//
//  matrices.h
//  alignHMM
//
//  Created by Zev Kronenberg on 7/29/15.
//  Copyright (c) 2015 Zev Kronenberg. All rights reserved.
//

#ifndef alignHMM_matrices_h
#define alignHMM_matrices_h
#include <string>
#include <iostream>
#include <cmath>
#include "phredUtils.h"
#include "float.h"

const uint MATCH_TO_MATCH_I = 0;
const uint MATCH_TO_INS_I   = 1;
const uint MATCH_TO_DEL_I   = 2;
const uint INDEL_TO_MATCH   = 3;
const uint INS_TO_INS       = 4;
const uint DEL_TO_DEL       = 5;

class alignHMM {

private:
    double *  transitions;
    double ** priorMatrix;
    double ** matchMatrix;
    double ** insertionMatrix;
    double ** deletionMatrix;

    phredUtils pu;
    
    uint NROW; /*!< Read length + 1     */
    uint NCOL; /*!< Haplotype length +1 */

    void dumpMat(double ** mat, int nrow, int ncol){
        for(int i = 0; i < nrow; i++){
            for(int j = 0; j < ncol; j++){
                std::cerr << mat[i][j] << "\t";
            }
            std::cerr << std::endl;
        }
    }
   
public:
    ~alignHMM(){
      delete transitions;
      for(uint i = 0; i < NROW; i++){
	delete[] matchMatrix[i];
	delete[] insertionMatrix[i];
	delete[] deletionMatrix[i];
	delete[] priorMatrix[i];
      }

    }

    alignHMM(int nrow, int ncol){

      if(nrow > ncol){
	std::cerr << "FATAL: haplotype smaller than read." << std::endl;
	exit(1);
      }
      
      NROW = nrow;
      NCOL = ncol;
      
      matchMatrix      = new double * [NROW];
      insertionMatrix  = new double * [NROW];
      deletionMatrix   = new double * [NROW];
      priorMatrix      = new double * [NROW];
      transitions      = new double [6];
      
      for(uint i = 0; i < NROW; i++){
	matchMatrix[i]      = new double [NCOL];
	insertionMatrix[i]  = new double [NCOL];
	deletionMatrix[i]   = new double [NCOL];
	priorMatrix[i]      = new double [NCOL];
      }
      for(uint i = 0; i < NROW; i++){
	for(uint j = 0; j < NCOL; j++){
	  matchMatrix[i][j]     = -INFINITY;
	  insertionMatrix[i][j] = -INFINITY;
	  deletionMatrix[i][j]  = -INFINITY;
	  priorMatrix[i][j]     = -INFINITY;
	}
      }
    }

    bool initPriors(std::string & haplotype,
                    std::string & readSeq,
                    std::string & readQual){



        
      for(uint i = 0; i < readSeq.size(); i++){
	char b = readSeq[i];
	char q = readQual[i];
	
	for(uint j = 0; j < haplotype.size(); j++){
	  char h = haplotype[j];
	  priorMatrix[i+1][j+1] = (b == h || b == 'N' || h == 'N') ?
	    pu.qualToProbLog10(q) : pu.qualToProbErrorLog10(q);
	}
      }
      return true;
    }

    bool initTransProbs(void){

      transitions[MATCH_TO_MATCH_I] = -2.744828e-05 ;
      transitions[MATCH_TO_INS_I]   = -4.5          ;
      transitions[MATCH_TO_DEL_I]   = -4.5          ; 
      transitions[INDEL_TO_MATCH]   = -0.04575749   ;
      transitions[INS_TO_INS]       = -1            ;
      transitions[DEL_TO_DEL]       = -1            ;
      
      return true;
        
    }
    
    void initializeDelMat(void){
      //        double initialValue = log10(1.0 / NCOL);
      double initialValue = log10(1.0 / NCOL);
        
        for(uint i = 0 ; i < NCOL; i++){
            deletionMatrix[0][i] = initialValue;
        }
    }

    // So very sorry (gross syntax)
    void updatecells(void){
        for (uint i = 1; i < NROW; i++) {
            for(uint j = 1; j < NCOL; j++){
	      
	      matchMatrix[i][j] = priorMatrix[i][j] + 
		pu.log10Add(
			    pu.log10Add(
					(matchMatrix[i - 1][j -1]     +transitions[MATCH_TO_MATCH_I]),
					(insertionMatrix[i - 1][j -1] +transitions[INDEL_TO_MATCH])
					),
			    (deletionMatrix[i - 1][j -1]  +transitions[INDEL_TO_MATCH])
			    );
	      
	                    
	      insertionMatrix[i][j] = pu.log10Add(
					       (matchMatrix[i - 1][j] + transitions[MATCH_TO_INS_I]),
					       (insertionMatrix[i - 1][j] +  transitions[INS_TO_INS])
					       );
	      
	      deletionMatrix[i][j] = pu.log10Add(
					      (matchMatrix[i][j - 1] +  transitions[MATCH_TO_DEL_I]),
					      (deletionMatrix[i][j-1] + transitions[DEL_TO_DEL]) 
					      ); 
            }
        }
    }
    
    double finalLikelihoodCalculation(){
      uint endI = NROW - 1;
      double finalProbSum = 
	pu.log10Add(
		    pu.log10Add(matchMatrix[endI][1], insertionMatrix[endI][1]),
		    pu.log10Add(matchMatrix[endI][1], deletionMatrix[endI][1])
		    );
      
      for(uint j = 2; j < NCOL; j++){ 
	finalProbSum  = 
	  pu.log10Add(finalProbSum,
		      pu.log10Add(
				  pu.log10Add(matchMatrix[endI][j] , insertionMatrix[endI][j]),
				  pu.log10Add(matchMatrix[endI][j] , deletionMatrix[endI][j])	  
				  )
		      );
      }
      return finalProbSum;
    }
    
    
    void dumpPrior(void){
        dumpMat(priorMatrix, NROW, NCOL);
    }
    void dumpTrans(void){


        std::cerr << "MATCH_TO_MATCH: " << transitions[MATCH_TO_MATCH_I] << std::endl;
        std::cerr << "MATCH_TO_INS: " << transitions[MATCH_TO_INS_I] << std::endl;
        std::cerr << "MATCH_TO_DEL: " <<  transitions[MATCH_TO_DEL_I] << std::endl;
        std::cerr << "INDEL_TO_MATCH: " << transitions[INDEL_TO_MATCH] << std::endl;
        std::cerr << "INS_TO_INS: " << transitions[INS_TO_INS] << std::endl;
        std::cerr << "DEL_TO_DEL: " << transitions[DEL_TO_DEL] << std::endl;
    
    }
    
    void dumpMatchMatrix(void){
        dumpMat(matchMatrix, NROW, NCOL);
    }
    void dumpDeletionMatrix(void){
        dumpMat(deletionMatrix, NROW, NCOL);
    }
    void dumpInsertionMatrix(void){
        dumpMat(insertionMatrix, NROW, NCOL);
    }
};

#endif
