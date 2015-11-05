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

const int TRANS_PROB_ARRAY_LENGTH = 6; /*!< Number of possible transitions */

const int MATCH_TO_MATCH_I = 0;
const int MATCH_TO_INS_I   = 1;
const int MATCH_TO_DEL_I   = 2;
const int INDEL_TO_MATCH   = 3;
const int INS_TO_INS       = 4;
const int DEL_TO_DEL       = 5;

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
      
        NROW = nrow;
        NCOL = ncol;
        
        matchMatrix      = new double * [NROW];
        insertionMatrix  = new double * [NROW];
        deletionMatrix   = new double * [NROW];
        priorMatrix      = new double * [NROW];
        transitions      = new double [TRANS_PROB_ARRAY_LENGTH];
        
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
        
        if(haplotype.size() < readSeq.size() ||
           haplotype.size() == 0 || readSeq.size() == 0){
            return false;
        }

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

      transitions[MATCH_TO_MATCH_I] = pu.phredToLog10(1);
      transitions[MATCH_TO_INS_I]   = pu.phredToLog10(45);
      transitions[MATCH_TO_DEL_I]   = pu.phredToLog10(45);
      transitions[INDEL_TO_MATCH]   = pu.phredToLog10(20);
      transitions[INS_TO_INS]       = pu.phredToLog10(30);
      transitions[DEL_TO_DEL]       = pu.phredToLog10(30);
      
      
      return true;
        
    }
    
    void initializeDelMat(void){
        double initialValue = log10(1.0 / NCOL);
        
        for(uint i = 0 ; i < NCOL; i++){
            deletionMatrix[0][i] = initialValue;
        }
    }
    
    void updatecells(void){
        for (uint i = 1; i < NROW;i++) {
            for(uint j = 1; j < NCOL;j++){
              
                double a = matchMatrix[i - 1][j -1]
                +transitions[MATCH_TO_MATCH_I];
                double b = insertionMatrix[i - 1][j -1]
                +transitions[INDEL_TO_MATCH];
                double c = deletionMatrix[i - 1][j -1]
                +transitions[INDEL_TO_MATCH];
                double t = pu.log10Add(a, b);
                matchMatrix[i][j] = priorMatrix[i][j] + pu.log10Add(t, c);
                
                a = matchMatrix[i-1][j]
                +transitions[MATCH_TO_INS_I];
                b = insertionMatrix[i-1][j]
                +transitions[INS_TO_INS];
                
                insertionMatrix[i][j] =
                pu.log10Add(a, b);
                
            
                a = matchMatrix[i][j-1]
                +transitions[MATCH_TO_DEL_I];
                b =  deletionMatrix[i][j-1]
                +transitions[DEL_TO_DEL];
                
                deletionMatrix[i][j]
                = pu.log10Add(a, b);
                
                
            }
        }
    }
    
    double finalLikelihoodCalculation(){
        uint endI = NROW -1;
        double finalProbSum =
        pu.log10Add(matchMatrix[endI][1]
        ,insertionMatrix[endI][1]);
        
        for(uint j =2; j < NCOL; j++){

            double tmp = pu.log10Add(matchMatrix[endI][j],
                                     insertionMatrix[endI][j]);
            
            finalProbSum  = pu.log10Add(finalProbSum, tmp);
            
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
