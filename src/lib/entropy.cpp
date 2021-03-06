//
//  entropy.cpp
//  wham
//
//  Created by Zev Kronenberg on 2/19/13.
//  Copyright (c) 2013 Zev Kronenberg. All rights reserved.
//

#include "entropy.h"

using namespace std;

void fastQ::setDNA(string dna){
  fastqSeq = dna;
}

double fastQ::entropy(int size){
  
  int length = fastqSeq.size();
  double window = (double) length - size;
  
  map <string, double> kmer_count;

  for(int i = 0; i < length - size; i++){
    
    string kmer = fastqSeq.substr(i, size);

    kmer_count[kmer]++;
  
  }

  double H = 0;

  map<string, double>::iterator it;
  for(it = kmer_count.begin(); it != kmer_count.end(); it++){
    double p = (it->second / window);
    H += p * log(1/p);
  }
  return H;
}
