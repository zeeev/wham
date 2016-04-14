#ifndef STATS_H
#define STATS_H


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of ints
 Function does   : calculates the mean
 Function returns: double
*/

double mean(std::vector<int> & data){

  double sum = 0;

  for(std::vector<int>::iterator it = data.begin(); 
      it != data.end(); it++){
    sum += (*it);
  }
  return sum / data.size();
}

double mean(std::vector<double> & data){

  double sum = 0;

  for(std::vector<double>::iterator it = data.begin(); 
      it != data.end(); it++){
    sum += (*it);
  }
  return sum / data.size();
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of doubles
 Function does   : calculates the var
 Function returns: double

*/

double var(std::vector<double> & data, double mu){
  double variance = 0;

  for(std::vector<double>::iterator it = data.begin(); 
      it != data.end(); it++){
    variance += pow((*it) - mu,2);
  }

  return variance / (data.size() - 1);
}

#endif
