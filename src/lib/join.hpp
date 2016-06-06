#ifndef JOIN_H
#define JOIN_H

#include <string>
#include <vector>
#include <iostream>


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of strings and separator

 Function does   : joins vector with separator

 Function returns: string

*/

std::string join(std::vector<std::string> & strings, std::string sep){

  std::string joined = *(strings.begin());

  for(std::vector<std::string>::iterator sit = strings.begin()+1;
      sit != strings.end(); sit++){

    joined = joined + sep + (*sit) ;
  }
  return joined;
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of doubles and a seperator

 Function does   : joins vector with separator

 Function returns: string

*/

std::string join(std::vector<double> & ints, std::string sep){

  std::stringstream ss;

  for(std::vector<double>::iterator sit = ints.begin(); sit != ints.end();
      sit++){
    ss << *sit << sep;
  }
  return ss.str();
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of ints and separator
 Function does   : joins vector with separator
 Function returns: string

*/

std::string join(std::vector<int> & ints, std::string sep){

  std::stringstream ss;

  for(std::vector<int>::iterator sit = ints.begin(); 
      sit != ints.end(); sit++){
    ss << *sit << sep;
  }
  return ss.str();
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of strings
 Function does   : joins vector with returns;
 Function returns: string
*/

std::string joinReturn(std::vector<std::string> strings){

  std::string joined = "";

  for(std::vector<std::string>::iterator sit = strings.begin(); 
      sit != strings.end();
      sit++){
    joined = joined + " " + (*sit) + "\n";
  }
  return joined;
}


#endif
