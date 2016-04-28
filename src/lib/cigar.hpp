#ifndef CIGAR_H
#define CIGAR_H

#include <vector>
#include <sstream>
#include <string>
#include "dataStructures.hpp"
#include "api/BamMultiReader.h"

#define BAM_CIGAR_SHIFT 4
#define BAM_CIGAR_MASK  ((1 << BAM_CIGAR_SHIFT) - 1)


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of BamTools::CigarOp
 Function does   : joins vector " ";
 Function returns: string
*/

std::string joinCig(std::vector<BamTools::CigarOp> strings){

  std::stringstream joined ;

  for(std::vector<BamTools::CigarOp>::iterator sit = strings.begin();
      sit != strings.end(); sit++){
    joined  << (*sit).Length << (*sit).Type;
  }
  return joined.str();
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of BamTools::BamTools::CigarOp, int length
 Function does   : looks for clips greater than or equal to length
 Function returns: bool
*/

bool areBothClipped(std::vector<BamTools::CigarOp> & ci){

  if(ci.front().Type == 'S' && ci.back().Type == 'S'){
    return true;
  }
  return false;
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of BamTools::BamTools::CigarOp, int length
 Function does   : looks for clips greater than or equal to length
 Function returns: bool
*/

bool IsLongClip(std::vector<BamTools::CigarOp> & ci,
		       unsigned int len){

  if(ci.front().Type == 'S' && ci.front().Length >= len){
    return true;
  }
  if(ci.back().Type == 'S' && ci.back().Length >= len){
    return true;
  }
  return false;
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector<cigDat>, string
 Function does   : loads cigar data into a string
 Function returns: NA

*/

void parseCigar(std::vector<cigDat> & parsedCigar, std::string cigar){

  unsigned int spot = 0;

  for(unsigned int i = 0; i < cigar.size(); i++){
    if(int(cigar[i] > 57)){
      cigDat tup;
      tup.Length = atoi(cigar.substr(spot, i-spot).c_str());
      tup.Type  = cigar[i];
      parsedCigar.push_back(tup);
    }
  }
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector<CigarOp>
 Function does   : calculates the number of matching bases
 Function returns: unsigned int, total number of matching bases
*/

int match(std::vector<BamTools::CigarOp> & co){

  int m = 0;

  for(std::vector<BamTools::CigarOp>::iterator it = co.begin();
      it != co.end(); it++){
    if(it->Type == 'M'){
      m += it->Length;
    }
  }
  return m;
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : int pointer, vector cigDat
 Function does   : finds end position
 Function returns: void
*/
void endPos(std::vector<cigDat> & cigs, int * pos){

    for(std::vector<cigDat>::iterator it = cigs.begin();
        it != cigs.end(); it++){

        switch( (*it).Type ){
        case 'M':
            {
                *pos += (*it).Length;
                break;
            }
        case 'X':
            {
                *pos += (*it).Length;
                break;
            }
        case 'D':
            {
                *pos += (*it).Length;
                break;
            }
        case '=':
            {
                *pos += (*it).Length;
                break;
            }
        case 'N':
            {
                *pos += (*it).Length;
                break;
            }
        default:
            break;
        }
        // WARNING: this needs to be double checked
        // *pos -= 1;
    }
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : SA tag string, and vector of saTag
 Function does   : parses string and converts positions to BAM
 Function returns: NA
*/

void parseSA(std::vector<saTag> & parsed,
             std::string tag,
             std::string type,
             std::map<std::string, int> & il){

    std::vector<std::string> sas = split(tag, ';');

  for(unsigned int i = 0 ; i < sas.size() -1 ; i++){

      saTag sDat;

      std::vector<std::string> sat = split (sas[i], ',');

      if(sat.size() != 6 && type == "SA"){
          std::cerr << "FATAL: failure to parse SA tag" << std::endl;
          exit(1);
      }

      sDat.seqid = il[sat[0]];

      if(type == "SA"){
          sDat.pos   = atoi(sat[1].c_str()) - 1;
          sDat.mapQ  = atoi(sat[4].c_str());
          sDat.NM    = atoi(sat[5].c_str());

          if(sat[2].compare("-") == 0){
          sDat.strand = true;
          }
          else{
              sDat.strand = false;
          }
          parseCigar(sDat.cig, sat[3]);
      }
      else if(type == "XP"){
          char strand = sat[1][0];
          sat[1].erase(0,1);
          sDat.pos   = atoi(sat[1].c_str()) - 1;
          if(strand == '-'){
              sDat.strand = true;
          }
          else{
              sDat.strand = false;
          }
          parseCigar(sDat.cig, sat[2]);
      }
      parsed.push_back(sDat);
  }

}


#endif
