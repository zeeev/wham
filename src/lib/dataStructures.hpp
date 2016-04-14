#ifndef DATASTRUCTURES_H
#define DATASTRUCTURES_H

#include <string>
#include <vector>
#include <map>

#include "api/BamMultiReader.h"

struct regionDat{
  int seqidIndex ;
  int start      ;
  int end        ;
};

struct readPair{
  int          flag;
  int          count;
  BamTools::BamAlignment al1;
  BamTools::BamAlignment al2;
};

struct cigDat{
  int  Length;
  char Type;
};

struct saTag{
  int seqid;
  int pos;
  bool strand;
  std::vector<cigDat> cig;
};


struct breakpoints{
  bool fail              ;
  bool two               ;
  char type              ;
  std::string refBase    ;
  int seqidIndexL        ;
  int seqidIndexR        ;
  std::string seqid      ;
  int merged             ;
  int refined            ;
  int five               ;
  int three              ;
  int svlen              ;
  int collapsed          ;
  int totalSupport       ;
  std::string id              ;
  std::string fives           ;
  std::string threes          ;
  std::vector<std::string> alleles ;
  std::vector<int>         supports;
  std::vector<std::vector<double> > genotypeLikelhoods ;
  std::vector<int>             genotypeIndex           ;
  std::vector<int>             nref                    ;
  std::vector<int>             nalt                    ;
  std::vector<std::string>          sml                ;
  std::vector<std::string>          smr                ;
  int posCIL;
  int posCIH;
  int endCIL;
  int endCIH;

  double lref;
  double lalt;

};


struct node;
struct edge;

struct edge{
  node * L;
  node * R;
  std::map<char,int> support;
};

struct node{
  int   seqid          ;
  int   pos            ;
  int   endSupport     ;
  int   beginSupport   ;
  bool  collapsed      ;
  std::vector <edge *> eds  ;
  std::map<std::string,int> sm   ;
};

struct graph{
  std::map< int, std::map<int, node *> > nodes;
  std::vector<edge *>   edges;
};


struct libraryStats{
  std::map<std::string, double> mus ; // mean of insert length for each indvdual across 1e6 reads
  std::map<std::string, double> sds ;  // standard deviation
  std::map<std::string, double> lq  ;  // 25% of data
  std::map<std::string, double> up  ;  // 75% of the data
  std::map<std::string, double> swm ;
  std::map<std::string, double> sws ;
  std::map<std::string, double> avgD;
  double overallDepth;
} ;


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : a node
 Function does   : visits all the edges and counts the support
 Function returns: int (count of support)
*/

int  getSupport(node * N){

  int support = 0;

  for(std::vector<edge *>::iterator ed = N->eds.begin() ;
      ed != N->eds.end(); ed++){
    support += (*ed)->support['L'];
    support += (*ed)->support['H'];
    support += (*ed)->support['S'];
    support += (*ed)->support['I'];
    support += (*ed)->support['D'];
    support += (*ed)->support['V'];
    support += (*ed)->support['M'];
    support += (*ed)->support['R'];
    support += (*ed)->support['X'];

  }
  return support;
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : left and right node
 Function does   : return true if the left is less than the right
 Function returns: bool
*/

bool sortNodesBySupport(node * L, node * R){
  return (getSupport(L) > getSupport(R)) ;
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : a left and right node
 Function does   : provides the rule for the sort
 Function returns: bool
*/

bool sortNodesByPos(node * L, node * R){
  return (L->pos < R->pos);
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of breakpoint calls
 Function does   : return true if the left is less than the right
 Function returns: bool
*/

bool sortBreak(breakpoints * L, breakpoints * R){

  if(L->seqidIndexL == R->seqidIndexL ){
    if(L->five < R->five){

      return true;
    }
  }
  else{
    if(L->seqidIndexL < R->seqidIndexL){

      return true;
    }
  }
  return false;
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : edge pointer
 Function does   : init
 Function returns: void

*/

void initEdge(edge * e){
  e->L = NULL;
  e->R = NULL;

  e->support['L'] = 0;
  e->support['H'] = 0;
  e->support['S'] = 0;
  e->support['I'] = 0;
  e->support['D'] = 0;
  e->support['V'] = 0;
  e->support['M'] = 0;
  e->support['R'] = 0;
  e->support['X'] = 0;
  e->support['K'] = 0;

}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : refID (int), left position (int), graph
 Function does   : determine if a node is in the graph
 Function returns: bool
*/

inline bool isInGraph(int refID, int pos, graph & lc){

  if(lc.nodes.find(refID) == lc.nodes.end()){
    return false;
  }

  if(lc.nodes[refID].find(pos) != lc.nodes[refID].end()){
    return true;
  }
  return false;
}

#endif
