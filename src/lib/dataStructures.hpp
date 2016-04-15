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
  std::map<std::string, double> low ;  // 25% of data
  std::map<std::string, double> upr ;  // 75% of the data
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


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : two node pointer
 Function does   : finds if nodes are directly connected
 Function returns: bool
*/

bool connectedNode(node * left, node * right){

  for(std::vector<edge *>::iterator l = left->eds.begin();
      l != left->eds.end(); l++){

    if(((*l)->L->pos == right->pos
	|| (*l)->R->pos == right->pos)
       && (*l)->R->seqid == right->seqid ){
      return true;
    }
  }
  return false;
}

//------------------------------- SUBROUTINE --------------------------------
/*
  Function input  : a vector of edge pointers
  Function does   : does a tree terversal and joins nodes
  Function returns: NA

*/


bool findEdge(std::vector<edge *> & eds, edge ** e, int pos){

  for(std::vector<edge *>::iterator it = eds.begin(); it != eds.end(); it++){
    if((*it)->L->pos == pos || (*it)->R->pos == pos ){
      (*e) = (*it);
      return true;
    }
  }
  return false;
}

//------------------------------- SUBROUTINE --------------------------------
/*
  Function input  : a vector of node pointers, and a postion
  Function does   : removes edges that contain a position
  Function returns: NA

*/


void removeEdges(std::vector<node *> & tree, int pos){

    for(std::vector<node *>::iterator rm = tree.begin();
        rm != tree.end(); rm++){

      std::vector<edge *> tmp;

      for(std::vector<edge *>::iterator e = (*rm)->eds.begin();
        e != (*rm)->eds.end(); e++){
      if( (*e)->L->pos != pos && (*e)->R->pos != pos  ){
        tmp.push_back((*e));
      }
      else{

      }
    }
    (*rm)->eds.clear();
    (*rm)->eds.insert((*rm)->eds.end(), tmp.begin(), tmp.end());
  }
}

//------------------------------- SUBROUTINE --------------------------------
/*
  Function input  : a vector of node pointers
  Function does   : does a tree terversal and joins nodes
  Function returns: NA

*/

void joinNodes(node * L, node * R, std::vector<node *> & tree){

  // quantifying which node has more support

  int lSupport = 0;
  int rSupport = 0;

  for(std::vector<edge *>::iterator le = L->eds.begin(); le != L->eds.end(); le++){
    lSupport += (*le)->support['D'] + (*le)->support['S'] + (*le)->support['I']
      +  (*le)->support['H'] +  (*le)->support['L'] + (*le)->support['X']
      + (*le)->support['M'] + (*le)->support['R'] + (*le)->support['V'] ;
  }
  for(std::vector<edge *>::iterator re = R->eds.begin(); re != R->eds.end(); re++){
    lSupport += (*re)->support['D'] + (*re)->support['S'] + (*re)->support['I']
      +  (*re)->support['H'] +  (*re)->support['L'] +  (*re)->support['X']
      +  (*re)->support['M'] +  (*re)->support['R'] +  (*re)->support['V'] ;
  }

  if(lSupport <= rSupport){

    L->collapsed = true;

    for(std::map<std::string,int>::iterator iz = L->sm.begin();
        iz != L->sm.end(); iz++){

      if(R->sm.find(iz->first) != R->sm.end()){
        R->sm[iz->first] += iz->second;
      }
      else{
        R->sm[iz->first] = iz->second;
      }
    }

    for(std::vector<edge *>::iterator lc =  L->eds.begin();
        lc != L->eds.end(); lc++){


      edge * e;

      int otherP = (*lc)->L->pos;

      if(L->pos  == otherP){
        otherP = (*lc)->R->pos;
      }
      if(findEdge(R->eds, &e,  otherP)){

        e->support['I'] += (*lc)->support['I'];
        e->support['D'] += (*lc)->support['D'];
        e->support['S'] += (*lc)->support['S'];
        e->support['H'] += (*lc)->support['H'];
        e->support['L'] += (*lc)->support['L'];
        e->support['R'] += (*lc)->support['R'];
        e->support['M'] += (*lc)->support['M'];
        e->support['V'] += (*lc)->support['V'];
        e->support['X'] += (*lc)->support['X'];

      }
      else{
        if((*lc)->L->pos == otherP){
          (*lc)->R = R;
        }
        else{
          (*lc)->L = R;
        }

        R->eds.push_back(*lc);
      }

    }
    removeEdges(tree, L->pos);

  }
  else{

    //    cerr << "Joining right" << endl;
    R->collapsed = true;

    for(std::map<std::string,int>::iterator iz = R->sm.begin();
        iz != R->sm.end(); iz++){

      if(L->sm.find(iz->first) != L->sm.end()){
        L->sm[iz->first] += iz->second;
      }
      else{
        L->sm[iz->first] = iz->second;
      }
    }

    for(std::vector<edge *>::iterator lc =  R->eds.begin();
        lc != R->eds.end(); lc++){

      edge * e;

      int otherP = (*lc)->L->pos;

      if(R->pos  == otherP){
        otherP = (*lc)->R->pos;
      }
      if(findEdge(L->eds, &e,  otherP)){
        e->support['I'] += (*lc)->support['I'];
        e->support['D'] += (*lc)->support['D'];
        e->support['S'] += (*lc)->support['S'];
        e->support['H'] += (*lc)->support['H'];
        e->support['L'] += (*lc)->support['L'];
        e->support['M'] += (*lc)->support['M'];
        e->support['X'] += (*lc)->support['X'];
        e->support['V'] += (*lc)->support['V'];
      }
      else{
        if((*lc)->L->pos == otherP){
          (*lc)->R = L;
        }
        else{
          (*lc)->L = L;
        }
        L->eds.push_back(*lc);

      }
    }
    removeEdges(tree, R->pos);
  }

}

//------------------------------- SUBROUTINE --------------------------------
/*
  Function input  : a vector of node pointers
  Function does   : does a tree terversal and joins nodes
  Function returns: NA
*/

void collapseTree(std::vector<node *> & tree){

    std::vector<node *> tmp;

    for(std::vector<node *>::iterator tr = tree.begin();
        tr != tree.end(); tr++){
        if((*tr)->collapsed){
            continue;
        }
        for(std::vector<node *>::iterator tt = tree.begin();
        tt != tree.end(); tt++){
            if((*tt)->collapsed){
                continue;
            }
            if( (*tr)->pos == (*tt)->pos ){
                continue;
            }

            if(abs( (*tr)->pos - (*tt)->pos ) < 20
               && ! connectedNode((*tr), (*tt))
               && ( (*tr)->seqid == (*tt)->seqid ) ){
                joinNodes((*tr), (*tt), tree);
            }
        }
    }
    for(std::vector<node *>::iterator tr = tree.begin();
        tr != tree.end(); tr++){
        if((*tr)->collapsed){
        }
        else{
            tmp.push_back((*tr));
        }
    }
    tree.clear();
    tree.insert(tree.end(), tmp.begin(), tmp.end());
}



#endif
