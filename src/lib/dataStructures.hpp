#ifndef DATASTRUCTURES_H
#define DATASTRUCTURES_H

#include <string>
#include <vector>
#include <map>
#include <string>

#include "api/BamMultiReader.h"
#include "fastahack/Fasta.h"
#include "ssw_cpp.h"

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
    int mapQ;
    int NM;
    std::vector<cigDat> cig;
};

struct node;
struct edge;

struct edge{
  node * L;
  node * R;
  std::map<char,int> support;
};

struct node{
  int   seqid                  ;
  int   pos                    ;
  std::vector <edge *>      eds;
  std::map<std::string, int> sm;
};

struct graph{
  std::map< int, std::map<int, node *> > nodes;
  std::vector<edge *>   edges;
};

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : a node
 Function does   : visits all the edges and counts the support
 Function returns: int (count of support)
*/

int getSupport(node * N){

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
        support += (*ed)->support['A'];
        support += (*ed)->support['Z'];
    }
    return support;
}


class breakpoint{
private:

    std::map<char, std::string > typeMap;
    std::vector<std::string>        refs;
    std::vector<std::string>    seqNames;
    std::map<std::string, int>       sms;

    char type;
    std::string typeName;

    int totalGraphWeight ;

    bool paired   ;
    bool bnd      ;
    bool masked   ;
    bool clustered;
    bool print    ;

    long int length;

    int nClustered       ;

    double totalCount    ;
    double tooFarCount   ;
    double tooCloseCount ;
    double internalDCount;
    double splitCount    ;
    double evertCount    ;
    double ssCount       ;
    double delCount      ;
    double insCount      ;
    double dupCount      ;
    double invCount      ;
    double traCount      ;
    double clusterFrac   ;
    double avgDist       ;

    void _count(node * n, std::map<edge *, int> & lu){
        for(std::vector<edge *>::iterator ed = n->eds.begin() ;
            ed != n->eds.end(); ed++){

            if(lu.find(*ed) != lu.end()){
                continue;
            }
            lu[*ed] = 1;

            if((*ed)->L->seqid == (*ed)->R->seqid){
                ssCount  += (*ed)->support['M'] ;
                ssCount  += (*ed)->support['A'] ;
                ssCount  += (*ed)->support['R'] ;
                splitCount  += (*ed)->support['S'] ;
                splitCount  += (*ed)->support['V'] ;
                splitCount  += (*ed)->support['Z'] ;
                tooFarCount += (*ed)->support['H'] ;
                tooFarCount += (*ed)->support['M'] ;
                evertCount  += (*ed)->support['X'] ;
                internalDCount += (*ed)->support['D'] ;
                delCount += (*ed)->support['H'] ;
                delCount += (*ed)->support['D'] ;
                delCount += (*ed)->support['S'] ;
                insCount += (*ed)->support['I'] ;
                insCount += (*ed)->support['R'] ;
                insCount += (*ed)->support['L'] ;
                invCount += (*ed)->support['M'] ;
                invCount += (*ed)->support['A'] ;
                invCount += (*ed)->support['R'] ;
                invCount += (*ed)->support['V'] ;
                dupCount += (*ed)->support['Z'] ;
                dupCount += (*ed)->support['S'] ;
                dupCount += (*ed)->support['X'] ;
            }
            else{
                traCount += (*ed)->support['H'] ;
                traCount += (*ed)->support['D'] ;
                traCount += (*ed)->support['S'] ;
                traCount += (*ed)->support['I'] ;
                traCount += (*ed)->support['R'] ;
                traCount += (*ed)->support['L'] ;
                traCount += (*ed)->support['M'] ;
                traCount += (*ed)->support['A'] ;
                traCount += (*ed)->support['V'] ;
                traCount += (*ed)->support['S'] ;
                traCount += (*ed)->support['X'] ;
                traCount += (*ed)->support['Z'] ;
            }
            totalCount += (*ed)->support['H'];
            totalCount += (*ed)->support['D'];
            totalCount += (*ed)->support['S'];
            totalCount += (*ed)->support['I'];
            totalCount += (*ed)->support['R'];
            totalCount += (*ed)->support['L'];
            totalCount += (*ed)->support['M'];
            totalCount += (*ed)->support['A'];
            totalCount += (*ed)->support['V'];
            totalCount += (*ed)->support['V'];
            totalCount += (*ed)->support['S'];
            totalCount += (*ed)->support['X'];
            totalCount += (*ed)->support['Z'];
        }
    }

    void _processInternal(void){
        if(nodeL->seqid != nodeR->seqid ){
            if(nodeL->seqid > nodeR->seqid){
                node * tmp;
                tmp   = nodeL;
                nodeL = nodeR;
                nodeR = tmp;
            }
        }
        else{
            if(nodeL->pos > nodeR->pos){
                node * tmp;
                tmp   = nodeL;
                nodeL = nodeR;
                nodeR = tmp;
            }
            this->length = this->nodeR->pos - nodeL->pos;
        }
    }

public:

    node * nodeL;
    node * nodeR;

    breakpoint(void): type('D')
                    , totalGraphWeight(0)
                    , paired(false)
                    , bnd(false)
                    , masked(false)
                    , clustered(false)
                    , print(true)
                    , length(0)
                    , nClustered(0)
                    , totalCount(0)
                    , tooFarCount(0)
                    , tooCloseCount(0)
                    , internalDCount(0)
                    , splitCount(0)
                    , evertCount(0)
                    , ssCount(0)
                    , delCount(0)
                    , insCount(0)
                    , dupCount(0)
                    , invCount(0)
                    , traCount(0)
                    , clusterFrac(0)
                    , avgDist(0)
                    , nodeL(NULL)
                    , nodeR(NULL){

        typeMap['D'] = "DEL";
        typeMap['U'] = "DUP";
        typeMap['T'] = "BND";
        typeMap['I'] = "INS";
        typeMap['V'] = "INV";

    }
    int getNClustered(void){
        return this->nClustered;
    }
    double getSMSupport(std::string & s){
        if(sms.find(s) == sms.end()){
            return 0;
        }
        else{
            return sms[s];
        }
    }
    double getDelCount(void){
        return this->delCount;
    }
    double getDupCount(void){
        return this->dupCount;
    }
    double getInvCount(void){
        return this->invCount;
    }
    double getTraCount(void){
        return this->traCount;
    }
    double getInsCount(void){
        return this->insCount;
    }
    double getClustFrac(void){
        return this->clusterFrac;
    }
    double getSameStrandCount(void){
        return this->ssCount;
    }
    double getTooFarCount(void){
        return this->tooFarCount;
    }
    double getTooCloseCount(void){
        return this->tooCloseCount;
    }
    double getInternalDelCount(void){
        return this->internalDCount;
    }
    double getEvertCount(void){
        return evertCount;
    }
    double getSplitReadCount(void){
        return splitCount;
    }

    char getType(void){
        return this->type;
    }
    double getAvgDist(void){
        return this->avgDist;
    }
    long int getLength(void){
        return length;
    }
    void setGoodPair(bool t){
        this->paired = t;
    }
    bool IsMasked(void){
        return this->masked;
    }
    bool IsPrint(void){
        return this->print;
    }
    void unsetMasked(void){
        this->masked = false;
    }
    void setMasked(void){
        this->masked = true;
    }
    void unSetPrint(void){
        this->print = false;
    }
    void setBadPair(void){
        this->paired = true;
    }
    void setTotalSupport(int t){
        this->totalCount = t;
    }
    bool IsBadPair(void){
        return this->paired;
    }
    int getTotalSupport(void){
        return this->totalCount;
    }

    bool add(node * n){
        if(nodeL != NULL && nodeR != NULL){
            return false;
        }
        if( nodeL == NULL ){
            nodeL = n;
            return true;
        }
        else{
            nodeR = n;
            _processInternal();
            return true;
        }
    }
    bool countSupportType(void){
        if(nodeL == NULL || nodeR == NULL){
            return false;
        }
        // make sure we dont double count edges
        std::map<edge *, int> lu;

        _count(nodeL, lu);
        _count(nodeR, lu);
        return true;
    }
    bool calcType(void){

        if(nodeL == NULL || nodeR == NULL){
            return false;
        }

        int maxN = delCount;

        if(dupCount > maxN){
            type = 'U';
            maxN = dupCount;
        }
        if(invCount > maxN){
            type = 'V';
            maxN = invCount;
        }
        if(insCount > maxN){
            type = 'I';
            maxN = insCount;
        }
        if(traCount > maxN){
            type = 'T';
            maxN = traCount;
        }
        typeName = typeMap[type];
        return true;
    }

    void getRefBases(std::string & reffile, std::map<int, string> & lookup){


        FastaReference rs;
        rs.open(reffile);

        refs.push_back(rs.getSubSequence(lookup[nodeL->seqid],
                                          nodeL->pos,  1));
        if(nodeL->seqid != nodeR->seqid){
            refs.push_back(rs.getSubSequence(lookup[nodeR->seqid],
                                              nodeR->pos,  1));
        }

        seqNames.push_back(lookup[nodeL->seqid]);
        seqNames.push_back(lookup[nodeR->seqid]);

    }


    friend std::ostream& operator<<(std::ostream& out,
                                    const breakpoint& foo);


    bool dupClusterCheck(void){
        return dupClusterCheck(nodeL, nodeR);
    }
    bool invClusterCheck(void){
        return invClusterCheck(nodeL, nodeR);
    }
    bool delClusterCheck(void){
        return delClusterCheck(nodeL, nodeR);
    }

    bool dupClusterCheck(node * left, node * right){

        if(left == NULL || right == NULL){
            return false;
        }

        if(type != 'U'){
            return false;
        }
        if(left->seqid != right->seqid ){
            return false;
        }

        if(left->pos > right->pos){
            node * tmp;
            tmp  = left;
            left = right;
            right = tmp;
        }

        int countY  = 0;
        int countN  = 0;
        long int distSum = 0;
        int distCount    = 0;

        for( std::vector<edge*>::iterator it = left->eds.begin();
             it != left->eds.end(); it++){
            // same strand too far or too close
            if( (*it)->support['X'] == 0 ){
                continue;
            }
            // don't want the same node or non leaf
            node * tmp = (*it)->L;
            if(tmp->eds.size() > 1 || tmp == left ){
                tmp = (*it)->R;
            }
            if(tmp->eds.size() > 1){
                continue;
            }

            long int distL = abs(left->pos  - tmp->pos);
            long int distR = abs(right->pos - tmp->pos);

            if(distR < distL){
                countY += (*it)->support['X'];
                distSum += distR;
                distCount += 1;
            }
            else{
                countN += (*it)->support['X'];
            }
        }
        for( std::vector<edge*>::iterator it = right->eds.begin();
             it != right->eds.end(); it++){
            // same strand too far or too close
            if( (*it)->support['X'] == 0 ){
                continue;
            }
            // don't want the same node or non leaf
            node * tmp = (*it)->L;
            if(tmp->eds.size() > 1 || tmp == nodeR ){
                tmp = (*it)->R;
            }
            if(tmp->eds.size() > 1){
                continue;
            }

            long int distL = abs(left->pos  - tmp->pos);
            long int distR = abs(right->pos - tmp->pos);

            if(distR < distL){
                countY += (*it)->support['X'];
                distSum += distR;
                distCount += 1;
            }
            else{
                countN += (*it)->support['X'];
            }
        }

        int sum = countY + countN;

        if(sum > 0){
            clusterFrac = double(countY) / double(sum);
            return true;
        }
        if(distCount > 0){
            avgDist = double(distSum) / double(distCount);
        }
        if(sum > 0 && distCount > 0){
            return true;
        }
        return false;
    }

    bool invClusterCheck(node * left, node * right){

        if(left == NULL || right == NULL){
            return false;
        }

        if(type != 'V'){
            return false;
        }
        if(left->seqid != right->seqid ){
            return false;
        }

        if(left->pos > right->pos){
            node * tmp;
            tmp = left;
            left = right;
            right = tmp;
        }

        int countY       = 0;
        int countN       = 0;
        long int distSum = 0;
        int distCount    = 0;

        for( std::vector<edge*>::iterator it = left->eds.begin();
             it != left->eds.end(); it++){
            // same strand too far or too close
            if((*it)->support['R'] == 0 && (*it)->support['A']){
                continue;
            }

            // don't want the same node or non leaf
            node * tmp = (*it)->L;
            if(tmp->eds.size() > 1 || tmp == left ){
                tmp = (*it)->R;
            }
            if(tmp->eds.size() > 1){
                continue;
            }

            if(tmp->pos < left->pos || tmp->pos > right->pos){
                continue;
            }

            long int distL = abs(left->pos  - tmp->pos);
            long int distR = abs(right->pos - tmp->pos);

            if(distR < distL){
                countY += (*it)->support['R'];
                countY += (*it)->support['A'];
                countY += (*it)->support['M'];

                distSum += distR;
                distCount += 1;
            }
            else{
                countN += (*it)->support['R'];
                countN += (*it)->support['A'];
                countN += (*it)->support['M'];
            }
        }
        for(std::vector<edge*>::iterator it = right->eds.begin();
            it != right->eds.end(); it++){

            // same strand too far or too close
            if((*it)->support['R'] == 0 && (*it)->support['A']){
                continue;
            }

            // don't want the same node or non leaf
            node * tmp = (*it)->L;
            if(tmp->eds.size() > 1 || tmp ==  left ){
                tmp = (*it)->R;
            }
            if(tmp->eds.size() > 1){
                continue;
            }

            if(tmp->pos < left->pos || tmp->pos > right->pos){
                continue;
            }


            // same strands need to be internal to the inversion
            if(tmp->pos > left->pos
               && tmp->pos < right->pos){

                long int distL = abs(left->pos - tmp->pos);
                long int distR = abs(right->pos - tmp->pos);

                if(distR > distL){
                    countY  += (*it)->support['R'];
                    countY  += (*it)->support['A'];
                    countY  += (*it)->support['M'];
                    distSum += distR;
                    distCount += 1;
             }
                else{
                    countN += (*it)->support['R'];
                    countN += (*it)->support['A'];
                    countN += (*it)->support['M'];
                }
            }
        }
        int sum = countY + countN;

        if(sum > 0){
            clusterFrac = double(countY) / double(sum);
        }
        if(distCount > 0){
            avgDist = double(distSum) / double(distCount);
        }
        if(sum > 0 && distCount > 0){
            return true;
        }
        return false;
    }

    void loadSMSupport(void){

        for(std::map<std::string, int>::iterator sm = nodeL->sm.begin();
            sm!= nodeL->sm.end(); sm++){

            if(sms.find(sm->first) == sms.end() ){
                sms[sm->first] = sm->second;
            }
            else{
                sms[sm->first] += sm->second;
            }
        }
        for(std::map<std::string, int>::iterator sm = nodeR->sm.begin();
            sm!=nodeR->sm.end(); sm++){

            if(sms.find(sm->first) == sms.end() ){
                sms[sm->first] = sm->second;
            }
            else{
                sms[sm->first] += sm->second;
            }

        }
    }


    bool delClusterCheck(node * left, node * right){

        if(left == NULL || right == NULL){
            return false;
        }

        if(type != 'D'){
            return false;
        }
        if(left->seqid != right->seqid ){
            return false;
        }

        if(left->pos > right->pos){
            node * tmp;
            tmp = left;
            left = right;
            right = tmp;
        }

        int countY       = 0;
        int countN       = 0;
        long int distSum = 0;
        int distCount    = 0;

        /* loop over 5' node (always sorted)
           - check if left pos = 5' node, if yes switch pos
           - count those nodes to the right of 3' node
         */

        for( std::vector<edge*>::iterator it = left->eds.begin();
             it != left->eds.end(); it++){

            if((*it)->support['H'] == 0){
                continue;
            }
            long int five = (*it)->L->pos;
            if((*it)->L == left){
                five = (*it)->R->pos;
            }

            if(five > right->pos ){
                countY += (*it)->support['H'];
                distSum += abs(five - right->pos);
                distCount += 1;
            }
        }

        /* loop over 3' node (always sorted)
           - check if right pos = 3' node, if yes switch pos
           - count those nodes to the left of 5' node
        */

        for( std::vector<edge*>::iterator it = right->eds.begin();
             it != right->eds.end(); it++){

            if((*it)->support['H']  == 0){
                continue;
            }
            long int three = (*it)->L->pos;
            if((*it)->L == right){
                three = (*it)->R->pos;
            }

            if( three < left->pos ){
                countY += (*it)->support['H'];
                distSum += abs(three - left->pos);
                distCount += 1;
            }
        }

        int sum = countY + countN;

        if(sum > 0){
            clusterFrac = double(countY) / double(sum);
        }
        if(distCount > 0){
            avgDist = double(distSum) / double(distCount);
        }

        if(sum > 0 && distCount > 0){
            nClustered = countY;
            return true;
        }
        return false;
    }
};

std::ostream& operator<<(std::ostream& out, const breakpoint & foo){

    long int start = foo.nodeL->pos ;
    long int end   = foo.nodeR->pos ;

    long int len = abs(foo.nodeR->pos - foo.nodeL->pos);

    if(foo.type == 'I'){
        end = start;
    }
    if(foo.type == 'D'){
        len = -1*len;
    }

    double sum = foo.delCount + foo.dupCount
        + foo.traCount + foo.invCount + foo.insCount;


    stringstream ss;
    int index = 0;
    for(std::map<std::string, int>::const_iterator sm = foo.sms.begin();
        sm != foo.sms.end(); sm++){
        if(index == 0){
            ss << sm->first;
        }
        else{
            ss << "," << sm->first;
        }
        index += 1;
    }


    if(foo.type != 'T' && (foo.nodeL->seqid == foo.nodeR->seqid)){
        out << foo.seqNames.front() << "\t"
            << start + 1            << "\t"
            << "."              << "\t"
            << foo.refs.front() << "\t"
            << "<" << foo.typeName << ">\t"
            << "."              << "\t"
            << "PASS"           << "\t"
            << "A="             << foo.totalCount
            << ";CIEND=-10,10;CIPOS=-10,10"
            << ";CF="           << foo.clusterFrac
            << ";CW="           << (foo.delCount / sum) << ","
            << (foo.dupCount / sum) << ","
            << (foo.invCount / sum) << ","
            << (foo.insCount / sum) << ","
            << (foo.traCount / sum)
            << ";D="            << foo.delCount
            << ";DI="           << foo.avgDist
            << ";END="          << end + 1
            << ";EV="           << foo.evertCount
            << ";I="            << foo.insCount
            << ";SR="           << foo.splitCount
            << ";SS="           << foo.ssCount
            << ";SVLEN="        << len
            << ";SVTYPE="       << foo.typeName
            << ";T="            << foo.traCount
            << ";TAGS="         << ss.str()
            << ";TF="           << foo.tooFarCount
            << ";U="            << foo.dupCount
            << ";V="            << foo.invCount;
    }

    return out;
}


struct libraryStats{
    std::map<std::string, double> mus ; // mean of insert length for each indvdual across 1e6 reads
    std::map<std::string, double> sds ; // standard deviation
    std::map<std::string, double> low ; // 25% of data
    std::map<std::string, double> upr ; // 75% of the data
    std::map<std::string, double> swm ;
    std::map<std::string, double> sws ;
    std::map<std::string, double> avgD;
    std::map<std::string, std::map<int, long int> > mqD;
    double overallDepth;
} ;

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


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : edge pointer
 Function does   : init
 Function returns: void

*/

void initEdge(edge * e){
  e->L = NULL;
  e->R = NULL;

  // mate too close
  e->support['L'] = 0;
  // mate too far
  e->support['H'] = 0;
  // split read
  e->support['S'] = 0;
  // split read different strand (split)
  e->support['V'] = 0;
  // insertion
  e->support['I'] = 0;
  // deletion
  e->support['D'] = 0;
  // same strand too far
  e->support['M'] = 0;
  // same strand too close
  e->support['R'] = 0;

  // everted read pairs
  e->support['X'] = 0;
  e->support['K'] = 0;
  // same strand
  e->support['A'] = 0;
  // split reads on wrong side of soft clip
  e->support['z'] = 0;

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

#endif
