/*

This program was created at:  Wed Aug 12 13:13:52 2015
This program was created by:  Zev Kronenberg

Contact: zev.kronenberg@gmail.com

Organization: Unviersity of Utah
    School of Medicine
    Salt Lake City, Utah


The MIT License (MIT)

Copyright (c) <2015> <Zev Kronenberg>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.


*/

#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "split.h"

using namespace std;

struct options{
       int maxDist     ;
       string filename ;
}globalOpts;

struct svDat{
  string seqid;
  int    five ;
  int    three;
  string id   ;
  string type ; 

  vector<string> lineDat  ;
 
  int ls;
  int rs;

  vector<string> lid;
  vector<string> rid;
  
  map<string, string> info;

};

int roundHalfwayDown(double d)
{
  return ceil(d - 0.5);
}

static const char *optString = "s:f:h";

//-------------------------------   OPTIONS   --------------------------------
int parseOpts(int argc, char** argv)
{
  int opt = 0;
  opt = getopt(argc, argv, optString);
  while(opt != -1){
    switch(opt){
    case 's':
      {
	globalOpts.maxDist  = atoi( ((string)optarg).c_str());
	break;
      }
    case 'f':
      {
	globalOpts.filename = optarg;
	break;
      }
    case 'h':
      {
	break;
      }
    case '?':
      {
	break;
      }
    }
    opt = getopt( argc, argv, optString ); 
  }
  return 1;
}
//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of strings and separator

 Function does   : joins vector with separator

 Function returns: string

*/

string join(vector<string> & strings, string sep){

  string joined = *(strings.begin());

  for(vector<string>::iterator sit = strings.begin()+1; sit != strings.end();
      sit++){

    joined = joined + sep + (*sit) ;
  }
  return joined;
}




//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector svDat pointers

 Function does   : merges calls

 Function returns: void

*/

void mergeAndDump(vector<svDat *> & svs){

  int fiveSum  = 0;
  int threeSum = 0;

  int lSupport = 0;
  int rSupport = 0;

  vector<string> ids;
  vector<string> merged;
  vector<string> refined;
  
  int COLLAPSED = 0;

    for(vector<svDat *>::iterator iz = svs.begin(); iz != svs.end(); iz++){
      fiveSum   += (*iz)->five;
      threeSum  += (*iz)->three;
      
      lSupport += (*iz)->ls;
      rSupport += (*iz)->rs;

      ids.push_back((*iz)->id);

      merged.push_back((*iz)->info["MERGED"]);
      refined.push_back((*iz)->info["REFINED"]);


    
    }
    
    if(svs.size() > 1){
      COLLAPSED = svs.size();
    }

  int  fiveAvg = roundHalfwayDown(double(fiveSum) / double(svs.size()));
  int  threeAvg = roundHalfwayDown(double(threeSum) / double(svs.size()));


  int svlen = threeAvg - fiveAvg;
  if(svs.front()->type.compare("DEL")){
    svlen = -svlen;
  }
  
  map<string, int> ridU;
  map<string, int> lidU;
  
  vector<string> ridV;
  vector<string> lidV;

  for(vector<svDat *>::iterator iz = svs.begin(); iz != svs.end(); iz++){
    for(vector<string>::iterator iy = (*iz)->lid.begin(); iy != (*iz)->lid.end(); iy++ ){
      if(lidU.find(*iy) != lidU.end()){
	lidU[*iy] += 1;
      }
      else{
	lidU[*iy]  = 0;
      }
    }
    for(vector<string>::iterator iy = (*iz)->rid.begin(); iy != (*iz)->rid.end(); iy++ ){
      if(ridU.find(*iy) != ridU.end()){
	ridU[*iy] += 1;
      }
      else{
	ridU[*iy] = 0;
      }
    }
  }

  for(map<string, int>::iterator iz = lidU.begin(); iz != lidU.end(); iz++){
    lidV.push_back(iz->first);
  }
  for(map<string, int>::iterator iz = lidU.begin(); iz != lidU.end(); iz++){
    ridV.push_back(iz->first);
  }
  
  stringstream bedpe;
  
  bedpe << svs.front()->seqid << "\t"
	<< (fiveAvg - 5) << "\t"
	<< (fiveAvg + 5) << "\t"
	<< svs.front()->seqid << "\t"
	<< (threeAvg - 5) << "\t"
	<< (threeAvg + 5) << "\t"
	<< svs.front()->type << ":" << join(ids, "-") << "\t"
	<< ".\t.\t.\t" 
	<< "SVTYPE=" << svs.front()->type 
	<< ";SVLEN=" << svlen << ";ID=" 
	<< join(ids, "-") << ";SUPPORT=" 
	<< lSupport << "," << rSupport 
	<< ";MERGED=" << join(merged, ",") << ";REFINED=" << join(refined, ",") 
	<< ";POS=" << fiveAvg << "," << threeAvg << ";LID=" << join(lidV, ",") 
	<< ";RID=" << join(ridV, ",") << ";" << "COLLAPSED=" << COLLAPSED <<  endl;

  cout << bedpe.str();

  for(vector<svDat *>::iterator iz = svs.begin(); iz != svs.end(); iz++){
    delete *iz;
  }

}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector svDat pointers

 Function does   : merges calls

 Function returns: void

*/

void mergeAndDump(vector<svDat *> & svs, bool){
  map<string, vector< svDat *> > typeGroup;
  for(vector<svDat *>::iterator iz = svs.begin(); iz != svs.end(); iz++){
    typeGroup[(*iz)->type].push_back(*iz);
  }

  for(map<string, vector< svDat *> >::iterator iz = typeGroup.begin(); iz != typeGroup.end(); iz++){
    mergeAndDump(iz->second);
  }
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : pointer to svDat and a stirng
 
 Function does   : parses WHAM-GRPAHING info into a svDat struct

 Function returns: void 

*/
void parseSV(svDat * sv, string & line){

  sv->lineDat =  split(line, "\t");
  vector<string> info = split(sv->lineDat[10], ";");
  for(vector<string>::iterator iz = info.begin(); iz != info.end(); iz++){
    vector<string> key_value = split(*iz, "=");
    sv->info[ key_value[0] ] = key_value[1];
  }
      
  vector<string> pos = split(sv->info["POS"], ',');

  vector<string> support = split(sv->info["SUPPORT"], ',');


  sv->type = sv->info["SVTYPE"];
  sv->id   = sv->info["ID"];
    
  sv->ls = atoi(support.front().c_str());
  sv->rs = atoi(support.back().c_str());

  sv->lid = split(sv->info["LID"], ',');
  sv->rid = split(sv->info["RID"], ',');
  
  sv->five = atoi(pos[0].c_str());
  sv->three = atoi(pos[1].c_str());
  sv->seqid = sv->lineDat[0];
  
}

//-------------------------------    MAIN     --------------------------------
/*
 Comments:
*/

int main( int argc, char** argv)
{
globalOpts.maxDist = 5;
int parse = parseOpts(argc, argv);

 vector<svDat *> SVbuffer;

 ifstream myfile (globalOpts.filename.c_str());
 
 string line;

 if(myfile.is_open()){
   while( getline (myfile, line) ){
     
     svDat * sv ;
     sv = new svDat;
     parseSV(sv, line);

     // nothing to compare against

     if(SVbuffer.size() == 0){
       SVbuffer.push_back(sv);
       continue;
     }

     // swithed seqids

     if(sv->seqid.compare(SVbuffer.back()->seqid) != 0){
       mergeAndDump(SVbuffer, true);
       SVbuffer.clear();
       SVbuffer.push_back(sv);
       continue;
     }
     
     // unsorted

     if(sv->five < SVbuffer.back()->five){
       cerr << endl;
       cerr << "FATAL: SVs are NOT sorted." << endl;
       cerr << "      " << SVbuffer.back()->seqid << " " << SVbuffer.back()->five;
       cerr << "      " << sv->seqid << " " << sv->five << endl;
       
       cerr << endl;
       exit(1);
     }

     if( (sv->five - SVbuffer.back()->five) < globalOpts.maxDist 
	 && (sv->three - SVbuffer.back()->three < globalOpts.maxDist )){
       SVbuffer.push_back(sv);
     }
     else{
       mergeAndDump(SVbuffer, true);
       SVbuffer.clear();
       SVbuffer.push_back(sv);
     }
   } 
   myfile.close();  
 }
 else{
   cerr << endl;
   cerr << "FATAL: Unable to open file." << endl;
   cerr << endl;
   exit(1);
 }


return 0;


}
