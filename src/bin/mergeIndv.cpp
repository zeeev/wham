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
#include <getopt.h>

using namespace std;

struct options{
       int maxDist     ;
       string filename ;
}globalOpts;

struct svDat{
  string seqid;
  string refBase;
  int    five ;
  int    three;
  string id   ;
  string type ; 
  int collapsed;

  vector<int>          fives;
  vector<int>         threes;
  vector<string> genotypeDat;
  vector<string> lineDat    ;
 
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

int NSAMP = 0;


void printVersion(void){
  cerr << "Version: " << VERSION << endl;
  cerr << "Contact: zev.kronenberg [at] gmail.com " << endl;
  cerr << "Notes  : -If you find a bug, please open a report on github!" << endl;
  cerr << endl;
}


void printHelp(void){
  //------------------------------- XXXXXXXXXX --------------------------------
  cerr << " Usage:  " << endl;
  cerr << "       mergeIndv -f my.combined.vcf -s 5" << endl;
  cerr << endl;
  cerr << " Required:  " << endl;
  //------------------------------- XXXXXXXXXX --------------------------------

  cerr << "          -f - <STRING> - A vcf file from WHAM-GRAPHENING" << endl;
  cerr << endl;
  cerr << " Optional:  " << endl;
  cerr << "          -s - <INT>   - Merge SVs with both breakpoints N BP away [5] " << endl;
  cerr << endl;
  printVersion();
}
//-------------------------------   OPTIONS   --------------------------------

void printHeader(void){
  
  stringstream header;

  header << "##fileformat=VCFv4.2" << endl;
  header << "##source=WHAM-GRAPHENING-mergeIndvs:" << VERSION << endl;
  header << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">" << endl;
  header << "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">" << endl;
  header << "##INFO=<ID=ID,Number=1,Type=String,Description=\"Unique hexadecimal identifier\">" << endl;
  header << "##INFO=<ID=SUPPORT,Number=2,Type=Integer,Description=\"Number of reads supporting POS and END breakpoints\">" << endl;
  header << "##INFO=<ID=MERGED,Number=1,Type=Integer,Description=\"SV breakpoints were joined without split read support 0=false 1=true\">" << endl;
  header << "##INFO=<ID=REFINED,Number=1,Type=Integer,Description=\"SV breakpoints were refined based on SW alignment 0=false 1=true\">" << endl;
  header << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">" << endl;
  header << "##INFO=<ID=POS,Number=2,Type=String,Description=\"POS and END\">" << endl;
  header << "##INFO=<ID=FIVE,Number=.,Type=Integer,Description=\"collapsed POS\">" << endl;
  header << "##INFO=<ID=THREE,Number=.,Type=Integer,Description=\"collapsed END\">" << endl;
  header << "##INFO=<ID=LID,Number=.,Type=String,Description=\"POS breakpoint support came from SM, independent of genotype\">" << endl;
  header << "##INFO=<ID=RID,Number=.,Type=String,Description=\"END breakpoint support came from SM, independent of genotype\">" << endl;
  header << "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">" << endl;
  header << "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">" << endl;
  header << "##INFO=<ID=COLLAPSED,Number=1,Type=Integer,Description=\"Number of SV calls merged into record\">" << endl;
  header << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
  header << "##FORMAT=<ID=GL,Number=.,Type=Float,Description=\"Genotype likelihoods comprised of comma separated floating point log10-scaled likelihoods for all possible genotypes given the set of alleles defined in the REF and ALT fields.\">" << endl;
  header << "##FORMAT=<ID=AS,Number=1,Type=Integer,Description=\"Number of reads that align better to ALT allele\">" << endl;
  header << "##FORMAT=<ID=RS,Number=1,Type=Integer,Description=\"Number of reads that align better to REF allele\">" << endl;
  header << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNONE" << endl;
  
  cout << header.str();
}

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
	printHelp();
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
 Function input  : vector of doubles

 Function does   : calculates the var

 Function returns: double

*/

double var(vector<int> & data, double mu){
  double variance = 0;

  for(vector<int>::iterator it = data.begin(); it != data.end(); it++){
    variance += pow((*it) - mu,2);
  }

  if(variance == 0){
    return 0;
  }

  return variance / (data.size() - 1);
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of ints

 Function does   : calculates the mean

 Function returns: double

*/

double mean(vector<int> & data){

  double sum = 0;

  for(vector<int>::iterator it = data.begin(); it != data.end(); it++){
    sum += (*it);
  }
  return sum / data.size();
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of strings and separator

 Function does   : joins vector with separator

 Function returns: string

*/

string join(vector<int> & ints, string sep){
  stringstream ss ;
  
  ss << ints.front();
  
  for(int i = 1 ; i < ints.size(); i++){
    ss << sep << ints[i];
  }

  return ss.str();

}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of strings and separator

 Function does   : joins vector with separator

 Function returns: string

*/

string join(vector<string> & strings, string sep){

  string joined = *(strings.begin());

  for(vector<string>::iterator sit = strings.begin()+1; 
      sit != strings.end(); sit++){

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

  long int fiveSum  = 0;
  long int threeSum = 0;

  int lSupport = 0;
  int rSupport = 0;

  vector<string> ids;
  vector<string> merged;
  vector<string> refined;

  char hex[8 + 1];
  for(int i = 0; i < 8; i++) {
    sprintf(hex + i, "%x", rand() % 16);
  }

  stringstream xx ;
  xx << hex;

  string newid = xx.str();

  if(svs.size() == 1){
    newid = svs.front()->id;
  }
   
  int COLLAPSED = 0;

  vector<int> FIVEPOS;
  vector<int> THREEPOS;

  for(vector<svDat *>::iterator iz = svs.begin(); iz != svs.end(); iz++){
    
    for(vector<int>::iterator it = (*iz)->fives.begin(); 
	it != (*iz)->fives.end(); it++){
      
      FIVEPOS.push_back(*it);
    }
    for(vector<int>::iterator it = (*iz)->threes.begin();
	it !=(*iz)->threes.end(); it++){
      
      THREEPOS.push_back(*it);
    }
    
    fiveSum   += (*iz)->five;
    threeSum  += (*iz)->three;
    
    lSupport += (*iz)->ls;
    rSupport += (*iz)->rs;
    
    ids.push_back((*iz)->id);
    
    merged.push_back((*iz)->info["MERGED"]);
    refined.push_back((*iz)->info["REFINED"]);
    
    if((*iz)->collapsed > 0){
      COLLAPSED += (*iz)->collapsed; 
    }
    else{
      COLLAPSED += 1;
    }
    
  }
  
  if(COLLAPSED == 1){
    COLLAPSED = 0;
  }
  
  double fiveRaw = mean(FIVEPOS);
  double threeRaw = mean(THREEPOS);

  double fiveSD  = sqrt(var(FIVEPOS,   fiveRaw));
  double threeSD = sqrt(var(THREEPOS, threeRaw));
    
  int  fiveAvg = roundHalfwayDown(fiveRaw);
  int  threeAvg = roundHalfwayDown(threeRaw);
  
  int fiveCIL = -roundHalfwayDown(fiveSD*2);
  int fiveCIH = -fiveCIL;

  if(fiveCIL > -10){
    fiveCIL = -10;
  }
  if(fiveCIH < 10){
    fiveCIH = 10;
  }

  int threeCIL = -roundHalfwayDown(threeSD*2);
  int threeCIH = -threeCIL;

  if(threeCIL > -10){
    threeCIL = -10;
  }
  if(threeCIH < 10){
    threeCIH = 10;
  }


  
  int svlen = threeAvg - fiveAvg;
  if(svs.front()->type.compare("DEL")==0){
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
  
  string refBase = svs.front()->refBase;
  
  if(COLLAPSED > 0){
    refBase = "N";
  }

  bedpe << svs.front()->seqid << "\t"
	<< (fiveAvg) << "\t"
	<< "WG:" << svs.front()->type << ":" << newid << "\t"
	<< refBase << "\t"
	<< "<" << svs.front()->type << ">\t"
	<< ".\t.\t"
	<< "SVTYPE=" << svs.front()->type 
	<< ";SVLEN=" << svlen << ";ID=" 
	<< join(ids, ",") << ";SUPPORT=" 
	<< lSupport << "," << rSupport 
	<< ";MERGED=" << join(merged, ",") << ";REFINED=" << join(refined, ",") 
        << ";END="   << threeAvg 
	<< ";POS="   << fiveAvg << "," << threeAvg << ";LID=" << join(lidV, ",") 
        << ";FIVE="  << join(FIVEPOS, ",")
        << ";THREE=" << join(THREEPOS, ",")
	<< ";RID=" << join(ridV, ",") << ";" 
	<< "CIPOS=" << fiveCIL << "," << fiveCIH << ";"
        << "CIEND=" << threeCIL << "," << threeCIH << ";"
	<< "COLLAPSED=" << COLLAPSED 
	<< "\tGT:GL:AS:RS\t.:.:.:." ;

  // merged SVs need to be re genotyped;
  

    
  cout << bedpe.str() << endl;

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
  vector<string> info = split(sv->lineDat[7], ";");
  sv->refBase = sv->lineDat[3];
  for(vector<string>::iterator iz = info.begin(); iz != info.end(); iz++){
    vector<string> key_value = split(*iz, "=");
    sv->info[ key_value[0] ] = key_value[1];
  }
      
  for(int gi = 9; gi < sv->lineDat.size() ; gi++){
    sv->genotypeDat.push_back(sv->lineDat[gi]);
  }

  vector<string> pos = split(sv->info["POS"], ',');

  vector<string> support = split(sv->info["SUPPORT"], ',');

  sv->type = sv->info["SVTYPE"];
  sv->id   = sv->info["ID"];
    
  sv->collapsed = atoi(sv->info["COLLAPSED"].c_str());
  sv->ls = atoi(support.front().c_str());
  sv->rs = atoi(support.back().c_str());

  sv->lid = split(sv->info["LID"], ',');
  sv->rid = split(sv->info["RID"], ',');
  
  sv->five = atoi(pos[0].c_str());
  sv->three = atoi(pos[1].c_str());
  sv->seqid = sv->lineDat[0];

  vector<std::string> fives;
  vector<std::string> threes;
  
  if(sv->info["FIVE"].find(",") != std::string::npos){
    fives = split(sv->info["FIVE"], ",");
  }
  else{
    fives.push_back(sv->info["FIVE"].c_str());
  }
  
  if(sv->info["THREE"].find(",") != std::string::npos){
    threes = split(sv->info["THREE"], ",");
  }
  else{
    threes.push_back(sv->info["THREE"].c_str());
  }
  for(int i = 0; i < fives.size(); i++){ 
    sv->fives.push_back(atoi(fives[i].c_str()));
    sv->threes.push_back(atoi(threes[i].c_str()));
  } 
}

//-------------------------------    MAIN     --------------------------------
/*
 Comments:
*/

int main( int argc, char** argv)
{
globalOpts.maxDist = 5;
parseOpts(argc, argv);

 vector<svDat *> SVbuffer;

 ifstream myfile (globalOpts.filename.c_str());
 
 string line;

 printHeader();


 if(myfile.is_open()){
   while( getline (myfile, line) ){
     
     if(line.substr(0,6).compare("#CHROM") == 0){
       vector<string> CHROM = split(line, "\t");
       NSAMP = CHROM.size() - 8;
     }

     if(line[0] == '#'){
       continue;
     }

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
   printHelp();
   cerr << endl;
   exit(1);
 }


return 0;


}
