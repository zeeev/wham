#include  "api/api_global.h"
#include  "api/BamReader.h"
#include <iostream>
#include <fstream>

using namespace std;
using namespace BamTools;

int main(int argc,  char * argv[]){

  if(argc > 2){
    cerr << "FATAL: too many options specified\n";
    cerr << "INFO: usage : anti-align myPaired-end.bam";
    return(1);
  }

  BamReader br;
  string fh = argv[1];

  if(! br.Open(fh)){
    cerr << "FATAL: could not open bam\n";
    cerr << "INFO: usage : anti-align myPaired-end.bam";
    return(1);
  }


  ofstream p1, p2;

  p1.open("PE1.fastq");
  p2.open("PE2.fastq");

  BamAlignment al;
  map<string, BamAlignment> reads;

  while(br.GetNextAlignment(al)){
    
    string rawBases  = al.QueryBases;
    string readName  = al.Name;
    
    if(al.IsMapped()){
      reads.erase(readName);
      continue;
    }

    // does the fastq contain Ns
    int    containsN = 0;
    for(int i = 0; i < rawBases.length(); i++){
      if (rawBases.at(i) == 'N'){
	containsN = 1;
	break;
      }
    }
    
    if(containsN == 1){
      reads.erase(readName);
      continue; 
    }
  
    map<string, BamAlignment>::iterator it;
    it = reads.find(readName);

    if(it != reads.end()){
      if(al.IsFirstMate()){
	p1 << "@" << readName << "/1"     << endl;
	p1 << rawBases                    << endl;
	p1 << "+"                         << endl;
	p1 << al.Qualities                << endl;
	p2 << "@" << reads[readName].Name << "/2" << endl;
	p2 << reads[readName].QueryBases  << endl;
	p2 << "+"                         << endl;
       	p2 << reads[readName].Qualities   << endl;
      }
      else{
	p2 << "@" << readName << "/2"      << endl;
        p2 << rawBases                     << endl;
	p2 << "+"                          << endl;
        p2 << al.Qualities                 << endl;
	p1 << "@" << reads[readName].Name  << "/1" << endl;
	p1 << reads[readName].QueryBases   << endl;
        p1 << "+"                          << endl;
	p1 << reads[readName].Qualities    << endl;
      }
      reads.erase(readName);
    }
    else{
      reads[readName] = al;
    }
  }
    p1.close();
    p2.close();
    br.Close();
    
    cerr << "INFO: antialign finished without errors" << endl;
    return 0;
}


