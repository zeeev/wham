#ifndef KMERUTILS
#define KMERUTILS

#include <unistd.h>
#include <stdint.h>
#define KMER_LEN 17

const uint64_t leftMask = 3;
const uint64_t dnA = 0; // 00
const uint64_t dnT = 1; // 01
const uint64_t dnG = 2; // 10
const uint64_t dnC = 3; // 11
const char blank[] = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
      char tests[] = "TCGAGAGGGGCCCAAAATTTTCCCCGGGCCCC";

char dna_lookup[4] = {'A', 'T', 'G', 'C'};

// ascii printing characters

uint8_t  ascii_lookup[130] = { 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,   // 9
                               5, 5, 5, 5, 5, 5, 5, 5, 5, 5,   // 19
                               5, 5, 5, 5, 5, 5, 5, 5, 5, 5,   // 29
                               5, 5, 5, 5, 5, 5, 5, 5, 5, 5,   // 39
                               5, 5, 5, 5, 5, 5, 5, 5, 5, 5,   // 49
                               5, 5, 5, 5, 5, 5, 5, 5, 5, 5,   // 59
                               5, 5, 5, 5, 5, 0, 5, 3, 5, 5,   // 69
                               5, 2, 5, 5, 5, 5, 5, 5, 5, 5,   // 79
                               5, 5, 5, 5, 1, 5, 5, 5, 5, 5,   // 89
                               5, 5, 5, 5, 5, 5, 5, 0, 5, 3,   // 99
                               5, 5, 5, 2, 5, 5, 5, 5, 5, 5,   // 109
                               5, 5, 5, 5, 5, 5, 1, 5, 5, 5,   // 119
                               5, 5, 5, 5, 5, 5, 5, 5, 5, 5 }; // 129
      

uint8_t BinToChar(const uint64_t bin, char * dna){

  uint64_t modBin = bin;

  for(uint8_t i = KMER_LEN; i > 0; i--){
    //    printf("translating bit: %d\t", (int)(leftMask & modBin));
    //    printf("tranlated bit: %c\n", dna_lookup[(leftMask & modBin)]);
    dna[i-1] = dna_lookup[(leftMask & modBin)];
    modBin >>= 2;
  }
  return 1;
}

uint64_t charArrayToBin(char * trace, uint32_t offset){
  uint64_t bin = 0;
  for(uint32_t i = offset ; i < (offset + KMER_LEN) ; i++){
//    printf("string length %d\t", (int)strlen(trace));
//    printf("i %d\t", i);
//    printf("offset %d\n", offset);

    if(dnA == ascii_lookup[(int)trace[i]]){
      bin |= dnA;
//      printf("base: %c\t", trace[i]);
//      printf("bin: %i\n", (int)bin);
    }
    else if (dnT == ascii_lookup[(int)trace[i]]){
      bin |= dnT;
      //      printf("base: %c\t", trace[i]);
      //      printf("bin: %i\n", (int)bin);
    }
    else if (dnG == ascii_lookup[(int)trace[i]]){
      bin |= dnG;
      //     printf("base: %c\t", trace[i]);
      //     printf("bin: %i\n", (int)bin);
    }
    else if (dnC == ascii_lookup[(int)trace[i]]){
      bin |= dnC;
      //      printf("base: %c\t", trace[i]);
      //      printf("bin: %i\n", (int)bin);
    }
    if(i < (offset + KMER_LEN -1)){
      bin <<= 2;
    }
  }
  return bin;
}

#endif
