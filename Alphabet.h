//----------------------------------------------------------------
#ifndef AlphabetH
#define AlphabetH
//----------------------------------------------------------------
#include <stdint.h>
//----------------------------------------------------------------
#define A           1
#define C           2
#define G           4
#define T           8
#define N           15
#define INS         16
#define DEL         32
//----------------------------------------------------------------
const int iAscii2BamSeq[256]=
{
    N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,
    N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,
    N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,
    N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,
    N,A,N,C,N,N,N,G,N,N,N,N,N,N,N,N,
    N,N,N,N,T,N,N,N,N,N,N,N,N,N,N,N,
    N,A,N,C,N,N,N,G,N,N,N,N,N,N,N,N,
    N,N,N,N,T,N,N,N,N,N,N,N,N,N,N,N,
    N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,
    N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,
    N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,
    N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,
    N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,
    N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,
    N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,
    N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N    
};
//----------------------------------------------------------------
#endif

