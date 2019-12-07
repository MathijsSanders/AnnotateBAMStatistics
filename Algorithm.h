//----------------------------------------------------------------
#ifndef AlgorithmH
#define AlgorithmH
//----------------------------------------------------------------
#include <string>
#include <vector>
#include <deque>
#include "AnnovarFile.h"
#include "BamFile.h"
#include "Mutex.h"
//----------------------------------------------------------------
using namespace std;
//----------------------------------------------------------------
typedef struct
{
    int     iBamFile;
    int     iPos;
    string  sTarget;
    
}JOB;
//----------------------------------------------------------------
typedef struct
{
    int     iSecDepth;
    int     iSecAltDepth;
    float   fSecAltBias;
    int     iSecHQDepth;
    int     iSecHQAltDepth;
    float   fSecHQAltBias;
    int     iDepth;
    int     iAltDepth;
    float   fAltBias;
    int     iHQDepth;
    int     iHQAltDepth;
    float   fHQAltBias;
    int     iXADepth;
    int     iXAAltDepth;

}STATISTICS;
//----------------------------------------------------------------
typedef struct
{
    int                 iAlt;
    int                 iLen;
    int                 iMinBaseScore;
    int                 iMinAlignmentScore;    
    int                 iMinMatchLength;

    PLP_DATA *          pPileupData;
    STATISTICS *        pStatistics;

}COMPUTE_STATISTICS;
//----------------------------------------------------------------
class CAlgorithm
{
    private:
        
        int             nJobsDone;

        CAnnovarFile    annovarFile;
        CMutex          mutex;
        
        deque<JOB>      jobs;
        
        static string & SampleNameFromPath(const string & sPath);

        static void     ComputeStatisticsSNP(COMPUTE_STATISTICS & pComputeStatistics);
        static void     ComputeStatisticsIndel(COMPUTE_STATISTICS & pComputeStatistics);
        static void *   Pileup(CAlgorithm * pAlgorithm);
        static void *   PileupRegions(CAlgorithm * pAlgorithm);
        
        void CreateOutputMatrix(void);
        void CreateJobs(void);
        void RunJobs(void);
        void WriteOutputMatrix(void);
        
    public:
        
        bool            bCountDuplicates;
        bool            bCountSupplementary;
        bool            bPileupRegions;
        bool            bVerbose;
        
        int             iMinBaseScore;
        int             iMinAlignmentScore;
        int             iMinMatchLength;
        int             nThreads;

        string          sAnnovarFileName;
        vector<string>  vBamFileNames;
        
        CAlgorithm(void);
        bool Run(void);
};
//----------------------------------------------------------------
#endif
