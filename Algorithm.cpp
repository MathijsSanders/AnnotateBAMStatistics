//----------------------------------------------------------------
// Name        : Algorithm.cpp
// Author      : Remco Hoogenboezem
// Version     :
// Copyright   :
// Description :
//----------------------------------------------------------------
#include <iostream>
#include <string>
#include <unordered_map>
#include <string.h>
#include <unistd.h>
#include "Alphabet.h"
#include "ProgressBar.h"
#include "Algorithm.h"
//----------------------------------------------------------------
string & CAlgorithm::SampleNameFromPath(const string & sPath)
{
    size_t          lastSlashPos;
    size_t          lastPeriodPos;
    static string   sSampleName;
    
    sSampleName=sPath;
    
    lastSlashPos=sSampleName.find_last_of("\\/");

    if(lastSlashPos!=string::npos)
    {
        sSampleName.erase(0,lastSlashPos+1);
    }

    lastPeriodPos=sSampleName.rfind('.');

    if(lastPeriodPos!=string::npos)
    {
        sSampleName.erase(lastPeriodPos);
    }
    
    return sSampleName;
}
//----------------------------------------------------------------
void CAlgorithm::ComputeStatisticsSNP(COMPUTE_STATISTICS & computeStatistics)
{
    bool    bIsAlt;
    bool    bIsHQ;

    int     nAlignments;
    int     iAlignment;
    int     nOperations;
    int     iOperation;
    int     iMatchLength;
    int     iMinMatchLength;
    int     iBase;
    int     iAlt;
    int     iMinBaseScore;
    int     iMinAlignmentScore;
    int     iQueryPos;
    int     iBaseScore;
    int     iStrand;
    int     iSecAltReadDepth[2];
    int     iSecHQAltReadDepth[2];
    int     iAltReadDepth[2];
    int     iHQAltReadDepth[2];
    int     iXADepth;
    int     iXAAltDepth;
    int     iSecDepth;
    int     iSecAltDepth;
    int     iSecHQDepth;
    int     iSecHQAltDepth;
    int     iDepth;
    int     iAltDepth;
    int     iHQDepth;
    int     iHQAltDepth;
    
    const bam_pileup1_t *   pAlignments;
    const bam1_t *          pAlignments_b;
    const bam_pileup1_t *   pMate;
    STATISTICS *            pStatistics;
    uint32_t *              pCIGAR;

    string                  sQueryName;

    unordered_map<string,const bam_pileup1_t *>             mAlignments;
    unordered_map<string,const bam_pileup1_t *>::iterator   itAlignments;
    unordered_map<string,const bam_pileup1_t *>::iterator   itAlignmentsEnd;
    
    //----------------------------------------------------------------
    //Initialize
    //----------------------------------------------------------------
    
    nAlignments=computeStatistics.pPileupData->nAlignments;
    pAlignments=computeStatistics.pPileupData->pAlignments;
    
    iAlt=computeStatistics.iAlt;
    iMinBaseScore=computeStatistics.iMinBaseScore;
    iMinAlignmentScore=computeStatistics.iMinAlignmentScore;
    iMinMatchLength=computeStatistics.iMinMatchLength;
    
    iSecAltReadDepth[0]=0; 
    iSecAltReadDepth[1]=0;
    iSecHQAltReadDepth[0]=0;
    iSecHQAltReadDepth[1]=0;
    iAltReadDepth[0]=0;
    iAltReadDepth[1]=0;
    iHQAltReadDepth[0]=0;
    iHQAltReadDepth[1]=0;
    
    iXADepth=0;
    iXAAltDepth=0;

    iSecDepth=0;
    iSecAltDepth=0;
    iSecHQDepth=0;
    iSecHQAltDepth=0;
    iDepth=0;
    iAltDepth=0;
    iHQDepth=0;
    iHQAltDepth=0;
    
    //----------------------------------------------------------------
    //Iterate over reads
    //----------------------------------------------------------------
    
    for(iAlignment=0;iAlignment<nAlignments;iAlignment++,pAlignments++)
    {
        //----------------------------------------------------------------
        //Skip del and refskip reads
        //----------------------------------------------------------------
        
        if((pAlignments->is_del)||(pAlignments->is_refskip))
        {
            continue;
        }

        //----------------------------------------------------------------
        //Skip if matchlength is too small
        //----------------------------------------------------------------

        pAlignments_b=pAlignments->b;
        
        nOperations=pAlignments_b->core.n_cigar;
        iMatchLength=0;
            
        pCIGAR=bam1_cigar(pAlignments_b);
            
        for(iOperation=0;iOperation<nOperations;iOperation++,pCIGAR++)
        {
            iMatchLength+=(bam_cigar_op(*pCIGAR)==BAM_CMATCH)*bam_cigar_oplen(*pCIGAR);
        }
            
        if(iMatchLength<iMinMatchLength)
        {
            continue;
        }
        
        //----------------------------------------------------------------
        //Skip if base is ambiguous
        //----------------------------------------------------------------
        
        iQueryPos=pAlignments->qpos;
        iBase=int(bam1_seqi(bam1_seq(pAlignments_b),iQueryPos));
        
        if(iBase==N)
        {
            continue;
        }
        
        //----------------------------------------------------------------
        //Compute read statistics
        //----------------------------------------------------------------
        
        bIsAlt=iBase==iAlt;
        bIsHQ=((iBaseScore=bam1_qual(pAlignments_b)[iQueryPos])>=iMinBaseScore) && (pAlignments_b->core.qual>=iMinAlignmentScore);
        
        iStrand=int(bam1_strand(pAlignments_b));
        
        iSecAltReadDepth[iStrand]+=bIsAlt;
        
        if(bIsHQ)
        {
            iSecHQAltReadDepth[iStrand]+=bIsAlt;
        }
        
        if((pAlignments_b->core.flag&BAM_FSECONDARY)==0)
        {
            iAltReadDepth[iStrand]+=bIsAlt;
            
            if(bIsHQ)
            {
                iHQAltReadDepth[iStrand]+=bIsAlt;
            }
        }
       
        if(bam_aux_get(pAlignments_b,"XA")!=NULL)
        {
            iXADepth++;
            iXAAltDepth+=bIsAlt;
        }
        
        //----------------------------------------------------------------
        //Accumulate fragments
        //----------------------------------------------------------------
        
        sQueryName=string(bam1_qname(pAlignments_b));
        itAlignments=mAlignments.find(sQueryName);
        
        if(itAlignments==mAlignments.end())
        {
            mAlignments[sQueryName]=pAlignments;
        }
        
        else
        {
            pMate=itAlignments->second;
            
            if(bam1_qual(pMate->b)[pMate->qpos]<iBaseScore) //Keep the mate with the best base score 
            {
                mAlignments[sQueryName]=pAlignments;
            }
        }
    }
    
    //----------------------------------------------------------------
    //Iterate over fragments
    //----------------------------------------------------------------
    
    itAlignmentsEnd=mAlignments.end();
    
    for(itAlignments=mAlignments.begin();itAlignments!=itAlignmentsEnd;itAlignments++)
    {
        //----------------------------------------------------------------
        //Compute fragment statistics
        //----------------------------------------------------------------

        iQueryPos=itAlignments->second->qpos;
        pAlignments_b=itAlignments->second->b;

        bIsAlt=int(bam1_seqi(bam1_seq(pAlignments_b),iQueryPos))==iAlt;        
        bIsHQ=(bam1_qual(pAlignments_b)[iQueryPos]>=iMinBaseScore) && (pAlignments_b->core.qual>=iMinAlignmentScore);
        
        iSecDepth++;
        iSecAltDepth+=bIsAlt;
        
        if(bIsHQ)
        {
            iSecHQDepth++;
            iSecHQAltDepth+=bIsAlt;
        }
        
        if((pAlignments_b->core.flag&BAM_FSECONDARY)==0)
        {
            iDepth++;
            iAltDepth+=bIsAlt;
            
            if(bIsHQ)
            {
                iHQDepth++;
                iHQAltDepth+=bIsAlt;
            }
        }
    }
    
    //----------------------------------------------------------------
    //Save results
    //----------------------------------------------------------------
    
    pStatistics=computeStatistics.pStatistics;
    
    pStatistics->iSecDepth=iSecDepth;
    pStatistics->iSecAltDepth=iSecAltDepth;
    pStatistics->fSecAltBias=float(double(iSecAltReadDepth[0])/double(max(iSecAltReadDepth[0]+iSecAltReadDepth[1],1)));
    pStatistics->iSecHQDepth=iSecHQDepth;
    pStatistics->iSecHQAltDepth=iSecHQAltDepth;
    pStatistics->fSecHQAltBias=float(double(iSecHQAltReadDepth[0])/double(max(iSecHQAltReadDepth[0]+iSecHQAltReadDepth[1],1)));
    pStatistics->iDepth=iDepth;
    pStatistics->iAltDepth=iAltDepth;
    pStatistics->fAltBias=float(double(iAltReadDepth[0])/double(max(iAltReadDepth[0]+iAltReadDepth[1],1)));
    pStatistics->iHQDepth=iHQDepth;
    pStatistics->iHQAltDepth=iHQAltDepth;
    pStatistics->fHQAltBias=float(double(iHQAltReadDepth[0])/double(max(iHQAltReadDepth[0]+iHQAltReadDepth[1],1)));
    pStatistics->iXADepth=iXADepth;
    pStatistics->iXAAltDepth=iXAAltDepth;
    
    //----------------------------------------------------------------
    //Done
    //----------------------------------------------------------------
}
//----------------------------------------------------------------
void CAlgorithm::ComputeStatisticsIndel(COMPUTE_STATISTICS & computeStatistics)
{
    bool    bIsAlt;
    bool    bIsHQ;
    
    int     nAlignments;
    int     iAlignment;
    int     nOperations;
    int     iOperation;
    int     iMatchLength;
    int     iMinMatchLength;
    int     iLen;
    int     iMinAlignmentScore;
    int     iAlignmentScore;
    int     iStrand;
    int     iSecAltReadDepth[2];
    int     iSecHQAltReadDepth[2];
    int     iAltReadDepth[2];
    int     iHQAltReadDepth[2];
    int     iXADepth;
    int     iXAAltDepth;
    int     iSecDepth;
    int     iSecAltDepth;
    int     iSecHQDepth;
    int     iSecHQAltDepth;
    int     iDepth;
    int     iAltDepth;
    int     iHQDepth;
    int     iHQAltDepth;
    
    const bam_pileup1_t *   pAlignments;
    const bam1_t *          pAlignments_b;
    const bam_pileup1_t *   pMate;
    STATISTICS *            pStatistics;
    uint32_t *              pCIGAR;

    string                  sQueryName;

    unordered_map<string,const bam_pileup1_t *>             mAlignments;
    unordered_map<string,const bam_pileup1_t *>::iterator   itAlignments;
    unordered_map<string,const bam_pileup1_t *>::iterator   itAlignmentsEnd;
    
    //----------------------------------------------------------------
    //Initialize
    //----------------------------------------------------------------
    
    nAlignments=computeStatistics.pPileupData->nAlignments;
    pAlignments=computeStatistics.pPileupData->pAlignments;
    
    iLen=computeStatistics.iLen;
    iMinAlignmentScore=computeStatistics.iMinAlignmentScore;
    iMinMatchLength=computeStatistics.iMinMatchLength;
    
    iSecAltReadDepth[0]=0; 
    iSecAltReadDepth[1]=0;
    iSecHQAltReadDepth[0]=0;
    iSecHQAltReadDepth[1]=0;
    iAltReadDepth[0]=0;
    iAltReadDepth[1]=0;
    iHQAltReadDepth[0]=0;
    iHQAltReadDepth[1]=0;
    
    iXADepth=0;
    iXAAltDepth=0;

    iSecDepth=0;
    iSecAltDepth=0;
    iSecHQDepth=0;
    iSecHQAltDepth=0;
    iDepth=0;
    iAltDepth=0;
    iHQDepth=0;
    iHQAltDepth=0;
    
    //----------------------------------------------------------------
    //Iterate over reads
    //----------------------------------------------------------------
    
    for(iAlignment=0;iAlignment<nAlignments;iAlignment++,pAlignments++)
    {
        //----------------------------------------------------------------
        //Skip del and refskip reads
        //----------------------------------------------------------------
        
        if((pAlignments->is_del)||(pAlignments->is_refskip))
        {
            continue;
        }

        //----------------------------------------------------------------
        //Skip if matchlength is too small
        //----------------------------------------------------------------

        pAlignments_b=pAlignments->b;
        
        nOperations=pAlignments_b->core.n_cigar;
        iMatchLength=0;
            
        pCIGAR=bam1_cigar(pAlignments_b);
            
        for(iOperation=0;iOperation<nOperations;iOperation++,pCIGAR++)
        {
            iMatchLength+=(bam_cigar_op(*pCIGAR)==BAM_CMATCH)*bam_cigar_oplen(*pCIGAR);
        }
            
        if(iMatchLength<iMinMatchLength)
        {
            continue;
        }
        
        //----------------------------------------------------------------
        //Compute read statistics
        //----------------------------------------------------------------
        
        bIsAlt=pAlignments->indel==iLen;
        bIsHQ=(iAlignmentScore=pAlignments_b->core.qual)>=iMinAlignmentScore;
        
        iStrand=int(bam1_strand(pAlignments_b));
        
        iSecAltReadDepth[iStrand]+=bIsAlt;
        
        if(bIsHQ)
        {
            iSecHQAltReadDepth[iStrand]+=bIsAlt;
        }
        
        if((pAlignments_b->core.flag&BAM_FSECONDARY)==0)
        {
            iAltReadDepth[iStrand]+=bIsAlt;
            
            if(bIsHQ)
            {
                iHQAltReadDepth[iStrand]+=bIsAlt;
            }
        }
       
        if(bam_aux_get(pAlignments_b,"XA")!=NULL)
        {
            iXADepth++;
            iXAAltDepth+=bIsAlt;
        }
        
        //----------------------------------------------------------------
        //Accumulate fragments
        //----------------------------------------------------------------
        
        sQueryName=string(bam1_qname(pAlignments_b));
        itAlignments=mAlignments.find(sQueryName);
        
        if(itAlignments==mAlignments.end())
        {
            mAlignments[sQueryName]=pAlignments;
        }
        
        else
        {
            pMate=itAlignments->second;
            
            if((pMate->indel==iLen)==bIsAlt)        
            {
                if(pMate->b->core.qual<iAlignmentScore)
                {
                    mAlignments[sQueryName]=pAlignments;
                }
            }
            
            else if(bIsAlt)
            {
                mAlignments[sQueryName]=pAlignments;
            }
        }
    }
    
    //----------------------------------------------------------------
    //Iterate over fragments
    //----------------------------------------------------------------
    
    itAlignmentsEnd=mAlignments.end();
    
    for(itAlignments=mAlignments.begin();itAlignments!=itAlignmentsEnd;itAlignments++)
    {
        //----------------------------------------------------------------
        //Compute fragment statistics
        //----------------------------------------------------------------

        pAlignments_b=itAlignments->second->b;

        bIsAlt=itAlignments->second->indel==iLen;
        bIsHQ=pAlignments_b->core.qual>=iMinAlignmentScore;
        
        iSecDepth++;
        iSecAltDepth+=bIsAlt;
        
        if(bIsHQ)
        {
            iSecHQDepth++;
            iSecHQAltDepth+=bIsAlt;
        }
        
        if((pAlignments_b->core.flag&BAM_FSECONDARY)==0)
        {
            iDepth++;
            iAltDepth+=bIsAlt;
            
            if(bIsHQ)
            {
                iHQDepth++;
                iHQAltDepth+=bIsAlt;
            }
        }
    }
    
    //----------------------------------------------------------------
    //Save results
    //----------------------------------------------------------------
    
    pStatistics=computeStatistics.pStatistics;
    
    pStatistics->iSecDepth=iSecDepth;
    pStatistics->iSecAltDepth=iSecAltDepth;
    pStatistics->fSecAltBias=float(double(iSecAltReadDepth[0])/double(max(iSecAltReadDepth[0]+iSecAltReadDepth[1],1)));
    pStatistics->iSecHQDepth=iSecHQDepth;
    pStatistics->iSecHQAltDepth=iSecHQAltDepth;
    pStatistics->fSecHQAltBias=float(double(iSecHQAltReadDepth[0])/double(max(iSecHQAltReadDepth[0]+iSecHQAltReadDepth[1],1)));
    pStatistics->iDepth=iDepth;
    pStatistics->iAltDepth=iAltDepth;
    pStatistics->fAltBias=float(double(iAltReadDepth[0])/double(max(iAltReadDepth[0]+iAltReadDepth[1],1)));
    pStatistics->iHQDepth=iHQDepth;
    pStatistics->iHQAltDepth=iHQAltDepth;
    pStatistics->fHQAltBias=float(double(iHQAltReadDepth[0])/double(max(iHQAltReadDepth[0]+iHQAltReadDepth[1],1)));
    pStatistics->iXADepth=iXADepth;
    pStatistics->iXAAltDepth=iXAAltDepth;
    
    //----------------------------------------------------------------
    //Done
    //----------------------------------------------------------------    
}
//----------------------------------------------------------------
void * CAlgorithm::Pileup(CAlgorithm * pAlgorithm)
{
    multimap<int,ENTRY> *   pAnnovarFileEntries;

    JOB                     job;
    COMPUTE_STATISTICS      computeStatistics;

    CBamFile                bamFile;
    
    pair<multimap<int,ENTRY>::iterator,multimap<int,ENTRY>::iterator>   equalRange;
    multimap<int,ENTRY>::iterator                                       itEqualRange;
    
    computeStatistics.iMinBaseScore=pAlgorithm->iMinBaseScore;
    computeStatistics.iMinAlignmentScore=pAlgorithm->iMinAlignmentScore;
    computeStatistics.iMinMatchLength=pAlgorithm->iMinMatchLength;

    pAlgorithm->mutex.Lock();
    
    for(;;)
    {
        //----------------------------------------------------------------
        //Get job from queue
        //----------------------------------------------------------------
        
        if(pAlgorithm->jobs.size()==0)
        {
            pAlgorithm->mutex.Unlock();
            pthread_exit(NULL);
        }

        job=pAlgorithm->jobs.front();
        pAlgorithm->jobs.pop_front();

        pAlgorithm->mutex.Unlock();

        //----------------------------------------------------------------
        //Open bam file
        //----------------------------------------------------------------
        
        if(bamFile.Open(pAlgorithm->vBamFileNames[job.iBamFile],0,pAlgorithm->bCountDuplicates,pAlgorithm->bCountSupplementary)==false)
        {
            cerr << "Error: Could not open bam file: " << pAlgorithm->vBamFileNames[job.iBamFile] << " (skipping job)" << endl;
            pAlgorithm->mutex.Lock();
            pAlgorithm->nJobsDone++;
            continue;
        }
        
        //----------------------------------------------------------------
        //Set target on bam file
        //----------------------------------------------------------------
        
        if(bamFile.SetRegion(job.sTarget)==false)
        {
            cerr << "Error: Could not set region: " << job.sTarget << " (skipping job)" << endl;
            pAlgorithm->mutex.Lock();            
            pAlgorithm->nJobsDone++;
            continue;
        }

        computeStatistics.pPileupData=bamFile.GetPileupData();
        
        //----------------------------------------------------------------
        //Get pointer to annovar file
        //----------------------------------------------------------------
        
        pAnnovarFileEntries=&pAlgorithm->annovarFile.mEntries.at(job.sTarget);
        
        //----------------------------------------------------------------
        //Pileup
        //----------------------------------------------------------------
        
        while(bamFile.PileupRegion())
        {
            equalRange=pAnnovarFileEntries->equal_range(computeStatistics.pPileupData->iTargetPos);
            
            for(itEqualRange=equalRange.first;itEqualRange!=equalRange.second;itEqualRange++)
            {
                computeStatistics.iAlt=itEqualRange->second.iAlt;
                computeStatistics.iLen=itEqualRange->second.iLen;
                computeStatistics.pStatistics=((STATISTICS*)itEqualRange->second.pUserData)+job.iBamFile;

                if((computeStatistics.iAlt!=INS)&&(computeStatistics.iAlt!=DEL))
                {
                    ComputeStatisticsSNP(computeStatistics);
                }
                
                else
                {
                    ComputeStatisticsIndel(computeStatistics);
                }
            }
        }
        
        //----------------------------------------------------------------
        //Job done
        //----------------------------------------------------------------
        
        pAlgorithm->mutex.Lock();
        pAlgorithm->nJobsDone++;
    }
}
//----------------------------------------------------------------
void * CAlgorithm::PileupRegions(CAlgorithm * pAlgorithm)
{
    multimap<int,ENTRY> *   pAnnovarFileEntries;

    JOB                     job;
    COMPUTE_STATISTICS      computeStatistics;

    CBamFile                bamFile;
    
    string                  sRegion;
    
    pair<multimap<int,ENTRY>::iterator,multimap<int,ENTRY>::iterator>   equalRange;
    multimap<int,ENTRY>::iterator                                       itEqualRange;
    
    computeStatistics.iMinBaseScore=pAlgorithm->iMinBaseScore;
    computeStatistics.iMinAlignmentScore=pAlgorithm->iMinAlignmentScore;
    computeStatistics.iMinMatchLength=pAlgorithm->iMinMatchLength;
    
    pAlgorithm->mutex.Lock();
    
    for(;;)
    {
        //----------------------------------------------------------------
        //Get job from queue
        //----------------------------------------------------------------
        
        if(pAlgorithm->jobs.size()==0)
        {
            pAlgorithm->mutex.Unlock();
            pthread_exit(NULL);
        }

        job=pAlgorithm->jobs.front();
        pAlgorithm->jobs.pop_front();

        pAlgorithm->mutex.Unlock();

        //----------------------------------------------------------------
        //Open bam file
        //----------------------------------------------------------------
        
        if(bamFile.Open(pAlgorithm->vBamFileNames[job.iBamFile],0,pAlgorithm->bCountDuplicates,pAlgorithm->bCountSupplementary)==false)
        {
            cerr << "Error: Could not open bam file: " << pAlgorithm->vBamFileNames[job.iBamFile] << " (skipping job)" << endl;
            pAlgorithm->mutex.Lock();
            pAlgorithm->nJobsDone++;
            continue;
        }
        
        //----------------------------------------------------------------
        //Set target on bam file
        //----------------------------------------------------------------
        
        sRegion=job.sTarget+string(":")+to_string(max(job.iPos-10,1))+string("-")+to_string(job.iPos+10);
        
        if(bamFile.SetRegion(sRegion)==false)
        {
            cerr << "Error: Could not set region: " << sRegion << " (skipping job)" << endl;
            pAlgorithm->mutex.Lock();
            pAlgorithm->nJobsDone++;
            continue;
        }

        computeStatistics.pPileupData=bamFile.GetPileupData();
        
        //----------------------------------------------------------------
        //Get pointer to annovar file
        //----------------------------------------------------------------
        
        pAnnovarFileEntries=&pAlgorithm->annovarFile.mEntries.at(job.sTarget);
        
        //----------------------------------------------------------------
        //Pileup
        //----------------------------------------------------------------
        
        while(bamFile.PileupRegion())
        {
            if(computeStatistics.pPileupData->iTargetPos==job.iPos) 
            {
                equalRange=pAnnovarFileEntries->equal_range(computeStatistics.pPileupData->iTargetPos);

                for(itEqualRange=equalRange.first;itEqualRange!=equalRange.second;itEqualRange++)
                {
                    computeStatistics.iAlt=itEqualRange->second.iAlt;
                    computeStatistics.iLen=itEqualRange->second.iLen;
                    computeStatistics.pStatistics=((STATISTICS*)itEqualRange->second.pUserData)+job.iBamFile;

                    if((computeStatistics.iAlt!=INS)&&(computeStatistics.iAlt!=DEL))
                    {
                        ComputeStatisticsSNP(computeStatistics);
                    }
                
                    else
                    {
                        ComputeStatisticsIndel(computeStatistics);
                    }
                }
            }
        }
        
        //----------------------------------------------------------------
        //Job done
        //----------------------------------------------------------------
        
        pAlgorithm->mutex.Lock();
        pAlgorithm->nJobsDone++;
    }
}
//----------------------------------------------------------------
void CAlgorithm::CreateOutputMatrix(void)
{
    int nEntries;
    int nBamFiles;

    STATISTICS * pStatistics;    
    
    map<string,multimap<int,ENTRY> >::iterator  itTarget;
    map<string,multimap<int,ENTRY> >::iterator  itTargetEnd;
    multimap<int,ENTRY>::iterator   itPos;
    multimap<int,ENTRY>::iterator   itPosEnd;
    
    nEntries=annovarFile.nEntries;
    nBamFiles=vBamFileNames.size();
    
    pStatistics=(STATISTICS*)memset(new STATISTICS[nEntries*nBamFiles],0,nEntries*nBamFiles*sizeof(STATISTICS));

    itTargetEnd=annovarFile.mEntries.end();
    
    for(itTarget=annovarFile.mEntries.begin();itTarget!=itTargetEnd;itTarget++)
    {
        itPosEnd=itTarget->second.end();
        
        for(itPos=itTarget->second.begin();itPos!=itPosEnd;itPos++,pStatistics+=nBamFiles)
        {
            itPos->second.pUserData=pStatistics;
        }
    }
}
//----------------------------------------------------------------
void CAlgorithm::CreateJobs(void)
{
    int nBamFiles;
    JOB job;
    
    map<string,multimap<int,ENTRY> >::iterator  itTarget;
    map<string,multimap<int,ENTRY> >::iterator  itTargetEnd;
    multimap<int,ENTRY>::iterator   itPos;
    multimap<int,ENTRY>::iterator   itPosEnd;

    nBamFiles=vBamFileNames.size();
    
    itTargetEnd=annovarFile.mEntries.end();
    
    if(bPileupRegions)  
    {
        for(job.iBamFile=0;job.iBamFile<nBamFiles;job.iBamFile++)
        {
            for(itTarget=annovarFile.mEntries.begin();itTarget!=itTargetEnd;itTarget++)
            {
                job.iPos=-1;
                job.sTarget=itTarget->first;
                
                itPosEnd=itTarget->second.end();

                for(itPos=itTarget->second.begin();itPos!=itPosEnd;itPos++)
                {
                    if(job.iPos!=itPos->first)  //Don't add duplicate jobs from multimap
                    {
                        job.iPos=itPos->first;
                        jobs.push_back(job);
                    }
                }
            }
        }
    }
    
    else                        
    {
        job.iPos=-1;
        
        for(job.iBamFile=0;job.iBamFile<nBamFiles;job.iBamFile++)
        {
            for(itTarget=annovarFile.mEntries.begin();itTarget!=itTargetEnd;itTarget++)
            {
                job.sTarget=itTarget->first;
                jobs.push_back(job);
            }
        }
    }
}
//----------------------------------------------------------------
void CAlgorithm::RunJobs(void)
{
    int         nJobs;
    int         iThread;
    
    pthread_t * pThreads;
    void *      pThreadReturn;
    
    nJobs=jobs.size();
    nJobsDone=0;

    //----------------------------------------------------------------
    //Create worker threads
    //----------------------------------------------------------------
    
    if(bVerbose) cerr << "Info: Create worker threads" << endl;
    
    pThreads=new pthread_t[nThreads];
    
    if(bPileupRegions==false)
    {
        for(iThread=0;iThread<nThreads;iThread++)
        {
            pthread_create(&pThreads[iThread],NULL,(void*(*)(void*))Pileup,(void*)this);  
        }
    }
    
    else
    {
        for(iThread=0;iThread<nThreads;iThread++)
        {
            pthread_create(&pThreads[iThread],NULL,(void*(*)(void*))PileupRegions,(void*)this);  
        }
    }
    
    //----------------------------------------------------------------
    //Wait for threads to finish
    //----------------------------------------------------------------
    
    if(bVerbose) cerr << "Info: Wait for worker threads to finish" << endl;
    
    for(iThread=0;iThread<nThreads;iThread++)
    {
        while(pthread_tryjoin_np(pThreads[iThread],(void**)&pThreadReturn))
        {
            if(bVerbose) cerr << cProgressBar[min((nJobsDone*100)/nJobs,100)] << '\r' << flush;
            sleep(1);
        }
    }
    
    if(nJobsDone!=nJobs)
    {
        cerr << "Error: Number of jobs done differs from number of jobs" << endl;
    }
    
    if(bVerbose) cerr << cProgressBar[min((nJobsDone*100)/nJobs,100)]  << endl;
    

    delete [] pThreads; 
}
//----------------------------------------------------------------
void CAlgorithm::WriteOutputMatrix(void)
{
    int nBamFiles;
    int iBamFile;
    
    STATISTICS * pStatistics;
    
    map<string,multimap<int,ENTRY> >::iterator  itTarget;
    map<string,multimap<int,ENTRY> >::iterator  itTargetEnd;
    multimap<int,ENTRY>::iterator   itPos;
    multimap<int,ENTRY>::iterator   itPosEnd;
    
    //----------------------------------------------------------------
    //Header
    //----------------------------------------------------------------
    
    cout << annovarFile.sHeader;
    
    nBamFiles=vBamFileNames.size();
    
    for(iBamFile=0;iBamFile<nBamFiles;iBamFile++)
    {
//        cout << "\t" << SampleNameFromPath(vBamFileNames[iBamFile]); 

        cout    << "\t___"
                "\t(" << SampleNameFromPath(vBamFileNames[iBamFile]) << ")sec depth"
                "\tsec alt depth"
                "\tsec alt freq"
                "\tsec alt bias"
                "\tsec HQ depth"
                "\tsec HQ alt depth"
                "\tsec HQ alt freq"
                "\tsec HQ alt bias"
                "\tsec HQ ratio"
                "\tdepth"
                "\talt depth"
                "\talt freq"
                "\talt bias"
                "\tHQ depth"
                "\tHQ alt depth"
                "\tHQ alt freq"
                "\tHQ alt bias"
                "\tHQ ratio"
                "\tXA depth"
                "\tXA alt depth";

    }    
    
    cout << '\n';

    //----------------------------------------------------------------
    //Entries
    //----------------------------------------------------------------
    
    itTargetEnd=annovarFile.mEntries.end();
    
    for(itTarget=annovarFile.mEntries.begin();itTarget!=itTargetEnd;itTarget++)
    {
        itPosEnd=itTarget->second.end();
        
        for(itPos=itTarget->second.begin();itPos!=itPosEnd;itPos++)
        {
            cout << itPos->second.sLine;

            pStatistics=(STATISTICS*)itPos->second.pUserData;

            for(iBamFile=0;iBamFile<nBamFiles;iBamFile++,pStatistics++)
            {
//                cout << '\t' << float(double(pStatistics->iSecAltDepth)/double(max(pStatistics->iSecDepth,1)));

                cout    << '\t'
                        << '\t' << pStatistics->iSecDepth
                        << '\t' << pStatistics->iSecAltDepth
                        << '\t' << float(double(pStatistics->iSecAltDepth)/double(max(pStatistics->iSecDepth,1)))
                        << '\t' << pStatistics->fSecAltBias
                        << '\t' << pStatistics->iSecHQDepth
                        << '\t' << pStatistics->iSecHQAltDepth
                        << '\t' << float(double(pStatistics->iSecHQAltDepth)/double(max(pStatistics->iSecHQDepth,1)))
                        << '\t' << pStatistics->fSecHQAltBias
                        << '\t' << float(double(pStatistics->iSecHQDepth)/double(max(pStatistics->iSecDepth,1)))
                        << '\t' << pStatistics->iDepth
                        << '\t' << pStatistics->iAltDepth
                        << '\t' << float(double(pStatistics->iAltDepth)/double(max(pStatistics->iDepth,1)))
                        << '\t' << pStatistics->fAltBias
                        << '\t' << pStatistics->iHQDepth
                        << '\t' << pStatistics->iHQAltDepth
                        << '\t' << float(double(pStatistics->iHQAltDepth)/double(max(pStatistics->iHQDepth,1)))
                        << '\t' << pStatistics->fHQAltBias
                        << '\t' << float(double(pStatistics->iHQDepth)/double(max(pStatistics->iDepth,1)))
                        << '\t' << pStatistics->iXADepth
                        << '\t' << pStatistics->iXAAltDepth;

            }
            
            cout << '\n';
        }
    }
}
//----------------------------------------------------------------
CAlgorithm::CAlgorithm(void)
{
    bCountDuplicates=false;
    bCountSupplementary=false;
    bPileupRegions=false;
    bVerbose=false;
    
    iMinBaseScore=30;
    iMinAlignmentScore=40;
    iMinMatchLength=10;
    nThreads=1;
}
//----------------------------------------------------------------
bool CAlgorithm::Run(void)
{
    //----------------------------------------------------------------
    //Open Annovar file
    //----------------------------------------------------------------
    
    if(bVerbose) cerr << "Info: Open annovar file" << endl;
    
    if(annovarFile.Open(sAnnovarFileName)==false)
    {
        return false;
    }
    
    //----------------------------------------------------------------
    //Create output matrix
    //----------------------------------------------------------------
    
    if(bVerbose) cerr << "Info: Create output matrix" << endl;
    
    CreateOutputMatrix();
    
    //----------------------------------------------------------------
    //Create jobs
    //----------------------------------------------------------------
    
    if(bVerbose) cerr << "Info: Create jobs" << endl;
    
    CreateJobs();

    //----------------------------------------------------------------
    //Create jobs
    //----------------------------------------------------------------
    
    if(bVerbose) cerr << "Info: Run jobs" << endl;
    
    RunJobs();
    
    //----------------------------------------------------------------
    //Output statistics
    //----------------------------------------------------------------
    
    if(bVerbose) cerr << "Info: Write output matrix" << endl;
    
    WriteOutputMatrix();
    
    //----------------------------------------------------------------
    //Clean up memory
    //----------------------------------------------------------------
    
    if(bVerbose) cerr << "Info: Clean up" << endl;
    
    delete [] (STATISTICS*)annovarFile.mEntries.begin()->second.begin()->second.pUserData;
    
    //----------------------------------------------------------------
    //Done
    //----------------------------------------------------------------
    
    if(bVerbose) cerr << "Info: Done" << endl;
    
    return true;
}
//----------------------------------------------------------------
