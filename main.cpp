//----------------------------------------------------------------
// Name        : main.cpp
// Author      : Remco Hoogenboezem
// Version     :
// Copyright   :
// Description : annotateBamStatistics
//----------------------------------------------------------------
#include <iostream>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>
#include "Algorithm.h"
//----------------------------------------------------------------
#define ANNOVAR_FILE                    'a'
#define BAM_FILES                       'b'
#define MIN_BASE_SCORE                  's'
#define MIN_ALIGNMENT_SCORE             'S'
#define MIN_MATCH_LENGTH                'm'
#define THREADS                         't'
#define COUNT_DUPLICATES                'd'
#define COUNT_SUPPLEMENTARY             'u'
#define PILEUP_REGIONS                  'r'
#define VERBOSE                         'v'
#define HELP                            'h'
#define SHORT_OPTIONS                   "a:b:s:S:m:t:durvh"
//----------------------------------------------------------------
struct option longOptions[] =
{
    {"annovar-file",required_argument,NULL,ANNOVAR_FILE},
    {"bam-files",required_argument,NULL,BAM_FILES},
    {"min-base-score",required_argument,NULL,MIN_BASE_SCORE},
    {"min-alignment-score",required_argument,NULL,MIN_ALIGNMENT_SCORE},
    {"min-match-length",required_argument,NULL,MIN_MATCH_LENGTH},
    {"threads",required_argument,NULL,THREADS},
    {"count-duplicates",no_argument,NULL,COUNT_DUPLICATES},
    {"count-supplementary",no_argument,NULL,COUNT_SUPPLEMENTARY},
    {"pileup-regions",no_argument,NULL,PILEUP_REGIONS},
    {"verbose",no_argument,NULL,VERBOSE},
    {"help",no_argument,NULL,HELP},
    {0, 0, 0, 0}
};
//----------------------------------------------------------------
vector<string> & Tokenize(const string & sString,const char * pDelimiter)
{
    int                     iStrLen;
    char *                  pString;
    char *                  pStringTemp;
    static vector<string>   vTokens;
    
    vTokens.clear();
    
    if((iStrLen=sString.size())>0)
    {
        pString=pStringTemp=(char*)memcpy(new char[iStrLen+1],sString.c_str(),iStrLen+1);

        while(pStringTemp)
        {
            vTokens.push_back(string(strsep(&pStringTemp,pDelimiter)));
        }
    
        delete [] pString;
    }
    
    return vTokens;
}
//----------------------------------------------------------------
int main(int argc,char ** argv)
{
    bool        bShowHelp;
    
    int         iOption;
    int         iOptionIndex;
    
    CAlgorithm  algorithm;    
    
    //----------------------------------------------------------------
    //Get options
    //----------------------------------------------------------------
    
    bShowHelp=(argc==1);
    
    while((iOption=getopt_long(argc,argv,SHORT_OPTIONS,longOptions,&iOptionIndex))>=0)
    {
        switch(iOption)
        {
            case ANNOVAR_FILE:
                
                algorithm.sAnnovarFileName=string(optarg);
                break;
                
            case BAM_FILES:
                
                algorithm.vBamFileNames=Tokenize(string(optarg),",");
                break;
                
            case MIN_BASE_SCORE:
                
                algorithm.iMinBaseScore=atoi(optarg);
                break;
                
            case MIN_ALIGNMENT_SCORE:
                
                algorithm.iMinAlignmentScore=atoi(optarg);
                break;
                
            case MIN_MATCH_LENGTH:
                
                algorithm.iMinMatchLength=atoi(optarg);
                break;
                
            case THREADS:
                
                algorithm.nThreads=atoi(optarg);
                break;
                
            case COUNT_DUPLICATES:
                
                algorithm.bCountDuplicates=true;
                break;
                
            case COUNT_SUPPLEMENTARY:
                
                algorithm.bCountSupplementary=true;
                break;
                
            case PILEUP_REGIONS:
                
                algorithm.bPileupRegions=true;
                break;
                
            case VERBOSE:
                
                algorithm.bVerbose=true;
                break;
                
            case HELP:
            default:
                
                bShowHelp=true;
                break;
        }
    }
    
    //----------------------------------------------------------------
    //Show help
    //----------------------------------------------------------------

    if(bShowHelp)
    {
        cerr << "AnnotateBAMStatistics [options] > annovarFileOut.txt (output always to std::out)"                                        << endl;
        cerr                                                                                                                                        << endl;
        cerr << "-a --annovar-file <text>           Single annovar file from table_annovar.pl script (required)"                                    << endl;
        cerr << "-b --bam-files <text>              One or more bam files (required)"                                                               << endl;
        cerr << "-s --min-base-score <int>          Minimum base score for filtered statistics (optional,default=30)"                               << endl;
        cerr << "-S --min-alignment-score <int>     Minimum alignment score for filtered statistics (optional,default=40)"                          << endl;
        cerr << "-m --min-match-length <int>        Minimum match length (optional default=10)"                                                     << endl;
        cerr << "-t --threads <int>                 Number of threads to use (optional default=1)"                                                  << endl;
        cerr << "-d --count-duplicates <void>       If specified duplicates are used in the statistics (optional default=false)"                    << endl;
        cerr << "-u --count-supplementary <void>    If specified supplementary alignments are used in the statistics (optional default=false)"      << endl;
        cerr << "-r --pileup-regions <void>         If specified individual regions from the annovar files are piled up (optional default=false)"   << endl;
        cerr << "-v --verbose <void>                If specified be verbose (optional default=false)"                                               << endl;
        cerr                                                                                                                                        << endl;
        
        return 0;
    }
    
    //----------------------------------------------------------------
    //Check options
    //----------------------------------------------------------------
    
    if(algorithm.sAnnovarFileName.empty())
    {
        cerr << "Error: Please specify an annovar file! (-h for help)" << endl;
        return 1;
    }
    
    if(algorithm.vBamFileNames.size()==0)
    {
        cerr << "Error: Please specify at least one bam file! (-h for help)" << endl;
        return 1;
    }
    
    //----------------------------------------------------------------
    //Lets go
    //----------------------------------------------------------------
    
    if(algorithm.Run()==false)
    {
        return 1;
    }
    
    //----------------------------------------------------------------
    //Done
    //----------------------------------------------------------------
    
    return 0;
}
//----------------------------------------------------------------
