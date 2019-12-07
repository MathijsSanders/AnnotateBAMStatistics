//----------------------------------------------------------------
// Name        : AnnovarFile.cpp
// Author      : Remco Hoogenboezem
// Version     :
// Copyright   :
// Description :
//----------------------------------------------------------------
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "Alphabet.h"
#include "AnnovarFile.h"
//----------------------------------------------------------------
bool CAnnovarFile::Open(const string & sFileName)
{
    char    cLine[1048576];
   
    FILE *  pAnnovarFile;
    
    char *  pLine;
    char *  pChr;
    char *  pStart;
    char *  pRef;
    char *  pAlt;
    
    ENTRY   entry;
    
    //----------------------------------------------------------------
    //Open annovar file
    //----------------------------------------------------------------
    
    if((pAnnovarFile=popen((string("bzcat -f ")+sFileName+string(" | zcat -f")).c_str(),"r"))==NULL)    //The annovar files can now also be compressed
    {
        cerr << "Error: Could not open annovar file: " << sFileName << endl;
        return false;
    }
    
    //----------------------------------------------------------------
    //Read and check annovar header
    //----------------------------------------------------------------
    
    if(fgets(cLine,1048576,pAnnovarFile)==NULL)
    {
        cerr << "Error: Could not read header from annovar file: " << sFileName << endl;
        pclose(pAnnovarFile);
        return false;
    }
    
    if(strstr(cLine,"Chr\tStart\tEnd\tRef\tAlt")==NULL)
    {
        cerr << "Error: Invalid annovar file header: " << sFileName << endl;
        pclose(pAnnovarFile);
        return false;
    }

    cLine[strcspn(cLine,"\n")]='\0';    //Remove new line
    sHeader=string(cLine);
    
    //----------------------------------------------------------------
    //Read annovar entries
    //----------------------------------------------------------------
    
    nEntries=0;
    mEntries.clear();
            
    entry.pUserData=NULL;
    
    while(fgets(cLine,1048576,pAnnovarFile)!=NULL)
    {
        cLine[strcspn(cLine,"\n")]='\0';
        entry.sLine=string(cLine);    
        
        pLine=cLine;
        pChr=strsep(&pLine,"\t");
        pStart=strsep(&pLine,"\t");
        strsep(&pLine,"\t");
        pRef=strsep(&pLine,"\t");
        pAlt=strsep(&pLine,"\t");
        
        if(pAlt==NULL)
        {
            cerr << "Error: Insufficient number of columns in  annovar file: " << sFileName << endl;
            pclose(pAnnovarFile);
            return false;
        }
        
        //----------------------------------------------------------------
        //Insertion
        //----------------------------------------------------------------
        
        if(pRef[0]=='-')
        {
            entry.iAlt=INS;
            entry.iLen=strlen(pAlt);
            mEntries[string(pChr)].insert(pair<int,ENTRY>(atoi(pStart)-1,entry));
        }

        //----------------------------------------------------------------
        //Deletion
        //----------------------------------------------------------------
        
        else if(pAlt[0]=='-')
        {
            entry.iAlt=DEL;
            entry.iLen=-strlen(pRef);
            mEntries[string(pChr)].insert(pair<int,ENTRY>(atoi(pStart)-2,entry));
        }

        //----------------------------------------------------------------
        //Substitution
        //----------------------------------------------------------------
        
        else
        {
            entry.iAlt=iAscii2BamSeq[pAlt[0]];
            entry.iLen=0;
            mEntries[string(pChr)].insert(pair<int,ENTRY>(atoi(pStart)-1,entry));            
        }

        nEntries++;
    }
    
    pclose(pAnnovarFile);    
    
    //----------------------------------------------------------------
    //Done
    //----------------------------------------------------------------

    return true;
}
//----------------------------------------------------------------