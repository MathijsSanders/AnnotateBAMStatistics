//----------------------------------------------------------------
#ifndef AnnovarFileH
#define	AnnovarFileH
//----------------------------------------------------------------
#include <string>
#include <map>
//----------------------------------------------------------------
using namespace std;
//----------------------------------------------------------------
typedef struct
{
    int     iAlt;
    int     iLen;
    void *  pUserData;
    string  sLine;
    
}ENTRY;
//----------------------------------------------------------------
class CAnnovarFile
{
    private:
        
    public:

        int                                 nEntries;
        string                              sHeader;
        map<string,multimap<int,ENTRY> >    mEntries;
        
        bool Open(const string & sFileName);
};
//----------------------------------------------------------------
#endif
