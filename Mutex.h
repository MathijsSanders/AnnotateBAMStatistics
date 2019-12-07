//----------------------------------------------------------------
#ifndef MutexH
#define	MutexH
//----------------------------------------------------------------
#include <pthread.h>
//----------------------------------------------------------------
class CMutex
{
    private:
        
        pthread_mutex_t mutex;
        
    public:
        
        CMutex(void);
        ~CMutex(void);
        
        void Lock(void);
        void Unlock(void);
};
//----------------------------------------------------------------
#endif