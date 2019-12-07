//----------------------------------------------------------------
// Name        : Mutex.cpp
// Author      : Remco Hoogenboezem
// Version     :
// Copyright   :
// Description :
//----------------------------------------------------------------
#include "Mutex.h"
//----------------------------------------------------------------
CMutex::CMutex(void)
{
    pthread_mutex_init(&mutex,NULL);
}
//----------------------------------------------------------------
CMutex::~CMutex(void)
{
    pthread_mutex_destroy(&mutex);
}
//----------------------------------------------------------------
void CMutex::Lock(void)
{
    pthread_mutex_lock(&mutex);
}
//----------------------------------------------------------------
void CMutex::Unlock(void)
{
    pthread_mutex_unlock(&mutex);
}
//----------------------------------------------------------------