#ifndef IODELEGATOR_H
#define IODELEGATOR_H

#include <fstream>
#include <vector>
#include <string>
#include <list>
#include <map>
#include <set>

/********************************************************
*
*
* in charge of data input
*
* @author: wlzhao Zhao
* @date:   Jul. 2013
* @email:  stonescx@gmail.com
*
*
********************************************************/

using namespace std;

class IODelegator
{
private:
    static const int LONGCHAR  = 2000;
    static const int FNLEN     = 1024;
    static const int precision = 3;

public:
    IODelegator();
    static double  *loadTEXTMat(const char *srcFn,  unsigned int &num, unsigned int &dim);
    static float   *loadFVECMat(const char *srcFn,  unsigned int &row, unsigned int &col);
    static void    test();

    ~IODelegator();
};

#endif
