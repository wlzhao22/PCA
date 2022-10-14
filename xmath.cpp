#include "xmath.h"

#include <iostream>
#include <cassert>
#include <cstring>
#include <cmath>

const float XMath::smallValue   = 0.00001f;

int XMath::sign(const float val)
{
    if(val >= 0)
        return 1;
    else
        return -1;
}


float XMath::xlog(const float base, const float val)
{
    float r = log(val)/log(base);
    return r;
}

void XMath::normVects(double *mat, const int dim, const int row)
{
    assert(mat);
    unsigned long i, j, rowloc = 0;
    double sum = 0;
    for(i = 0; i < row; i++)
    {
        rowloc = i*dim;
        for(j = 0; j < dim; j++)
        {
            sum += mat[rowloc+j]*mat[rowloc+j];
        }
        sum = sqrt(sum);
        if(sum > 0)
        {
            for(j = 0; j < dim; j++)
            {
                mat[rowloc+j] = mat[rowloc+j]/sum;
            }
        }
        sum = 0;
    }

    return ;
}


void XMath::normVects(float *vects, const unsigned int d0,
                      const unsigned int n0, float *lens)
{
    assert(vects);
    assert(lens);
    unsigned int i = 0, j = 0, loc = 0;
    float len = 0;
    for(i = 0; i < n0; i++)
    {
        len = 0;
        for(j = 0; j < d0; j++)
        {
            len += vects[loc+j]*vects[loc+j];
        }
        len     = sqrt(len);
        lens[i] = len;
        if(len > XMath::smallValue)
        {
            for(j = 0; j < d0; j++)
            {
                vects[loc+j] = vects[loc+j]/len;
            }
        }
        else
        {
            lens[i] = 0;
        }
        loc += d0;
    }
    return ;
}

bool XMath::transpose(double *mat, const int col, const int row)
{
    assert(mat);
    double *tmpmat = new double[col*row];
    memcpy(tmpmat, mat, col*row*sizeof(double));
    long loc, c = 0;
    for(unsigned int i = 0; i < row; i++)
    {
        for(unsigned int j = 0; j < col; j++)
        {
            loc      = j*row+i;
            mat[loc] = tmpmat[c];
            c++;
        }
    }

    delete [] tmpmat;
    return true;
}


void XMath::print_matrix(const char *msg,  ostream *out_strm,
                         int num, int dim, const double* a)
{
    unsigned int i = 0, j = 0;
    cout<<msg<<" "<<dim<<"x"<<num<<endl;
    (*out_strm)<<dim<<" "<<num<<endl;
    for( i = 0; i < num; i++ )
    {
        for( j = 0; j < dim; j++ )
            (*out_strm)<<a[i*dim+j]<<" ";
        (*out_strm)<<endl;
    }
    return ;
}

void XMath::test()
{

    cout<<endl;
}
