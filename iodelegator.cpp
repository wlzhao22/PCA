#include "iodelegator.h"

#include "vstring.h"

#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <cstdio>
#include <cmath>
#include <map>

using namespace std;

double *IODelegator::loadTEXTMat(const char *srcFn, unsigned int &row, unsigned int &col)
{
    assert(srcFn);
    int vals[2] = {0};

    ifstream *inStrm = new ifstream();
    inStrm->open(srcFn, ios::in);
    if(inStrm->fail())
    {
        cout<<"Fail to read "<<srcFn<<endl;
        delete inStrm;
        exit(0);
    }

    (*inStrm)>>vals[0];
    (*inStrm)>>vals[1];
    assert(vals[0] > 0 && vals[1] > 0);

    row = (unsigned int)vals[0];
    col = (unsigned int)vals[1];
    double *mat = new double[row*col];
    unsigned int irow = 0, idim = 0, loc = 0;

    for(irow = 0; irow < row; irow++)
    {
        loc = irow*col;
        for(idim = 0; idim < col; idim++)
        {
            (*inStrm) >>mat[loc + idim];
        }
    }
    inStrm->close();
    delete inStrm;
    return mat;
}

float *IODelegator::loadFVECMat(const char *srcFn, unsigned int &row, unsigned int &col)
{
    float *matrix = NULL, *pMat = NULL;
    unsigned int i = 0;
    ifstream *inStrm = new ifstream(srcFn, ios::in|ios::binary);

    if(!inStrm->is_open())
    {
        row = col = 0;
        cout<<"File "<<srcFn<<" cannot be found!\n";
        return NULL;
    }
    long bg = inStrm->tellg();
    inStrm->read((char*)&col, sizeof(unsigned int));
    inStrm->seekg(0, ios::end);
    long sz = inStrm->tellg() - bg;
    inStrm->close();
    row = sz/(sizeof(unsigned int) + sizeof(float)*col);

    if(row <= 0)
    {
        cout<<"The num. of rows of fvecs has been indicated as '0'!\n";
        return NULL;
    }

    matrix = new float[row*col];
    pMat   = matrix;
    float *tmpVect = new float[col];

    inStrm->open(srcFn, ios::in|ios::binary);
    while(!inStrm->eof() && i < row)
    {
        pMat = matrix + i*col;
        inStrm->read((char*)&col, sizeof(unsigned int));
        inStrm->read((char*)tmpVect, sizeof(float)*col);
        memcpy(pMat, tmpVect, sizeof(float)*col);
        i++;
    }
    inStrm->close();
    delete [] tmpVect;
    tmpVect = NULL;

    return matrix;
}

void IODelegator::test()
{
    const char *itmsfn = "/home/wlzhao/datasets/cvlad/cvlad_imgnet_hesaff_fsv.txt";
    const char *srcfn  = "/home/wlzhao/datasets/vlad/hesaff/vlad_flickr1m";
    const char *fvecs  = "/home/wlzhao/datasets/vocab/cvlad/fvlad_trn16_dense_nrrsift_pca_";

    const char *dstfn  = "/home/wlzhao/datasets/paris/vocab/flickr60k_voc_rsift_inria_64.txt";
    const char *lstfn  = "/home/wlzhao/datasets/ksh/etc/trn_lst.txt";
    const char *dstlst = "/home/wlzhao/datasets/ksh/etc/ksh_trn1k2.txt";

    char src[1024], dst[1024];
    unsigned int r, c, idx = 0;

    return ;
}
