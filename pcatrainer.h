#ifndef PCATRAINER_H
#define PCATRAINER_H

/**
PCA analysis

@author:    Wan-Lei Zhao
@institute: Xiamen University
@date:      Mar-23-2015
@updated:   Oct-13-2022

given input matrix in nRow and nCol,  if nRow >= nCol
the PCA mapping matrix is trained with covariance matrix
otherwise when nRow < nCol, the PCA mapping matrix is learned
by Gram method.

**/

#include <ostream>

using namespace std;

struct IndexItem {

public:
   IndexItem(int id, float val0){
      Id  = id;
      val = val0;
   }
   int Id;
   float val;

    static bool LGcomparer(const IndexItem *a,const IndexItem *b)
    {
        return (a->val > b->val);
    }

    static bool LLcomparer(const IndexItem *a,const IndexItem *b)
    {
        return (a->val < b->val);
    }

};

class PCATrainer
{
    public:
        PCATrainer(){}
        virtual ~PCATrainer(){}
        static bool pca_learn(const char *srcfn, const char *dstfn);
        static bool pcaLearn(const char *srcfn, const char *dstfn);
        static void batchPCALearn(const char *srcdir, const char *dstdir, const char *fn);
        static bool pcaByCovSVD(double *mat,  const unsigned int nRow, const unsigned int nCol, const char *dstfn);
        static bool pcaByCovEig(double *mat,  const unsigned int nRow, const unsigned int nCol, const char *dstfn);
        static bool pcaByGramSVD(double *mat, const unsigned int nRow, const unsigned int nCol, const char *dstfn);
        static bool pcaByGramEig(double *mat, const unsigned int nRow, const unsigned int nCol, const char *dstfn);

        static void print_matrix(const char *msg,  ostream *out_strm, int m,
                                 int n, const double* a, const double *eigV,
                                 const double *avg, int lda);
        static void print_matrix(const char *msg,  ostream *out_strm, int m,
                                 int n, const double* mat);

        static void printFvecs(ostream *outStrm, int row,
                            int col, const double* a, const double *eigV,
                            const double *avg, int lda);

        static void test_svd();
        static void test();
};

#endif
