#include "pcatrainer.h"
#include "iodelegator.h"
#include "vstring.h"
#include "timer.h"

#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <vector>
#include <cstdio>
#include <cmath>

using namespace std;
#include "xmath.h"

#ifdef __cplusplus
extern "C" {
#endif
///call svd from lapack
void dgesvd_(char* jobu, char* jobvt,  int* m,     int* n,   double* a,
             int* lda,   double* s,    double* u,  int* ldu, double* vt,
             int* ldvt,  double* work, int* lwork, int* info);

void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda,
             double* w, double* work, int* lwork, int* info );
#define AAA

#ifdef __cplusplus
}
#endif

bool PCATrainer::pca_learn(const char *srcfn, const char *dstfn)
{
    cout<<srcfn<<endl;
    unsigned int num  = 0, dim = 0, i, j, loc;
    cout<<"Loading matrix .......................... ";
    double *mat  = IODelegator::loadTEXTMat(srcfn, num, dim);
    cout<<num<<"x"<<dim<<endl;

    cout<<"Allocate memory ......................... ";
    double *u    = new double[num*num];
    double *vt   = new double[dim*dim];

    double *s    = new double[dim];
    double *avg  = new double[dim];
    double *work = NULL;
    double  wkopt = 0;
    cout<<"succeed\n";

    memset(avg, 0, sizeof(double)*dim);
    ofstream dst_strm;
    dst_strm.open(dstfn, ios::out);
    if(!dst_strm.is_open())
    {
        cerr<<"Destine file cannot for write!\n";
        exit(1);
    }
    else
    {
        dst_strm.close();
    }

    for(i = 0; i < num; i++)
    {
        loc   = dim*i;
        for(j = 0; j < dim; j++)
        {
            avg[j] += mat[loc+j]/num;
        }
    }

    for(i = 0; i < num; i++)
    {
        loc   = dim*i;
        for(j = 0; j < dim; j++)
        {
            mat[loc+j] = mat[loc+j] - avg[j];
        }
    }

    double rnum = sqrt(num - 1.0f);
    for(i = 0; i < num; i++)
    {
        loc = dim*i;
        for(j = 0; j < dim; j++)
        {
            mat[loc+j] = mat[loc+j]/rnum;
        }
    }

    ///m: col; n: row
    int m   = num, lda = num,  n = dim;
    int ldu = num, ldvt = dim, info = 0, lwork = -1;

    /**perform SVD on mat^t**/
    cout<<"Perform PCA training .................... ";
    XMath::transpose(mat, dim, num);
    dgesvd_("All", "All", &m, &n, mat, &lda, s, u,
            &ldu, vt, &ldvt, &wkopt, &lwork, &info);

    lwork = (int)wkopt;
    work  = new double[lwork];
    dgesvd_("All", "All", &m, &n, mat, &lda, s, u,
            &ldu, vt, &ldvt, work, &lwork, &info);
    /********end of SVD******/

    cout<<"done\n";

    if(info > 0)
    {
        cerr<<"The algorithm computing SVD failed to converge.\n";
        exit(1);
    }

    dst_strm.open(dstfn, ios::out);
    XMath::transpose(vt, n, n);// vt is already in row-wise
    if(!dst_strm.is_open())
    {
        cerr<<"Destine file cannot for write!\n";
        exit(1);
    }
    else
    {
        ///print_matrix("Eigen values: ", &cout, 1, n, s, 1);
        print_matrix("", &dst_strm, n, n, vt, s, avg, n);
        dst_strm.close();
    }

    delete [] s;
    delete [] u;
    delete [] vt;
    delete [] avg;
    delete [] mat;
    delete [] work;
    return true;
}

bool PCATrainer::pcaLearn(const char *srcfn, const char *dstfn)
{
    cout<<"Loading matrix .......................... ";
    unsigned int num = 0, dim = 0;
    double *mat  = IODelegator::loadTEXTMat(srcfn, num, dim);
    cout<<num<<"x"<<dim<<endl;
    assert(num > 0);
    assert(dim > 0);
    Timer *mytm = new Timer();
    mytm->start();
    if(num >= dim)
    {
        PCATrainer::pcaByCovSVD(mat, num, dim, dstfn);
    }
    else
    {
        PCATrainer::pcaByGramEig(mat, num, dim, dstfn);
    }
    mytm->end(true);
    delete mytm;
    delete [] mat;

    return true;
}

bool PCATrainer::pcaByCovSVD(double *mat, const unsigned int nRow, const unsigned int nCol, const char *dstfn)
{
    unsigned int num  = nRow, dim = nCol, i, j, loc;
    cout<<"Allocate memory ......................... ";
    double *u    = new double[num*num];
    double *vt   = new double[dim*dim];

    double *eigV = new double[dim];
    double *avg  = new double[dim];
    double *work = NULL;
    double  wkopt = 0;
    cout<<"succeed\n";

    memset(avg, 0, sizeof(double)*dim);
    ofstream dst_strm;
    dst_strm.open(dstfn, ios::out);
    if(!dst_strm.is_open())
    {
        cerr<<"Destine file cannot for write!\n";
        exit(1);
    }
    else
    {
        dst_strm.close();
    }

    for(i = 0; i < num; i++)
    {
        loc   = dim*i;
        for(j = 0; j < dim; j++)
        {
            avg[j] += mat[loc+j]/num;
        }
    }

    for(i = 0; i < num; i++)
    {
        loc   = dim*i;
        for(j = 0; j < dim; j++)
        {
            mat[loc+j] = mat[loc+j] - avg[j];
        }
    }

    double rnum = sqrt(num - 1.0f);
    for(i = 0; i < num; i++)
    {
        loc = dim*i;
        for(j = 0; j < dim; j++)
        {
            mat[loc+j] = mat[loc+j]/rnum;
        }
    }
    XMath::transpose(mat, dim, num);

    ///m: col; n: row
    int m   = num, lda = num,  n = dim;
    int ldu = num, ldvt = dim, info = 0, lwork = -1;

    /**perform SVD on mat^t**/
    cout<<"Perform PCA training (by COV-SVD) ....... ";

    dgesvd_("All", "All", &m, &n, mat, &lda, eigV, u,
            &ldu, vt, &ldvt, &wkopt, &lwork, &info);


    lwork = (int)wkopt;
    work  = new double[lwork];
    dgesvd_("All", "All", &m, &n, mat, &lda, eigV, u,
            &ldu, vt, &ldvt, work, &lwork, &info);
    /********end of SVD******/

    cout<<"done\n";

    if(info > 0)
    {
        cerr<<"The algorithm computing SVD failed to converge.\n";
        exit(1);
    }


    XMath::transpose(vt, n, n); /// vt is already in row-wise
    if(VString::endWith(dstfn, ".txt"))
    {
        dst_strm.open(dstfn, ios::out);
    }
    else if(VString::endWith(dstfn, ".fvecs"))
    {
        dst_strm.open(dstfn, ios::out|ios::binary);
    }
    if(!dst_strm.is_open())
    {
        cerr<<"Destine file cannot for write!\n";
        exit(1);
    }
    else
    {
        /*** redundant*
        for(i = 0; i < dim; i++)
        {
            eigV[i] = sqrt(eigV[i]);
        }
        /**/
        if(VString::endWith(dstfn, ".txt"))
        {
            print_matrix("", &dst_strm, dim, dim, vt, eigV, avg, dim);
        }
        else if(VString::endWith(dstfn, ".fvecs"))
        {
            printFvecs(&dst_strm, dim, dim, vt, eigV, avg, dim);
        }
        dst_strm.close();
    }

    delete [] eigV;
    delete [] u;
    delete [] vt;
    delete [] avg;
    delete [] work;
    return true;
}

bool PCATrainer::pcaByCovEig(double *mat, const unsigned int nRow, const unsigned int nCol, const char *dstfn)
{
    unsigned int num  = nRow, dim = nCol, i, j, k, loc;
    cout<<"Allocate memory ......................... ";
    double *cov  = new double[dim*dim];
    double *eigV = new double[dim];
    double *avg  = new double[dim];
    double *work = NULL;
    double  wkopt = 0;
    cout<<"succeed\n";

    memset(avg, 0, sizeof(double)*dim);
    ofstream dst_strm;
    dst_strm.open(dstfn, ios::out);
    if(!dst_strm.is_open())
    {
        cerr<<"Destine file cannot for write!\n";
        exit(1);
    }
    else
    {
        dst_strm.close();
    }

    for(i = 0; i < num; i++)
    {
        loc   = dim*i;
        for(j = 0; j < dim; j++)
        {
            avg[j] += mat[loc+j]/num;
        }
    }

    for(i = 0; i < num; i++)
    {
        loc   = dim*i;
        for(j = 0; j < dim; j++)
        {
            mat[loc+j] = mat[loc+j] - avg[j];
        }
    }

    double rnum = num + 0.0f;
    memset(cov, 0, sizeof(double)*dim*dim);
    for(i = 0; i < dim; i++)
    {
        loc = dim*i;
        for(j = 0; j < dim; j++)
        {
            for(k = 0; k < num; k++)
            {
                cov[loc+j] += mat[k*dim+i]*mat[k*dim+j];
            }
        }
    }

    for(i = 0; i < dim; i++)
    {
        loc = dim*i;
        for(j = 0; j < dim; j++)
        {
            cov[loc+j] = cov[loc+j]/rnum;
        }
    }

    /**perform eigenvalue decomposition on cov**/
    int lda = dim, m = dim, info = 0, lwork = -1;

    cout<<"Perform PCA training (By COV-Eig) .......... ";
    dsyev_("Vectors", "Upper", &m, cov, &lda, eigV, &wkopt, &lwork, &info);

    lwork = (int)round(wkopt);
    work = new double[lwork*sizeof(double)];
    dsyev_("Vectors", "Upper", &m, cov, &lda, eigV, work, &lwork, &info);
    cout<<"done\n";
    if(info > 0)
    {
        cerr<<"The algorithm computing eigenvectors failed to converge.\n";
        exit(1);
    }
    /********end of eigenvalue decomposition******/

    if(info > 0)
    {
        cerr<<"The algorithm computing eigenvalue failed to converge.\n";
        exit(1);
    }

    vector<IndexItem*>::iterator vit;
    vector<IndexItem*> EigsOrder;
    IndexItem *crntItm = NULL;
    for(i = 0; i < dim; i++)
    {
        float temp = (float)eigV[i];
        if(temp > 0.0000001f)
        {
            crntItm = new IndexItem(i, temp);
            EigsOrder.push_back(crntItm);
        }
    }
    stable_sort(EigsOrder.begin(), EigsOrder.end(), IndexItem::LGcomparer);

    m = EigsOrder.size();
    ///cout<<"Size: "<<m<<"x"<<dim<<endl;
    double *pcaMat = new double[m*dim];
    memset(pcaMat, 0, sizeof(double)*m*dim);
    for(vit = EigsOrder.begin(), i = 0; vit != EigsOrder.end(); vit++, i++)
    {
        crntItm = *vit;
        loc = crntItm->Id*dim;
        memcpy(pcaMat+i*dim, cov+loc, sizeof(double)*dim);
        eigV[i] = crntItm->val;
    }

    for(vit = EigsOrder.begin(); vit != EigsOrder.end(); vit++)
    {
        crntItm = *vit;
        delete crntItm;
    }
    EigsOrder.clear();

    if(VString::endWith(dstfn, ".txt"))
    {
        dst_strm.open(dstfn, ios::out);
    }
    else if(VString::endWith(dstfn, ".fvecs"))
    {
        dst_strm.open(dstfn, ios::out|ios::binary);
    }
    if(!dst_strm.is_open())
    {
        cerr<<"Destine file cannot for write!\n";
        exit(1);
    }
    else
    {
        if(VString::endWith(dstfn, ".txt"))
        {
            print_matrix("", &dst_strm, m, dim, pcaMat, eigV, avg, dim);
        }
        else if(VString::endWith(dstfn, ".fvecs"))
        {
            printFvecs(&dst_strm, m, dim, pcaMat, eigV, avg, dim);
        }
        dst_strm.close();
    }

    delete [] eigV;
    delete [] cov;
    delete [] avg;
    delete [] work;
    delete [] pcaMat;
    return true;
}

bool PCATrainer::pcaByGramSVD(double *mat, const unsigned int nRow, const unsigned int nCol, const char *dstfn)
{
    unsigned int num  = nRow, dim = nCol, i, j, k, loc;

    cout<<"Allocate memory ......................... ";
    double *u    = new double[num*num];
    double *vt   = new double[dim*dim];
    double *gram = new double[dim*num];
    double *avg  = new double[dim];
    double *eigV = new double[num];
    double *work = NULL;
    cout<<"succeed\n";

    memset(avg,  0, sizeof(double)*dim);
    memset(u,    0, sizeof(double)*num*num);
    memset(vt,   0, sizeof(double)*dim*dim);
    memset(eigV, 0, sizeof(double)*num);

    ofstream dst_strm;
    dst_strm.open(dstfn, ios::out);
    if(!dst_strm.is_open())
    {
        cerr<<"Destine file cannot open for write!\n";
        exit(1);
    }
    else
    {
        dst_strm.close();
    }

    for(i = 0; i < num; i++)
    {
        loc   = dim*i;
        for(j = 0; j < dim; j++)
        {
            avg[j] += mat[loc+j]/num;
        }
    }

    for(i = 0; i < num; i++)
    {
        loc   = dim*i;
        for(j = 0; j < dim; j++)
        {
            mat[loc+j] = mat[loc+j] - avg[j];
        }
    }

    ///m: col; n: row
    int m   = num, lda = num,  n = dim;
    int ldu = num, ldvt = dim, info = 0, lwork = -1;
    double wkopt;

    cout<<"Perform PCA training (By Gram) .......... " << endl;
    dgesvd_("All", "All", &m, &n, mat, &lda, eigV, u,
            &ldu, vt, &ldvt, &wkopt, &lwork, &info);

    if(info > 0)
    {
        cerr<<"The algorithm computing eigenvectors failed to converge at its first step.\n";
        exit(1);
    }

    lwork = (int)wkopt;
    work  = new double[lwork];
    dgesvd_("All", "All", &m, &n, mat, &lda, eigV, u,
            &ldu, vt, &ldvt, work, &lwork, &info);
    /********end of SVD******/

    if(info > 0)
    {
        cerr<<"The algorithm computing eigenvectors failed to converge at its second step.\n";
        exit(1);
    }

    if(VString::endWith(dstfn, ".txt"))
    {
        dst_strm.open(dstfn, ios::out);
    }
    else if(VString::endWith(dstfn, ".fvecs"))
    {
        dst_strm.open(dstfn, ios::out|ios::binary);
    }
    if(!dst_strm.is_open())
    {
        cerr<<"Destine file cannot for write!\n";
        exit(1);
    }
    else
    {
        memset(gram, 0, sizeof(double)*num*dim);
        ///XMath::transpose(u,   num, num);
        XMath::transpose(mat, dim, num);
        for(i = 0; i < dim; i++)
        {
            for(j = 0; j < num; j++)
            {
                for(k = 0; k < num; k++)
                {
                    gram[i*num+j] += mat[i*num+k]*u[j*num+k];
                }
            }
        }

        for(i = 0; i < dim; i++)
        {
            eigV[i] = sqrt(eigV[i]);
        }

        XMath::transpose(gram, dim, num);
        XMath::normVects(gram, dim, num);
        if(VString::endWith(dstfn, ".txt"))
        {
            print_matrix("", &dst_strm, num, dim, gram, eigV, avg, m);
        }
        else if(VString::endWith(dstfn, ".fvecs"))
        {
            printFvecs(&dst_strm, num, dim, gram, eigV, avg, m);
        }
        dst_strm.close();
    }

    /**
    ofstream *outStrm = new ofstream("/home/wlzhao/datasets/hi.txt", ios::out);
    for(j = 0; j < num; j++)
    {
      for(k = 0; k < num; k++)
      {
          (*outStrm)<<u[num*j+k]<<" ";
      }
      (*outStrm)<<endl;
    }
    outStrm->close();
    /**/

    delete [] u;
    delete [] avg;
    delete [] vt;
    delete [] work;
    delete [] gram;
    delete [] eigV;
    return true;
}


bool PCATrainer::pcaByGramEig(double *mat, const unsigned int nRow, const unsigned int nCol, const char *dstfn)
{
    unsigned int num  = nRow, dim = nCol, i, j, k, loc, ri, rj;

    cout<<"Allocate memory ......................... ";
    double *gram = new double[num*num];
    double *vt = new double[num*dim];
    double *avg  = new double[dim];
    double *work = NULL;
    double *pcaMat = NULL;
    cout<<"succeed\n";

    memset(avg,  0, sizeof(double)*dim);
    memset(gram, 0, sizeof(double)*num*num);
    memset(vt,   0, sizeof(double)*num*dim);

    ofstream dst_strm;
    dst_strm.open(dstfn, ios::out);
    if(!dst_strm.is_open())
    {
        cerr<<"Destine file cannot for write!\n";
        exit(1);
    }
    else
    {
        dst_strm.close();
    }

    for(i = 0; i < num; i++)
    {
        loc   = dim*i;
        for(j = 0; j < dim; j++)
        {
            avg[j] += mat[loc+j]/num;
        }
    }

    for(i = 0; i < num; i++)
    {
        loc   = dim*i;
        for(j = 0; j < dim; j++)
        {
            mat[loc+j] = mat[loc+j] - avg[j];
        }
    }

    /**Compute Gram matrix**/
    for(i = 0; i < num; i++)
    {
        loc = i*num;
        ri = i*dim;
        for(j = 0; j < num; j++)
        {
            rj = j*dim;
            for(k = 0; k < dim; k++)
            {
                gram[loc+j] += mat[ri+k]*mat[rj+k];
            }
        }
    }

    int lda = num, m = num, info = 0, lwork = -1;
    double wkopt = 0;
    double *eigV = new double[m];

    cout<<"Perform PCA training (By Gram) .......... ";
    dsyev_("Vectors", "Upper", &m, gram, &lda, eigV, &wkopt, &lwork, &info);

    lwork = (int)round(wkopt);
    work = new double[lwork*sizeof(double)];
    dsyev_("Vectors", "Upper", &m, gram, &lda, eigV, work, &lwork, &info);
    cout<<"done\n";
    if(info > 0)
    {
        cerr<<"The algorithm computing eigenvectors failed to converge.\n";
        exit(1);
    }

    vector<IndexItem*>::iterator vit;
    vector<IndexItem*> EigsOrder;
    IndexItem *crntItm = NULL;
    for(i = 0; i < num; i++)
    {
        float temp = (float)eigV[i];
        if(temp > 0.000000f)
        {
            crntItm = new IndexItem(i, temp);
            EigsOrder.push_back(crntItm);
        }
    }
    stable_sort(EigsOrder.begin(), EigsOrder.end(), IndexItem::LGcomparer);

    if(VString::endWith(dstfn, ".txt"))
    {
        dst_strm.open(dstfn, ios::out);
    }
    else if(VString::endWith(dstfn, ".fvecs"))
    {
        dst_strm.open(dstfn, ios::out|ios::binary);
    }
    if(!dst_strm.is_open())
    {
        cerr<<"Destine file cannot for write!\n";
        exit(1);
    }
    else
    {
        cout<<"=============== Post-processing ========= "<<endl;
        cout<<"Multiply with original matrix ........... ";
        for(i = 0; i < dim; i++)
        {
            loc = i*num;
            for(j = 0; j < num; j++)
            {
                for(k = 0; k < num; k++)
                {
                    vt[loc+j] += mat[k*dim+i]*gram[j*num+k];
                }
            }
        }
        cout<<"done\n";
        cout<<"Transpose mapping matrix ................ ";
        XMath::transpose(vt, num, dim);
        XMath::normVects(vt, dim, num);
        cout<<"done\n";
        cout<<"Sort mapping matrix ..................... ";
        m = EigsOrder.size();
        pcaMat = new double[m*dim];
        memset(pcaMat, 0, sizeof(double)*m*dim);
        for(vit = EigsOrder.begin(), i = 0; vit != EigsOrder.end(); vit++, i++)
        {
            crntItm = *vit;
            loc = crntItm->Id*dim;
            memcpy(pcaMat+i*dim, vt+loc, sizeof(double)*dim);
            eigV[i] = crntItm->val;
        }
        cout<<"done\n";

        if(VString::endWith(dstfn, ".txt"))
        {
            print_matrix("", &dst_strm, m, dim, pcaMat, eigV, avg, m);
        }
        else if(VString::endWith(dstfn, ".fvecs"))
        {
            printFvecs(&dst_strm, m, dim, pcaMat, eigV, avg, m);
        }
        dst_strm.close();
    }

    delete [] gram;
    delete [] avg;
    delete [] vt;
    delete [] work;
    delete [] eigV;
    if(pcaMat != NULL)
        delete [] pcaMat;
    for(vit = EigsOrder.begin(); vit != EigsOrder.end(); vit++)
    {
        crntItm = *vit;
        delete crntItm;
    }
    EigsOrder.clear();
    return true;
}

void PCATrainer::print_matrix(const char *msg,  ostream *out_strm, int row,
                              int col, const double* a, const double *eigV,
                              const double *avg, int lda)
{
    int i, j , loc = 0;

    cout<<"Matrix size ............................. "<<row<<"x"<<col<<endl;

    (*out_strm)<<row<<" "<<col<<endl;
    for( j = 0; j < col; j++ )
    {
        (*out_strm)<<avg[j]<<" ";
    }
    (*out_strm)<<endl;

    for( j = 0; j < row; j++ )
    {
        ///assert(eigV[j] > 0);
        (*out_strm)<<(float)eigV[j]<<" ";
    }

    (*out_strm)<<endl;

    for( i = 0; i < row; i++ )
    {
        loc = i*col;
        for( j = 0; j < col; j++ )
        {
            (*out_strm)<<(float)a[loc+j]<<" ";
        }
        (*out_strm)<<endl;
    }
    return ;
}

void PCATrainer::print_matrix(const char *msg,  ostream *out_strm, int m,
                              int n, const double* mat)
{
    int i, j , loc = 0;
    cout<<msg<<" "<<m<<"x"<<n<<endl;
    (*out_strm)<<m<<" "<<n<<endl;

    for( i = 0; i < m; i++ )
    {
        loc = i*n;
        for( j = 0; j < n; j++ )
        {
            (*out_strm)<<mat[loc+j]<<" ";
        }
        (*out_strm)<<endl;
    }
    return ;
}

void PCATrainer::printFvecs(ostream *outStrm, int row,
                            int col, const double* a, const double *eigV,
                            const double *avg, int lda)
{
    int i, j , loc = 0;

    cout<<"Matrix size ............................. "<<row<<"x"<<col<<endl;
    float *mean = new float[col];
    for(j = 0; j < col; j++ )
    {
        mean[j] = (float)avg[j];
    }
    outStrm->write((char*)&col, sizeof(int));
    outStrm->write((char*)mean, sizeof(float)*col);

    for(j = 0; j < col; j++ )
    {
        mean[j] = (float)eigV[j];
    }
    outStrm->write((char*)&col, sizeof(int));
    outStrm->write((char*)mean, sizeof(float)*col);

    for(i = 0; i < row; i++ )
    {
        loc = i*col;
        for(j = 0; j < col; j++ )
        {
            mean[j] = (float)a[loc+j];
        }
        outStrm->write((char*)&col, sizeof(int));
        outStrm->write((char*)mean, sizeof(float)*col);
    }
    delete [] mean;
    return ;

}

void PCATrainer::batchPCALearn(const char *srcdir, const char *dstdir, const char *fn)
{
    char srcFn[1024], dstFn[1024];
    for(int i = 0; i < 8; i++)
    {
        sprintf(srcFn, "%s%s_mat_%d.txt",   srcdir, fn, i);
        sprintf(dstFn, "%s%s_pca_%d.fvecs", dstdir, fn, i);
        cout<<srcFn<<endl;
        PCATrainer::pcaLearn(srcFn, dstFn);
    }
}

void PCATrainer::test()
{
    const char *srcfn    = "/home/wlzhao/datasets/bignn/sift1m/sift4k.txt";
    const char *dstfn    = "/home/wlzhao/datasets/bignn/sift1m/pcasift.txt";

    PCATrainer::pcaLearn(srcfn, dstfn);

}
