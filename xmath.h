#ifndef XMATH_H
#define XMATH_H

#include <ostream>

using namespace std;

class XMath
{
    public:
        static const float smallValue;
        
    public:
        XMath(){}
        virtual ~XMath(){}
        static int    sign(const float val);
        static float  xlog(const float base, const float val);
        static bool   transpose(double *mat, const int col, const int row);
        static void   print_matrix(const char *msg,  ostream *out_strm,
                                 const int num, const int dim, const double* a);

        static void   normVects(double *mat,  const int dim, const int row);
        static void   normVects(float *vects, const unsigned int d0, const unsigned int n0, float *lens);
        static void   test();

};

#endif
