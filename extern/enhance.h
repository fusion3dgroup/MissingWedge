
#include <stdio.h>
#include <cv.h>

struct OPTS{

    double maxItr;
    double gamma;
    double beta;
    double relchg_tol;
    char real_sol;
    char normalize;

};

/* Function Declarations */
void enhance(IplImage *image, int n, double m_iter, double m_gamma, double m_beta, double m_tol, double m_aTV, double m_aL1, CvMat *U);
