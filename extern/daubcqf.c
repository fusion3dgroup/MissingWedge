#include <cv.h>
#include "daubcqf.h"

void daubcqf(int N, CvMat * hH){

    int K = N/2;
    int a=1, p=1, q=1, j=0;

    int h_0[2], p_0[2];
    h_0[0] = 1;
    h_0[1] = 1;

    for(j = 0; j < K-1; j++){
        a = -a*0.25*(j+K-1)/j;
        //h_0 =
        p_0[0] = 1;
        p_0[1] = -1;


    }

    //q = sort(roots(q));

    CvMat *qQ = cvCreateMat(1,6,CV_32FC1);
    //cvSolvePoly()
    double qt = cvmGet(qQ,1,1);

    //cvFilter2D(qQ, hH,)

    hH = cvCreateMat(1,6,CV_32FC1);
    cvmSet(hH,1,1,1) ;
    cvmSet(hH,1,2,1);


}
