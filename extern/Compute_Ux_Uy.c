//#include "mex.h"
// #include <math.h>
#include <cv.h>

/*
 INPUT: U
 OUTPUT: Ux, Uy

 Equivalent to the following code:
         Ux = [diff(U,1,2), U(:,1) - U(:,n)];
         Uy = [diff(U,1,1); U(1,:)-U(m,:)];

*/
void Compute_Ux_Uy(CvMat* U, CvMat* Ux, CvMat* Uy)
{

    /* matlab : Ux = [diff(U,1,2), U(:,1) - U(:,n)];
                     diff_blank                     */
    int c,r;

    Ux = cvCreateMat(U->rows,U->cols, CV_32F);
    //CvMat* UX_blank = cvCreateMat(Ux->rows,Ux->cols, CV_32F);
    //int inp = U->rows;
    for(c = 0; c <U->cols; c++){
        for(r = 0; r < U->rows-1; r++){
            float elements_01 = cvmGet(U, r,c);
            float elements_02 = cvmGet(U, r+1, c);
            float result = elements_02 - elements_01;
            cvmSet(Ux,r,c,result);

            float elements_03 = cvmGet(U,0,c);
            float elements_06 = cvmGet(U,U->rows-1,c);
            float result_02 = elements_03 - elements_06;
            cvmSet(Ux,U->rows-1,c,result_02);
            }
    }


    /*  matlab : Uy = [diff(U,1,1); U(1,:)-U(m,:)];*/
     Uy = cvCreateMat(U->rows,U->cols, CV_32F);
    for(c = 0; c <U->cols-1; c++){
        for(r = 0; r < U->rows; r++){
            float elements_01 = cvmGet(U, r,c);
            float elements_02 = cvmGet(U, r, c+1);
            float result = elements_02 - elements_01;
            cvmSet(Uy,r,c, result);

            float elements_03 = cvmGet(U,r,0);
            float elements_04 = cvmGet(U,r,U->cols-1);
            float result_02 = elements_03 - elements_04;
            cvmSet(Uy,r,U->cols-1, result_02);
        }
     }

}

