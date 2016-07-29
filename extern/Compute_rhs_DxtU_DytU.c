//#include "mex.h"

#include <cv.h>
/*
 INPUT: Wx, Wy, bx, by, tau
 OUTPUT: RHS

 Equivalent to the following code:

    RHS = tau*(DxtU(Wx-bx)+DytU(Wy-by));

    % compute D'_x(U)
    function dxtu = DxtU(U)
        dxtu = [U(:,end)-U(:, 1) U(:,1:end-1)-U(:,2:end)];
    end

    % compute D'_y(U)
    function dytu = DytU(U)
        dytu = [U(end,:)-U(1, :); U(1:end-1,:)-U(2:end,:)];
    end

*/

void Compute_rhs_DxtU_DytU(CvMat*Wx, CvMat* Wy, CvMat* bx, CvMat* by, double tau, CvMat* rhs)
{

    /* matlab : RHS = tau*(DxtU(Wx-bx)+DytU(Wy-by));
                              DxtU         DyTU        */

    CvMat* DxtU_INDEX = cvCreateMat(Wx->rows, Wx->cols, CV_32F);
    CvMat* DytU_INDEX = cvCreateMat(Wy->rows, Wy->cols, CV_32F);

    cvmSub(Wx, bx, DxtU_INDEX);
    cvmSub(Wy, by, DytU_INDEX);

     /* matlab  : dxtu = [U(:,end)-U(:, 1) U(:,1:end-1)-U(:,2:end)];
                            DxtU_INDEX_01         DxtU_INDEX_02     */
    CvMat *DxtU_INDEX_01 = cvCreateMat(1,DxtU_INDEX->cols,CV_32F);
    int c, r;
    for(c = 0; c<DxtU_INDEX->cols; c++){
        /* U(:,end)*/
        float elements_01 = cvmGet(DxtU_INDEX,DxtU_INDEX->rows-1,c);
        /* U(:, 1) */
        float elements_02 = cvmGet(DxtU_INDEX,0,c);
        float result = elements_01 - elements_02;
        cvmSet(DxtU_INDEX_01, 0,c,result);
    }
    CvMat* DxtU_INDEX_02 = cvCreateMat(DxtU_INDEX->rows -1,DxtU_INDEX->cols, CV_32F);
    for(c = 0; c<DxtU_INDEX->cols; c++){
        for(r = 0; r<DxtU_INDEX->rows -1 ; r++){
         /* U(:,1:end-1) */
        float elements_01 = cvmGet(DxtU_INDEX,r,c);
         /* U(:,2:end) */
        float elements_02 = cvmGet(DxtU_INDEX, r+1,c);
        float result = elements_01 - elements_02;
        cvmSet(DxtU_INDEX_02, r,c,result);

        }
    }
    CvMat* DxtU = cvCreateMat(DxtU_INDEX->rows,DxtU_INDEX->cols, CV_32F);
    for(c = 0; c<DxtU_INDEX->cols; c++){
        for(r = 0; r<DxtU_INDEX->rows -1 ; r++){
            if(r == 0){
                float elements = cvmGet(DxtU_INDEX_01, 0, c);
                cvmSet(DxtU, 0,c,elements);
            }
            else{
                float elements_02 = cvmGet(DxtU_INDEX_02, r-1, c);
                cvmSet(DxtU, r,c,elements_02);
            }
        }
    }

    /* matlab :  dytu = [U(end,:)-U(1, :); U(1:end-1,:)-U(2:end,:)];
                            DytU_INDEX_01       DytU_INDEX_02      */

    CvMat *DytU_INDEX_01 = cvCreateMat(DytU_INDEX->rows,1,CV_32F);
    for(r = 0; r<DytU_INDEX->rows; r++){
        /* U(end,:)*/
        float elements_01 = cvmGet(DytU_INDEX,r,DytU_INDEX->cols-1);
        /* U(1, :) */
        float elements_02 = cvmGet(DytU_INDEX,r,0);
        float result = elements_01 - elements_02;
        cvmSet(DytU_INDEX_01, r,0,result);
    }
    CvMat* DytU_INDEX_02 = cvCreateMat(DytU_INDEX->rows,DytU_INDEX->cols-1, CV_32F);
    for(c = 0; c<DxtU_INDEX->cols-1; c++){
        for(r = 0; r<DytU_INDEX->rows ; r++){
         /* U(1:end-1,:) */
        float elements_01 = cvmGet(DytU_INDEX,r,c);
         /* U(2:end,:) */
        float elements_02 = cvmGet(DytU_INDEX, r,c+1);
        float result = elements_01 - elements_02;
        cvmSet(DytU_INDEX_02, r,c,result);

        }
    }
    CvMat* DytU = cvCreateMat(DytU_INDEX->rows,DytU_INDEX->cols, CV_32F);
    for(c = 0; c<DytU_INDEX->cols-1; c++){
        for(r = 0; r<DytU_INDEX->rows; r++){
            if(c == 0){
                float elements = cvmGet(DytU_INDEX_01, r, 0);
                cvmSet(DytU, r,0,elements);
            }
            else{
                float elements_02 = cvmGet(DytU_INDEX_02, r, c-1);
                cvmSet(DytU, r,c,elements_02);
            }
        }
    }
    rhs = cvCreateMat(DxtU->rows, DxtU->cols, CV_32F);
    cvmAdd(DxtU,DytU,rhs);

    for(c = 0; c<rhs->cols; c++){
        for(r = 0; r<rhs->rows -1 ; r++){
            float elements = cvmGet(rhs,r,c);
            float result = tau * elements;
            cvmSet(rhs,r,c,result);
        }
     }


    cvReleaseMat(&DxtU_INDEX);
    cvReleaseMat(&DytU_INDEX);
    cvReleaseMat(&DxtU_INDEX_01);
    cvReleaseMat(&DxtU_INDEX_02);
    cvReleaseMat(&DytU_INDEX_01);
    cvReleaseMat(&DytU_INDEX_02);

}
