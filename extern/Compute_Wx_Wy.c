#include <cv.h>
#include <complex.h>

/*
INPUT: Ux, Uy, bx, by, tau
OUTPUT: Wx, Wy

Equivalent to the following MATLAB code

    tic;
    UUx = Ux + bx; UUy = Uy + by;
    V = sqrt(UUx.*conj(UUx) + UUy.*conj(UUy));
    V = max(V - tau, 0) ./ max(V,eps);
    Wx = V.*UUx; Wy = V.*UUy;
    toc;

 */

void Compute_Wx_Wy(CvMat* Ux, CvMat* Uy, CvMat* bx, CvMat* by,double tau, CvMat* Wx, CvMat* Wy){

   /* matalb :  UUx = Ux + bx; UUy = Uy + by; */
    CvMat* UUx = cvCreateMat(Ux->rows, Ux->cols,CV_32F);
    CvMat* UUy = cvCreateMat(Uy->rows, Uy->cols,CV_32F);
    CvMat* V = cvCreateMat(UUx->rows, UUx->cols,CV_32F);

    Wx = cvCreateMat(V->rows, V->cols,CV_32F);
    Wy = cvCreateMat(V->rows, V->cols,CV_32F);

    cvmAdd(Ux, bx, UUx);
    cvmAdd(Uy, by, UUy);

    /* matlab : V = sqrt(UUx.*conj(UUx) + UUy.*conj(UUy)); */
    int c, r;
     for(c = 0; c<UUx->cols; c++){
        for(r = 0; r<UUx->rows; r++){
            float elements_uux = cvmGet(UUx, r, c);
            float elements_uuy = cvmGet(UUy, r, c);
            float result = sqrt(elements_uux * conj(elements_uux) + elements_uuy * conj(elements_uuy));
            cvmSet(V,r,c,result);
        }
    }
    /* matlab : V = max(V - tau, 0) ./ max(V,eps); */
    float eps = 2.2204e-16;
     for(c = 0; c<V->cols; c++){
        for(r = 0; r<V->rows; r++){
            float elements = cvmGet(V,r,c);
            float elements_a = 0, elements_b = 0;
            if(elements - tau < 0){ elements_a = 0;}
            else{elements_a = elements;}
            if(elements < eps) { elements_b = eps;}
            else{elements_b = elements;}
            float result = elements_a / elements_b;
            cvmSet(V,r,c, result);
        }
    }
    /* matlab  :   Wx = V.*UUx; Wy = V.*UUy;*/
    cvmMul(V, UUx, Wx);
    cvmMul(V, UUy, Wy);

    cvReleaseMat(&UUx);
    cvReleaseMat(&UUy);
    //cvReleaseMat(&UU);
    cvReleaseMat(&V);

}




