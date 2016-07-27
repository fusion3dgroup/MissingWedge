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

    // Ux = [diff(U,1,2), U(:,1) - U(:,n)];
    CvMat *diff_blank, *X1;
    int i,j;
    for(j = 1; j <=U->cols; j++){
        for(i = 1; i <= U->rows-1; i++){
            if(i <= U->rows-1){
            cvmSet(diff_blank, i, j, cvmGet(U,i+1,j) - cvmGet(U,i,j));
            }
            cvmSet(X1, i, j,  cvmGet(U, i, 1) - cvmGet(U, i, U->cols));
        }
     }

    Ux = cvCreateMat(diff_blank->rows + X1->rows, diff_blank->cols + X1->cols, CV_32FC1);
    for(j = 1; j<diff_blank->cols; j++){
        for(i = 1; i<diff_blank->rows; i++){
            cvmSet(Ux, i, j, cvmGet(diff_blank, i, j));
         }
    }
    for(j = diff_blank->cols; j<X1->cols; j++){
        for(i = diff_blank->rows; i<X1->rows; i++){
            cvmSet(Ux, i, j, cvmGet(X1, i, j));
         }
    }

    //         Uy = [diff(U,1,1); U(1,:)-U(m,:)];
    cvSetZero(diff_blank);
    cvSetZero(X1);

    for(j = 1; j <=U->cols-1; j++){
        for(i = 1; i <= U->rows; i++){
            if(j<U->cols-1){
            cvmSet(diff_blank, i, j, cvmGet(U,i,j+1) - cvmGet(U,i,j));
            }
            cvmSet(X1, i, j,  cvmGet(U, 1, j) - cvmGet(U, U->rows, j));
        }
     }

    Uy = cvCreateMat(diff_blank->rows + X1->rows, diff_blank->cols + X1->cols, CV_32FC1);
    for(j = 1; j<diff_blank->cols; j++){
        for(i = 1; i<diff_blank->rows; i++){
            cvmSet(Uy, i, j, cvmGet(diff_blank, i, j));
         }
    }
    for(j = diff_blank->cols; j<X1->cols; j++){
        for(i = diff_blank->rows; i<X1->rows; i++){
            cvmSet(Uy, i, j, cvmGet(X1, i, j));
         }
    }



    cvReleaseMat(&X1);
    cvReleaseMat(&diff_blank);

}

/* The gateway routine. */
#if 0
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    double *Ur, *Ui, *Uxr, *Uyr, *Uxi, *Uyi;
    double xr, yr, xi, yi;
    bool bComplex;
    mwSize rows, cols;
    mwSize i, j, pos;

    /* check for the proper number of arguments */
    if(nrhs != 1)
      mexErrMsgTxt("Exactly one input argument is required.");
    if(nlhs != 2)
      mexErrMsgTxt("Need two output arguments.");
    /*Check that both inputs are row vectors*/
    rows = mxGetM(prhs[0]);
    cols = mxGetN(prhs[0]);
    bComplex = mxIsComplex(prhs[0]);

    /* get pointers to the real and imaginary parts of the inputs */
    Ur = mxGetPr(prhs[0]);
    Ui = mxGetPi(prhs[0]);

    /* create Ux and Uy */
    if (bComplex)
    {
        plhs[0] = mxCreateDoubleMatrix(rows, cols, mxCOMPLEX);
        plhs[1] = mxCreateDoubleMatrix(rows, cols, mxCOMPLEX);
    }
    else
    {
        plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
        plhs[1] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    }

    Uxr = mxGetPr(plhs[0]); Uyr = mxGetPr(plhs[1]); // Vr = mxGetPr(plhs[2]);
    if (bComplex)
        Uxi = mxGetPi(plhs[0]); Uyi = mxGetPi(plhs[1]);

    /* MAIN COMPUTATION */
    if (bComplex)
    {
        /* Complex COMPUTATION */
        pos = 0;
        for (j = 0; j < cols-1; j++)
        {
            for (i = 0; i < rows-1; i++)
            {
                xr = Ur[pos+rows]-Ur[pos];
                xi = Ui[pos+rows]-Ui[pos];
                yr = Ur[pos+1]-Ur[pos];
                yi = Ui[pos+1]-Ui[pos];
                Uxr[pos] = xr; Uxi[pos] = xi;
                Uyr[pos] = yr; Uyi[pos] = yi;
                pos++;
            }
            /* for i = rows-1 */
            {
                xr = Ur[pos+rows]-Ur[pos];
                xi = Ui[pos+rows]-Ui[pos];
                yr = Ur[pos-rows+1]-Ur[pos];
                yi = Ui[pos-rows+1]-Ui[pos];
                Uxr[pos] = xr; Uxi[pos] = xi;
                Uyr[pos] = yr; Uyi[pos] = yi;
                pos++;
            }
        }
        /* for j = cols - 1 */
            for (i = 0; i < rows-1; i++)
            {
                xr = Ur[i]-Ur[pos];
                xi = Ui[i]-Ui[pos];
                yr = Ur[pos+1]-Ur[pos];
                yi = Ui[pos+1]-Ui[pos];
                Uxr[pos] = xr; Uxi[pos] = xi;
                Uyr[pos] = yr; Uyi[pos] = yi;
                pos++;
            }
            /* for i = rows-1 */
            {
                xr = Ur[rows-1]-Ur[pos];
                xi = Ui[rows-1]-Ui[pos];
                yr = Ur[pos-rows+1]-Ur[pos];
                yi = Ui[pos-rows+1]-Ui[pos];
                Uxr[pos] = xr; Uxi[pos] = xi;
                Uyr[pos] = yr; Uyi[pos] = yi;
            }
    }
    else
    {
    /* REAL-only COMPUTATION */
        pos = 0;
        for (j = 0; j < cols-1; j++)
        {
            for (i = 0; i < rows-1; i++)
            {
                xr = Ur[pos+rows]-Ur[pos];
                yr = Ur[pos+1]-Ur[pos];
                Uxr[pos] = xr;
                Uyr[pos] = yr;
                pos++;
            }
            /* for i = rows-1 */
            {
                xr = Ur[pos+rows]-Ur[pos];
                yr = Ur[pos-rows+1]-Ur[pos];
                Uxr[pos] = xr;
                Uyr[pos] = yr;
                pos++;
            }
        }
        /* for j = cols - 1 */
            for (i = 0; i < rows-1; i++)
            {
                xr = Ur[i]-Ur[pos];
                yr = Ur[pos+1]-Ur[pos];
                Uxr[pos] = xr;
                Uyr[pos] = yr;
                pos++;
            }
            /* for i = rows-1 */
            {
                xr = Ur[rows-1]-Ur[pos];
                yr = Ur[pos-rows+1]-Ur[pos];
                Uxr[pos] = xr;
                Uyr[pos] = yr;
            }
    }

    return;
}
#endif // 0




