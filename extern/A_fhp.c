
#include <cv.h>
#include <math.h>


int round_y(double a)
{
    int round_y;

    if(a < 0)  round_y = (int)(a * 100 + 0.5);
    else   round_y  = (int)(a * 100 + -0.5);

}
void A_fhp(IplImage *x, CvMat *OMEGA, CvMat* y)
{
    int w = x->width;
    int h = x->height;
    int n,p,o,i,j;

    if (w > h){
        n = round_y(sqrt(w));
    }
    else n = round_y(sqrt(h));

    IplImage *x_dst = cvCreateImage(cvSize(n,n), 8,1);
    IplImage *yc = cvCreateImage(cvGetSize(x_dst), 8,1);

    cvResize(x,x_dst,0);
    cvDFT(x_dst,x_dst,CV_DXT_FORWARD, 0);

    for(p=1; p< x_dst->width;p++){
        for(o=1; o< x_dst->height;o++){
           yc->imageData[o+p*yc->widthStep] = 1/n * x_dst->imageData[o+p*x_dst->widthStep];

        }
    }

   //IplImage to CvMat
    CvMat matHader, *pSrcMat;
    pSrcMat = cvGetMat(yc,&matHader,0,0);

    /*
        y = [yc(1,1); sqrt(2)*real(yc(OMEGA)) sqrt(2)*imag(yc(OMEGA))];
    */
    cvmSet(y,1,1,cvmGet(yc,1,1));

    for(i=1; i< OMEGA->cols;i++){
        for(j=1; j< OMEGA->rows;j++){
            double element = cvmGet(OMEGA,i,j);
            double c=element/i+j ;
            double elemet_2 = cvmGet(pSrcMat,i,c) * sqrt(2); /*real*/

            cvmSet(y,1,2,elemet_2);
        }
    }

    for(i=1; i< OMEGA->cols;i++){
        for(j=1; j< OMEGA->rows;j++){
            double element = cvmGet(OMEGA,i,j);
            double c=element/i+j ;
            double elemet_2 = cvmGet(pSrcMat,i,c) * sqrt(2); /*imag*/

            cvmSet(y,1,3,elemet_2);
        }
    }


    cvReleaseImage(&x_dst);
    cvReleaseImage(&yc);
    cvReleaseMat(&pSrcMat);
}
