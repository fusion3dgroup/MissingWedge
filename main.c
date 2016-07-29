#include <stdio.h>
#include <stdlib.h>

#include <cv.h>
#include "extern/enhance.h"

int main()
{

    int m = 256, n = 256;

    //Read the .png image
    IplImage *image = cvLoadImage("sample.png");

    double m_aTV = 1e-4, m_aL1 = 1e-6, m_iter = 2, m_gamma = 1.0, m_beta = 10, m_tol = 1e-10;

    IplImage *dstImageU = cvCreateImage(cvGetSize(image), IPL_DEPTH_8U, 1);

   //IplImage to CvMat
    //CvMat matHeader, *F_Mat;
    //F_Mat = cvGetMat(image,&matHeader,0,0);

    enhance(image, n, m_iter, m_gamma, m_beta, m_tol, m_aTV, m_aL1, &dstImageU);




    cvReleaseImage(&image);
    cvReleaseImage(&dstImageU);

    return 0;
}
