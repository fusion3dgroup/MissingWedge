#include <stdio.h>
#include <stdlib.h>

#include <cv.h>
#include "extern/enhance.h"

int main()
{

    int m = 256, n = 256;

    //Read the .png image
    IplImage *image = cvLoadImage("sample.png");

    double m_aTV = 1e-4, m_aL1 = 1e-6, m_iter = 500, m_gamma = 1.0, m_beta = 10, m_tol = 1e-10;

    IplImage *dstImageU = cvCreateImage(cvGetSize(image), IPL_DEPTH_8U, 1);
    enhance(image, n, m_iter, m_gamma, m_beta, m_tol, m_aTV, m_aL1, &dstImageU);


    cvShowImage("image", image);
    cvWaitKey(0);

    cvReleaseImage(&image);
    cvReleaseImage(&dstImageU);
    cvDestroyWindow("image");

    return 0;
}
