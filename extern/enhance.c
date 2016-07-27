/* Include files */
#include "enhance.h"


void enhance(IplImage *image, int n, double m_iter, double m_gamma, double m_beta,
                    double m_tol, double m_aTV, double m_aL1, CvMat* U)
{
    struct OPTS opts_para;

  /*
    matlab code : constraint = 1;  x = reshape(I,n*n,1); L   = 199;
                  [M,Mh,mi,mhi] = LineMask_limitedangleRange(L,n);
  */
  IplImage *x = cvCreateImage(cvSize(n,n), 8, 1);
  cvResize(image,x, 0);

   int constraint = 1;
   int L = 199;

   CvMat * M, *Mh, *mi, *mhi;

   LineMask_limitedandleTange(L,n, &M, &Mh, &mi, &mhi);

    CvMat *OMEGA = cvCreateMat(mhi->rows,mhi->cols,CV_32FC1);
    cvmCopy(mhi, OMEGA);

    int k;
    if(OMEGA->rows > OMEGA->cols){   //number of vector element?
        k= OMEGA->rows;
    }
    else k= OMEGA->cols;

    CvMat *b;
    A_fhp(x, OMEGA, &b);

    CvMat *picks = cvCreateMat(OMEGA->rows, OMEGA->cols, CV_32FC1);
    cvmCopy(OMEGA,picks);

    int cM, cN, p, o;
    IplImage *FB = cvCreateImage(cvGetSize(image), IPL_DEPTH_8U, 1);
    IplImage *F = cvCreateImage(cvGetSize(image), IPL_DEPTH_8U, 1);
    cvCopy(image,F,0);

    cM = F->width;
    cN = F->height;

    cvFFT(F,F, 0,0);
    for(p = 0; p<F->width; p++){
        for(o=0; o<F->height; o++){
            FB->imageData[o+p*FB->widthStep] = F->imageData[o+p*F->widthStep] * sqrt(cM*cN);
        }
    }

   //IplImage to CvMat
    CvMat matHader, *FB_Mat;
    FB_Mat = cvGetMat(FB,&matHader,0,0);

    //B=FB(Picks); edit
    //WT = []; W=[];

    double aTV = m_aTV;
    double aL1 = m_aL1;

    CvMat *wav, *W, *WT;
    daubcqf(4, &wav);

    midwt(x,wav,&W);
    mdwt(x,wav, &WT);

    opts_para.maxItr = m_iter;
    opts_para.gamma = m_gamma;
    opts_para.beta = m_beta;
    opts_para.relchg_tol = 1;
    opts_para.normalize = 1;

    CvMat * pick ,*B;
    cvSetZero(pick);        //pick = false(m,n);
    //pick(picks)=true;

    double min, max;
    cvMinMaxLoc(image,&min,&max, 0,0,0);

    double range = max -min; //difference between the maximun and the minimum of a I
    int m = 256;

    RecPF_constraint(m,n,aTV, aL1,pick,B,2, opts_para, WT, W,range,image,constraint, &U);


    cvReleaseMat(&picks);
    cvReleaseMat(&FB_Mat);
    cvReleaseMat(&wav);
    cvReleaseMat(&W);
    cvReleaseMat(&WT);
    cvReleaseImage(&FB);
    cvReleaseImage(&F);

  }

