/* Include files */
#include "enhance.h"
#include <math.h>
#define pi 3.141592

#define round_ban(a) (floor(a)*pow(10,0)+0.5)/pow(10,0)

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
   int N = n;

   //CvMat * M, *Mh, *mi,*mhi;
   //LineMask_limitedandleTange(L,n, &M, &Mh, &mi, &mhi);

    #if 1
    //CvMat *mhi = cvCreateMat(17974,1,CV_32F);
    //int i;
    //for(i=0; i<mhi->rows; i++){
    //    float elements = i+194;
    //    cvmSet(mhi, i, 0, elements);
    //}
        /* matlab : thc = linspace(0, pi-pi/L, L); */

    CvMat *thc = cvCreateMat(1,L,CV_32F);

    int i,j;
    cvmSet(thc,0,0,0.);

    for(i = 1; i< L; i++){
        float elements = i*(pi-pi/L)/(L-1);
        cvmSet(thc,0,i,elements);
    }

    /*matlab : M = zeros(N);*/
    CvMat *M = cvCreateMat(N,N,CV_32F);
    cvmSetZero(M);

    /* matlab : (-N/2+1 : N/2-1) */
    int ll = 0;
    CvMat *blank_aa = cvCreateMat(1, N, CV_32F);
    CvMat *yr = cvCreateMat(1,N,CV_32F);
    CvMat *xc = cvCreateMat(1,N,CV_32F);

    for(ll = 0; ll < N; ll++){
        float elements= N*(-1.0)/2.0 + ll*1.0 + 1.0;
        cvmSet(blank_aa, 0, ll, elements);
    }

    /* matlab : for loop */
    int nn=0;


    for(ll = 0; ll < L; ll++){
        if(cvmGet(thc,0,ll)>= pi/3 && cvmGet(thc,0,ll)<= (pi/2+pi/18) ||
            cvmGet(thc,0,ll) >=(2*pi/3-pi/18) && cvmGet(thc,0,ll)<= 2*pi/3){

            /* matlab : xc = round(cot(thc(ll))*(-N/2+1:N/2-1))+N/2+1; */
            float element_B = cvmGet(thc,0,ll);
                  element_B = tan(element_B);
                //cotangent(0)..
                if(element_B == 0) element_B = 0.1;
                   element_B = 1.0/element_B;

            for(i=0; i<blank_aa->cols-1; i++){
                int element_C = round_ban(element_B * cvmGet(blank_aa,0, i))+N/2+1;
                cvmSet(xc,0,i,element_C);
            }
            for(nn = 0; nn<N-1; nn++){
                /* matlab : M(nn+1,xc(nn)) = 1;*/
                int elements_A = cvmGet(xc,0,nn);
                cvmSet(M,nn+1,elements_A,1.0);
            }

        }//if
        else{
            if(cvmGet(thc,0,ll) <= pi/4 || cvmGet(thc,0,ll)>3*pi/4){
                /* matlab : yr = round(tan(thc(ll))*(-N/2+1:N/2-1))+N/2+1; */
                float element_B = cvmGet(thc,0,ll);
                      element_B = tan(element_B);

                for(i=0; i<blank_aa->cols-1; i++){
                    int element_C = round_ban(element_B * cvmGet(blank_aa,0, i))+N/2+1;
                    cvmSet(yr,0,i,element_C);
                }
                for(nn = 0; nn<N-1; nn++){
                    /* matlab : M(yr(nn), nn+1) = 1;*/
                    int elements_A = cvmGet(yr,0,nn);
                    cvmSet(M,elements_A,nn+1,1.0);
                }
            }
            else{
                /* matlab : xc = round(cot(thc(ll))*(-N/2+1:N/2-1))+N/2+1; */
                float element_B = cvmGet(thc,0,ll);
                      element_B = tan(element_B);
                    //cotangent(0)..
                    if(element_B == 0) element_B = 0.1;
                       element_B = 1.0/element_B;

                for(i=0; i<blank_aa->cols-1; i++){
                    int element_C = round_ban(element_B * cvmGet(blank_aa,0, i))+N/2+1;
                    cvmSet(xc,0,i,element_C);
                }
                for(nn = 0; nn<N-1; nn++){
                    /* matlab : M(nn+1,xc(nn)) = 1;*/
                    int elements_A = cvmGet(xc,0,nn);
                    cvmSet(M,nn+1,elements_A,1.0);
                }
            }
        }//else
    }//for


    //upper half plane mask
    CvMat *Mh = cvCreateMat(N,N,CV_32FC1);
    cvmSetZero(Mh);
    cvmCopy(M,Mh);

    /* matlab : Mh =(N/2+2:N,:) = 0; */
    int r,c;
    for(c = 0; c<N; c++){
        for(r = N/2+2; r<N; r++){
            cvmSet(Mh, r, c, 0);
        }
    }
    /* matlab : Mh =(N/2+1,N/2+1:N) = 0; */
    for(c = N/2; c<N; c++){
        cvmSet(Mh, N/2+1, c, 0);
    }

    /*M = ifftshift(M); */
    //cvFFT(M,M,1,1); //0:fft, 1:ifft
    //cvifftshift(M);

    /* matlab : mMH = find(M)*/
    int count = 0;
    for(c = 0; c< M->cols; c++){
        for(r = 0; r< M->rows; r++){
            if(cvmGet(M,r,c) != 0){
                count=count+1;
            }
        }
    }
    int a= count;
    CvMat *mi = cvCreateMat(a+1,1,CV_32F);
    count = 0;
    for(c = 0; c< M->cols; c++){
        for(r = 0; r< M->rows; r++){
            if(cvmGet(M,r,c) != 0){
                float elements = cvmGet(M,r,c);
                cvmSet(mi,count,0, elements);    //find the nonzero elements in matrix
                count=count+1;
            }
        }
    }
    /*Mh = ifftshift(Mh); */
    //cvifftshift(Mh,Mh,1,1);

    /* matlab : mMH = find(M)*/
    count = 0;
    for(c = 0; c< Mh->cols; c++){
        for(r = 0; r< Mh->rows; r++){
            if(cvmGet(Mh,r,c) != 0){
                count=count+1;
            }
        }
    }

    a= count;
    CvMat *mhi = cvCreateMat(a+1,1,CV_32F);
    count = 0;
    for(c = 0; c< Mh->cols; c++){
        for(r = 0; r< Mh->rows; r++){
            if(cvmGet(Mh,r,c) != 0){
                cvmSet(mhi,count, 0, cvmGet(Mh,r,c));    //find the nonzero elements in matrix
                count=count+1;
            }
        }
    }

    #endif // 0

    CvMat *OMEGA = cvCreateMat(mhi->rows+1,1,CV_32F);
    /* matlab : OMEGA = mhi*/
   // cvmCopy(mhi, OMEGA);

    /* k=length(OMEGA);*/
    //int k = OMEGA->rows;

    /* matlab : OMEGA = [1;OMEGA] */
    cvmSet(OMEGA, 0, 0, 1.);
    for(i=0; i<mhi->rows; i++){
        float elements = cvmGet(mhi, i,0);
        cvmSet(OMEGA, i+1, 0, elements);
    }


    CvMat *b;
    //A_fhp(x, OMEGA, &b);

    CvMat *picks = cvCreateMat(OMEGA->rows, OMEGA->cols, CV_32FC1);
    cvmCopy(OMEGA,picks);

    /* matalb F = I; */
    IplImage *F = cvCreateImage(cvGetSize(image), IPL_DEPTH_8U, 1);
    cvCopy(image,F,0);


   //IplImage to CvMat
    CvMat *F_Mat = cvCreateMat(F->width, F->height, CV_32F);
//    int c,r;
    for(c = 0; c<F->height; c++){
        for(r=0; r<F->width; r++){
            float elemnets = cvGetReal2D(F,r,c);
            cvmSet(F_Mat, c,r, elemnets);
        }
    }

    CvMat *FB = cvCreateMat(F_Mat->rows, F_Mat->cols, CV_32F);

    int m = F_Mat->cols;
    n = F_Mat->rows;

    /* matlab : FB = fft2(F)/sqrt(m*n); */
    cvFFT(F_Mat, F_Mat, CV_DXT_FORWARD ,0);


    for(c = 0; c<F_Mat->cols; c++){
        for(r=0; r<F_Mat->rows; r++){
            float elements= cvmGet(F_Mat,r,c)/sqrt(m*n);
            cvmSet(FB,r,c,elements);
        }
    }


    /* matlab : B=FB(Picks); */
    CvMat *B = cvCreateMat(picks->rows, picks->cols, CV_32F);

    for(r = 0; r< picks->rows; r++){
        int pick_blank = cvmGet(picks,r,0);
        int quotient = pick_blank / FB->cols;
        int remainder = pick_blank % FB->cols;
        float elements = cvmGet(FB,quotient,remainder);
        cvmSet(B, r, 0, elements);
    }


    /* matlab : WT = []; W=[]; */
     CvMat *wav, *W, *WT;

     /*matlab : aTV = m_aTV; aL1 = m_aL1; */
    double aTV = m_aTV;
    double aL1 = m_aL1;

    //daubcqf(4, &wav);
    wav = cvCreateMat(1,8,CV_32F);
    cvmSet(wav,0,0,0.4830);
    cvmSet(wav,0,1,0.8365);
    cvmSet(wav,0,2,0.2241);
    cvmSet(wav,0,3,-0.1294);
    cvmSet(wav,0,4,0.1294);
    cvmSet(wav,0,5,0.2241);
    cvmSet(wav,0,6,-0.8365);
    cvmSet(wav,0,7,0.4830);


    //midwt(x,wav,&W);
    midwt(x,wav);
    //mdwt(x,wav, &WT);
    mdwt(x,wav);

    opts_para.maxItr = m_iter;
    opts_para.gamma = m_gamma;
    opts_para.beta = m_beta;
    opts_para.relchg_tol = 1;
    opts_para.normalize = 1;


    /* matlab : pick = false(m,n); */
    CvMat * pick = cvCreateMat(m, n, CV_32F);
    cvSetZero(pick);        //
    /* matlab : pick(picks)=true; To Logical*/

    /* matlab : range(I(:)); */
    double min, max;
    cvMinMaxLoc(image,&min,&max, 0,0,0);

    double range = max -min; //difference between the maximun and the minimum of a I

    RecPF_constraint(m,n,aTV, aL1,pick,B,2, opts_para, WT, W,range,image,constraint, &U);


    cvReleaseMat(&picks);
    cvReleaseMat(&FB);
    cvReleaseMat(&wav);
    cvReleaseMat(&W);
    cvReleaseMat(&WT);
    cvReleaseImage(&FB);
    cvReleaseImage(&F);

    cvReleaseMat(&thc);
    cvReleaseMat(&xc);
    cvReleaseMat(&yr);
    cvReleaseMat(&blank_aa);

  }
