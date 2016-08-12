/* Include files */
#include "enhance.h"
//#include "LineMask_limitedangleRange.h"
#include <math.h>

#define pi 3.14159265358979323846

#define round_ban(a) (floor(a)*pow(10,0)+0.5)/pow(10,0)
#define SQUARE(x) x*x

void enhance(IplImage *image, int n, double m_iter, double m_gamma, double m_beta,
                    double m_tol, double m_aTV, double m_aL1, CvMat* U)
{
    struct OPTS opts_para;

   /*
    matlab code : constraint = 1;  x = reshape(I,n*n,1); L   = 199;
                  [M,Mh,mi,mhi] = LineMask_limitedangleRange(L,n);
  */
  IplImage *x_ipl = cvCreateImage(cvSize(n,n), 8, 1);
  cvResize(image,x_ipl, 0);

   int constraint = 1;
   int L = 199;
   int N = n;
   int r,c;

    CvMat *x_mat = cvCreateMat(x_ipl->height,x_ipl->width,CV_32F);
   /* matlab : org = double(i)*/
    for(c = 0; c<x_ipl->height; c++){
        for(r = 0; r<x_ipl->width; r++){
            float elements = cvGetReal2D(x_ipl,r,c) / 255.0;
            cvmSet(x_mat,r,c,elements);
        }
    }
    /* matlab : x = reshape(I,n*n,1); */
    CvMat *x = cvCreateMat(n*n,1,CV_32F);
    for(c = 0; c<x_mat->cols; c++){
        for(r = 0; r<x_mat->rows; r++){
                float elements = cvmGet(x_mat,r,c);
                int aa = r;
                aa = aa + c*x_mat->rows;
                cvmSet(x,aa,0,elements);
        }
    }

    #if 1
    //LineMask_limitedandleTange(L,n, mhi, M, Mh, mi);
    /* matlab : thc = linspace(0, pi-pi/L, L); */

    CvMat *thc = cvCreateMat(1,L,CV_32F);

    int i,j;
    cvmSet(thc,0,0,0.0);

    for(i = 1; i< L; i++){
        float elements = i*(pi-pi/L)/(L-1);
        cvmSet(thc,0,i,elements);
    }


    /*matlab : M = zeros(N);*/
    CvMat *M = cvCreateMat(N,N,CV_32F);
    cvmSetZero(M);

    /* matlab : (-N/2+1 : N/2-1) */
    int ll = 0;
    CvMat *blank_aa = cvCreateMat(1, N-1, CV_32F);
    CvMat *yr = cvCreateMat(1,N-1,CV_32F);
    CvMat *xc = cvCreateMat(1,N-1,CV_32F);


    for(ll = 0; ll < blank_aa->cols; ll++){
        float elements= N*(-1.0)/2.0+1.0 + ll*1.0;
        cvmSet(blank_aa, 0, ll, elements);
    }

    /* matlab : for loop */
    int nn=0;


    for(ll = 0; ll < L; ll++){
        if(((cvmGet(thc,0,ll) >= pi/3) && (cvmGet(thc,0,ll)<= (pi/3+pi/18))) ||
            ((cvmGet(thc,0,ll) >= (2*pi/3-pi/18)) && (cvmGet(thc,0,ll)<= 2*pi/3))){

            /* matlab : xc = round(cot(thc(ll))*(-N/2+1:N/2-1))+N/2+1; */
            float element_A = cvmGet(thc,0,ll);
            float element_B = 1.0/tan(element_A);
                //cotangent(0)..
                //if(element_B == 0) {element_B = 0.0;}
            for(i=0; i<blank_aa->cols; i++){
                int element_C = round(element_B * cvmGet(blank_aa,0, i))+N/2+1;
                cvmSet(xc,0,i,element_C);
            }
            for(nn = 0; nn<N-1; nn++){
                /* matlab : M(nn+1,xc(nn)) = 1;*/
                int elements_D = cvmGet(xc,0,nn);
                cvmSet(M,nn+1,elements_D-1,0.0);
            }

        }//if

        else{
            if((cvmGet(thc,0,ll) <= pi/4) || (cvmGet(thc,0,ll) > 3*pi/4)){
                /* matlab : yr = round(tan(thc(ll))*(-N/2+1:N/2-1))+N/2+1; */
                float element_A = cvmGet(thc,0,ll);
                float element_B = tan(element_A);

                for(i=0; i<blank_aa->cols; i++){
                    int element_C = round(element_B * cvmGet(blank_aa,0, i))+N/2+1;
                    cvmSet(yr,0,i,element_C);
                }
                for(nn = 0; nn<N-1; nn++){
                    /* matlab : M(yr(nn), nn+1) = 1;*/
                    int elements_D = cvmGet(yr,0,nn);
                    cvmSet(M,elements_D-1,nn+1,1.0);
                }
             }

             else{
                /* matlab : xc = round(cot(thc(ll))*(-N/2+1:N/2-1))+N/2+1; */
                float element_A = cvmGet(thc,0,ll);
                float  element_B = 1.0/tan(element_A);
                    //cotangent(0)..
                    //if(element_B == 0) {element_B = 0.0;}
                     // else{ element_B = 1.0/element_B;}

                for(i=0; i<blank_aa->cols; i++){
                    int element_C = round(element_B * cvmGet(blank_aa,0, i))+N/2+1;
                    cvmSet(xc,0,i,element_C);
                }
                for(nn = 0; nn<N-1; nn++){
                    /* matlab : M(nn+1,xc(nn)) = 1;*/
                    int elements_D = cvmGet(xc,0,nn);
                    cvmSet(M,nn+1,elements_D-1,1.0);
                }
             }
        }//else
    }//for


    cvReleaseMat(&yr);
    cvReleaseMat(&xc);
    cvReleaseMat(&blank_aa);
    cvReleaseMat(&thc);

    //upper half plane mask
    CvMat *Mh = cvCreateMat(N,N,CV_32FC1);
    cvmSetZero(Mh);
    cvmCopy(M,Mh);



    /* matlab : Mh(N/2+2:N,:) = 0; */
    for(c = 0; c<N; c++){
        for(r = N/2+1; r<N; r++){
            cvmSet(Mh, r, c, 0.0);
        }
    }
    /* matlab : Mh(N/2+1,N/2+1:N) = 0; */
    for(c = N/2; c<N; c++){
        cvmSet(Mh, N/2, c, 0.0);
    }

    /*M = ifftshift(M); */
    CvMat *M_blank = cvCreateMat(Mh->rows, Mh->cols, CV_32F);

    for(c = 0; c< M->cols; c++){
        for(r = 0; r< M->rows; r++){
            if(c < M->cols/2 && r < M->rows/2){ // first quadrant
                float elements = cvmGet(M, r, c);
                cvmSet(M_blank,r+M->rows/2,c+M->cols/2,elements);
            }
            if(c >= M->cols/2 && r < M->rows/2){    //second
                float elements = cvmGet(M, r, c);
                cvmSet(M_blank,r+M->rows/2,c-M->cols/2,elements);
            }
            if(c < M->cols/2 && r >= M->rows/2){    //third
                float elements = cvmGet(M, r, c);
                cvmSet(M_blank,r-M->rows/2,c+M->cols/2,elements);
            }
            if(c >= M->cols/2 && r >= M->rows/2){    //third
                float elements = cvmGet(M, r, c);
                cvmSet(M_blank,r-M->rows/2,c-M->cols/2,elements);
            }
        }
    }
    /* matlab : mi = find(M)*/
    int count = 0;
    for(c = 0; c< M_blank->cols; c++){
        for(r = 0; r< M_blank->rows; r++){
            if(cvmGet(M_blank,r,c) != 0){
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
                float elements = c * Mh->rows + r+1;
                cvmSet(mi,count,0, elements);    //find the nonzero elements in matrix
                count=count+1;
            }
        }
    }

    /*Mh = ifftshift(Mh); */
    CvMat *Mh_blank = cvCreateMat(Mh->rows, Mh->cols, CV_32F);

    for(c = 0; c< Mh->cols; c++){
        for(r = 0; r< Mh->rows; r++){
            if(c < Mh->cols/2 && r < Mh->rows/2){ // first quadrant
                float elements = cvmGet(Mh, r, c);
                cvmSet(Mh_blank,r+Mh->rows/2,c+Mh->cols/2,elements);
            }
            if(c >= Mh->cols/2 && r < Mh->rows/2){    //second
                float elements = cvmGet(Mh, r, c);
                cvmSet(Mh_blank,r+Mh->rows/2,c-Mh->cols/2,elements);
            }
            if(c < Mh->cols/2 && r >= Mh->rows/2){    //third
                float elements = cvmGet(Mh, r, c);
                cvmSet(Mh_blank,r-Mh->rows/2,c+Mh->cols/2,elements);
            }
            if(c >= Mh->cols/2 && r >= Mh->rows/2){    //third
                float elements = cvmGet(Mh, r, c);
                cvmSet(Mh_blank,r-Mh->rows/2,c-Mh->cols/2,elements);
            }
        }
    }

    /* matlab : mhi = find(Mh)*/
    count = 0;
    for(c = 0; c< Mh_blank->cols; c++){
        for(r = 0; r< Mh_blank->rows; r++){
            if(cvmGet(Mh_blank,r,c) != 0){
                count=count+1;
            }
        }
    }
    a= count;
    CvMat *mhi = cvCreateMat(a,1,CV_32F);
    count = 0;
    for(c = 0; c< Mh_blank->cols; c++){
        for(r = 0; r< Mh_blank->rows; r++){
            if(cvmGet(Mh_blank,r,c) != 0){
                float elements = c * Mh_blank->rows + r+1;
                cvmSet(mhi,count, 0, elements);    //find the nonzero elements in matrix
                count=count+1;
            }
        }
    }
    cvReleaseMat(&M);
    cvReleaseMat(&Mh);
    cvReleaseMat(&M_blank);
    cvReleaseMat(&Mh_blank);

    #endif // 0

    CvMat *OMEGA = cvCreateMat(mhi->rows+1,1,CV_32F);

    /* k=length(OMEGA);*/
    //int k = OMEGA->rows;

    /* matlab : OMEGA = mhi
                OMEGA = [1;OMEGA] */
    cvmSet(OMEGA, 0, 0, 1.);
    for(r=0; r<mhi->rows; r++){
        float elements = cvmGet(mhi, r,0);
        cvmSet(OMEGA, r+1, 0, elements);
    }
    cvReleaseMat(&mi);
    cvReleaseMat(&mhi);
    //CvMat *b;
    //A_fhp(x, OMEGA, &b);

    CvMat *picks = cvCreateMat(OMEGA->rows, OMEGA->cols, CV_32F);
    cvmCopy(OMEGA,picks);

   cvReleaseMat(&OMEGA);

   /* matlab : F = I */
    CvMat *F_Mat = cvCreateMat(x_mat->rows, x_mat->cols, CV_32F);
    cvmCopy(x_mat,F_Mat);

    int m = F_Mat->rows;
    n = F_Mat->cols;

    /* matlab : FB = fft2(F)/sqrt(m*n); */
    CvMat *FB = cvCreateMat(F_Mat->rows, F_Mat->cols, CV_32F);
    CvMat* FB_imag =cvCreateMat(F_Mat->rows,F_Mat->cols,CV_32FC1);
    CvMat *imaginaryInput = cvCreateMat(F_Mat->rows,F_Mat->cols,CV_32F);
    CvMat *F_Mat_real = cvCreateMat(F_Mat->rows,F_Mat->cols,CV_32FC1);
    CvMat *F_Mat_FFT = cvCreateMat(F_Mat->rows,F_Mat->cols,CV_32FC2);

    cvmSetZero(imaginaryInput);
    cvMerge(F_Mat, imaginaryInput,NULL, NULL, F_Mat_FFT);

    cvDFT(F_Mat_FFT, F_Mat_FFT, CV_DXT_FORWARD ,0);
    cvSplit(F_Mat_FFT,F_Mat_real,imaginaryInput, NULL, NULL);

    cvReleaseMat(&F_Mat_FFT);


    for(c = 0; c<F_Mat->cols; c++){
        for(r=0; r<F_Mat->rows; r++){
            float elements_real= cvmGet(F_Mat_real,r,c)/sqrt(m*n);
            float elements_imag = cvmGet(imaginaryInput,r,c)/sqrt(m*n);
            cvmSet(FB,r,c,elements_real);
            cvmSet(FB_imag,r,c,elements_imag);
        }
    }
    cvReleaseMat(&F_Mat);
    cvReleaseMat(&imaginaryInput);
    cvReleaseMat(&F_Mat_real);


    /* matlab : B=FB(Picks); */
    CvMat *B = cvCreateMat(picks->rows, picks->cols, CV_32F);
    CvMat *B_imag = cvCreateMat(picks->rows, picks->cols, CV_32F);
    for(r = 0; r< picks->rows; r++){
        int pick_blank = cvmGet(picks,r,0) -1 ;
        int quotient = pick_blank/ FB->cols;
        int remainder = pick_blank % FB->cols;
        float elements = cvmGet(FB,remainder,quotient);
        float elements_imag = cvmGet(FB_imag,remainder,quotient);
        cvmSet(B, r, 0, elements);
        cvmSet(B_imag,r,0,elements_imag);
    }
    cvReleaseMat(&FB);
    cvReleaseMat(&FB_imag);
    /* matlab : WT = []; W=[]; */

     /*matlab : aTV = m_aTV; aL1 = m_aL1; */
    double aTV = m_aTV;
    double aL1 = m_aL1;

    //daubcqf(4, &wav);
    CvMat *wav = cvCreateMat(1,4,CV_32F);
    cvmSet(wav,0,0,0.4830);
    cvmSet(wav,0,1,0.8365);
    cvmSet(wav,0,2,0.2241);
    cvmSet(wav,0,3,-0.1294);
    //cvmSet(wav,0,4,0.1294);
    //cvmSet(wav,0,5,0.2241);
    //cvmSet(wav,0,6,-0.8365);
    //cvmSet(wav,0,7,0.4830);


    CvMat *W_mat = cvCreateMat(x->rows, x->cols, CV_32F);
    CvMat *WT_mat = cvCreateMat(x->rows, x->cols, CV_32F);

    opts_para.maxItr = m_iter;
    opts_para.gamma = m_gamma;
    opts_para.beta = m_beta;
    opts_para.relchg_tol = m_tol;
    opts_para.real_sol = 1;
    opts_para.normalize = 1;


    /* matlab : pick = false(m,n); */
    CvMat * pick = cvCreateMat(m, n, CV_32F);
    cvSetZero(pick);

    /* matlab : pick(picks)=true; To Logical*/
    for(r = 0; r<picks->rows; r++){
        int elements = cvmGet(picks, r, 0) - 1;
        int quotient = elements / pick->rows;
        int remainder = elements % pick->rows;
        cvmSet(pick, remainder,quotient, 1);
    }
    cvReleaseMat(&picks);
   /* matlab : range(I(:)); */
    double min, max;
    cvMinMaxLoc(x_mat,&min,&max, 0,0,0);

    int range = max -min; //difference between the maximun and the minimum of a I



    RecPF_constraint(m,n,aTV, aL1,pick,B,B_imag,2, opts_para, WT_mat, W_mat,range,x_mat,constraint,wav);

    cvReleaseImage(&x_ipl);
    cvReleaseMat(&x_mat);
    cvReleaseMat(&x);

    cvReleaseMat(&F_Mat);

    cvReleaseMat(&B);
    cvReleaseMat(&B_imag);

    cvReleaseMat(&wav);

    cvReleaseMat(&W_mat);
    cvReleaseMat(&WT_mat);
    cvReleaseMat(&pick);


  }
