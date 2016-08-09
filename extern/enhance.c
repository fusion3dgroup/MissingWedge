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
/* debug ------------------------------------------------------------*/
#if 0
    FILE * pf;
    pf = fopen("out_dst.txt", "w");
    for(r = 0; r<image->height; r++){
        for(c = 0; c<image->width; c++){
            float elements = cvGetReal2D(image,r,c);
            fprintf(pf, "%f ", elements,0);
            if(c == image->height-1){
                fprintf(pf,"\n ",0);
            }
        }
    }
            int col_size = image->height;
            int row_size = image->width;
    fclose(pf);
#endif // 0
/*---------------------------------------------------------------*/

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
/* debug ------------------------------------------------------------*/
#if 0
    FILE * pf;
    pf = fopen("out_dst.txt", "w");
    for(r = 0; r<x_ipl->height; r++){
        for(c = 0; c<x_ipl->width; c++){
            float elements = (float)cvGetReal2D(x_ipl,r,c);
            //printf("%f ",elements);
            fprintf(pf, "%f ", elements,0);

            if(c == x_ipl->height-1){
                fprintf(pf,"\n ",0);
            }
        }
    }
            //int col_size = freq->cols;
           // int row_size = freq->rows;
    fclose(pf);
#endif // 0
/*---------------------------------------------------------------*/


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
/* debug ------------------------------------------------------------*/
#if 0
    FILE * pf;
    pf = fopen("out_dst.txt", "w");
    for(r = 0; r<M->rows; r++){
        for(c = 0; c<M->cols; c++){
            float elements = cvmGet(M,r,c);
            //printf("%f ",elements);
            fprintf(pf, "%f ", elements,0);

            if(c == M->cols-1){
                fprintf(pf,"\n ",0);
            }
        }
    }
            int col_size = M->cols;
            int row_size = M->rows;
    fclose(pf);
#endif // 0
/*---------------------------------------------------------------*/
/* debug ------------------------------------------------------------*/
#if 0
    //find(M)
    FILE * pf;
    pf = fopen("out_dst.txt", "w");

    int d_count = 0;
    for(c = 0; c< M->cols; c++){
        for(r = 0; r< M->rows; r++){
            if(cvmGet(M,r,c) != 0){
                d_count=d_count+1;
            }
        }
    }

    int d_a= d_count;
    CvMat *d_mi = cvCreateMat(d_a,1,CV_32F);
    d_count= 0;
    for(c = 0; c< M->cols; c++){
        for(r = 0; r< M->rows; r++){
            if(cvmGet(M,r,c) != 0.0){
                float elements = c * M->rows + r+1;
                cvmSet(d_mi,d_count, 0, elements);    //find the nonzero elements in matrix
                d_count=d_count+1;
            }
        }
    }
    d_count = 0;

    for(c = 0; c<d_mi->cols; c++){
        for(r = 0; r<d_mi->rows; r++){
            float elements = cvmGet(d_mi,r,c);
            //printf("%f ",elements);
            fprintf(pf, "%f ", elements,0);

            if(c == d_mi->cols-1){
                fprintf(pf,"\n ",0);
            }
        }
    }
            int col_size = d_mi->cols;
            int row_size = d_mi->rows;

    fclose(pf);
#endif // 0
/*---------------------------------------------------------------*/

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

/* debug ------------------------------------------------------------*/
#if 0
    FILE * pf;
    pf = fopen("out_dst.txt", "w");
    for(r = 0; r<mhi->rows; r++){
        for(c = 0; c<mhi->cols; c++){
            float elements = cvmGet(mhi,r,c);
            //printf("%f ",elements);
            fprintf(pf, "%f ", elements,0);

            if(c == mhi->cols-1){
                fprintf(pf,"\n ",0);
            }
        }
    }
            int col_size = mhi->cols;
            int row_size = mhi->rows;
    fclose(pf);
#endif // 0
/*---------------------------------------------------------------*/

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

    //CvMat *b;
    //A_fhp(x, OMEGA, &b);

    CvMat *picks = cvCreateMat(OMEGA->rows, OMEGA->cols, CV_32F);
    cvmCopy(OMEGA,picks);

   /* matlab : F = I */
    CvMat *F_Mat = cvCreateMat(x_mat->rows, x_mat->cols, CV_32F);
    cvmCopy(x_mat,F_Mat);

    int m = F_Mat->rows;
    n = F_Mat->cols;

/* debug ------------------------------------------------------------*/
#if 0
    FILE * pf;
    pf = fopen("out_dst.txt", "w");
    for(r = 0; r<F_Mat->rows; r++){
        for(c = 0; c<F_Mat->cols; c++){
            float elements = cvmGet(F_Mat,r,c);
            //printf("%f ",elements);
            fprintf(pf, "%f ", elements,0);

            if(c == F_Mat->cols-1){
                fprintf(pf,"\n ",0);
            }
        }
    }
            int col_size = F_Mat->cols;
            int row_size = F_Mat->rows;
    fclose(pf);
#endif // 0
/*---------------------------------------------------------------*/

    /* matlab : FB = fft2(F)/sqrt(m*n); */
    CvMat *FB = cvCreateMat(F_Mat->rows, F_Mat->cols, CV_32F);

    CvMat *imaginaryInput = cvCreateMat(F_Mat->rows,F_Mat->cols,CV_32F);
    CvMat *imageinaryOutput = cvCreateMat(F_Mat->rows,F_Mat->cols,CV_32F);

    CvMat *realOutput = cvCreateMat(F_Mat->rows,F_Mat->cols,CV_32F);

    CvMat *F_Mat_blank = cvCreateMat(F_Mat->rows,F_Mat->cols,CV_32FC1);

    CvMat *F_Mat_FFT = cvCreateMat(F_Mat->rows,F_Mat->cols,CV_32FC2);
    CvMat *F_Mat_blank_out = cvCreateMat(F_Mat->rows,F_Mat->cols,CV_32FC2);


    cvmSetZero(imaginaryInput);
    cvMerge(F_Mat, imaginaryInput,NULL, NULL, F_Mat_FFT);

    cvDFT(F_Mat_FFT, F_Mat_blank_out, CV_DXT_FORWARD ,0);
    cvSplit(F_Mat_blank_out,F_Mat_blank,imageinaryOutput, NULL, NULL);

    cvmAdd(F_Mat_blank,imageinaryOutput,F_Mat_blank);

    cvReleaseMat(&realOutput);
    cvReleaseMat(&imaginaryInput);
    cvReleaseMat(&imageinaryOutput);
    cvReleaseMat(&F_Mat_FFT);
    cvReleaseMat(&F_Mat_blank_out);

    //cvDFT(F_Mat_blank, F_Mat_blank,CV_DXT_INVERSE_SCALE,0);
/*
        CvMat *fft2_U_inv_blank = cvCreateMat(F_Mat_blank->rows, F_Mat_blank->cols, CV_32FC2);
        CvMat *U_blank = cvCreateMat(F_Mat_blank->rows, F_Mat_blank->cols, CV_32FC2);


        cvMerge(F_Mat_blank, imageinaryOutput,NULL, NULL, fft2_U_inv_blank);
        cvDFT(fft2_U_inv_blank, U_blank, CV_DXT_INVERSE_SCALE, 0);
        cvSplit(U_blank,F_Mat_blank,imageinaryOutput, NULL, NULL);
    //cvReleaseMat(&imageinaryOutput);
*/

/* debug ------------------------------------------------------------*/
#if 0
    FILE * pf;
    pf = fopen("out_dst.txt", "w");
    for(r = 0; r<F_Mat_blank->rows; r++){
        for(c = 0; c<F_Mat_blank->cols; c++){
            float elements = cvmGet(F_Mat_blank,r,c);

            //printf("%f ",elements);
            fprintf(pf, "%f ", elements,0);

            if(c == F_Mat_blank->cols-1){
                fprintf(pf,"\n ",0);
            }
        }
    }
            int col_size = F_Mat_blank->cols;
            int row_size = F_Mat_blank->rows;
    fclose(pf);
#endif // 0
/*---------------------------------------------------------------*/


    for(c = 0; c<F_Mat->cols; c++){
        for(r=0; r<F_Mat->rows; r++){
            float elements= cvmGet(F_Mat_blank,r,c)/sqrt(m*n);
            cvmSet(FB,r,c,elements);
        }
    }
/* debug ------------------------------------------------------------*/
#if 0
    FILE * pf;
    pf = fopen("out_dst.txt", "w");
    for(r = 0; r<FB->rows; r++){
        for(c = 0; c<FB->cols; c++){
            float elements = cvmGet(FB,r,c);
            //printf("%f ",elements);
            fprintf(pf, "%f ", elements,0);

            if(c == FB->cols-1){
                fprintf(pf,"\n ",0);
            }
        }
    }
            int col_size = FB->cols;
            int row_size = FB->rows;
    fclose(pf);
#endif // 0
/*---------------------------------------------------------------*/

    /* matlab : B=FB(Picks); */
    CvMat *B = cvCreateMat(picks->rows, picks->cols, CV_32F);

    for(r = 0; r< picks->rows; r++){
        int pick_blank = cvmGet(picks,r,0) -1 ;
        int quotient = pick_blank/ FB->cols;
        int remainder = pick_blank % FB->cols;
        float elements = cvmGet(FB,remainder,quotient);
        cvmSet(B, r, 0, elements);
    }
/* debug ------------------------------------------------------------*/
#if 0
    FILE * pf;
    pf = fopen("out_dst.txt", "w");
    for(r = 0; r<B->rows; r++){
        for(c = 0; c<B->cols; c++){
            float elements = cvmGet(B,r,c);
            //printf("%f ",elements);
            fprintf(pf, "%f ", elements,0);

            if(c == B->cols-1){
                fprintf(pf,"\n ",0);
            }
        }
    }
            int col_size = B->cols;
            int row_size = B->rows;
    fclose(pf);
#endif // 0
/*---------------------------------------------------------------*/

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
#if 0
    /* matlab idwt */
    double d_x[x->cols * x->rows], d_wav[wav->rows*wav->cols];
    double W[x->cols * x->rows];

     /* matlab dwt */
    double d_x_dwt[x->cols * x->rows], d_wav_dwt[wav->rows * wav->cols];
    double WT[x->cols * x->rows];

    for(r = 0; r< x->rows * x->cols; r++){
         d_x[r] = cvmGet(x,r,0);
         d_x_dwt[r] = cvmGet(x,r,0);
    }
    for(c = 0; c<wav->cols; c++){
         d_wav[c] = cvmGet(wav,0,c);
         d_wav_dwt[c] = cvmGet(wav,0,c);
    }
    //y :input
    midwt(W,d_wav,d_x);
    //x : input
    mdwt(d_x_dwt,d_wav_dwt,WT);

    /* double -> cvmat */
    CvMat *W_mat = cvCreateMat(x->rows, x->cols, CV_32F);
    CvMat *WT_mat = cvCreateMat(x->rows, x->cols, CV_32F);

    //int kk =0, count_aa = 0;
    for(r = 0; r<(x->rows*x->cols); r++){
        cvmSet(W_mat,r,0,W[r]);
        cvmSet(WT_mat,r,0,WT[r]);
    }
#endif // 0
/* debug ------------------------------------------------------------*/
#if 0
    FILE * pf;
    pf = fopen("out_dst.txt", "w");
    for(r = 0; r<WT_mat->rows; r++){
        for(c = 0; c<WT_mat->cols; c++){
            float elements = cvmGet(WT_mat,r,c);
            //printf("%f ",elements);
            fprintf(pf, "%f ", elements,0);

            if(c == WT_mat->cols-1){
                fprintf(pf,"\n ",0);
            }
        }
    }
            int col_size = WT_mat->cols;
            int row_size = WT_mat->rows;
    fclose(pf);
#endif // 0
/*---------------------------------------------------------------*/


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

/* debug ------------------------------------------------------------*/
#if 0
    FILE * pf;
    pf = fopen("out_dst.txt", "w");
    for(r = 0; r<pick->rows; r++){
        for(c = 0; c<pick->cols; c++){
            float elements = cvmGet(pick,r,c);
            //printf("%f ",elements);
            fprintf(pf, "%f ", elements,0);

            if(c == pick->cols-1){
                fprintf(pf,"\n ",0);
            }
        }
    }
            int col_size = pick->cols;
            int row_size = pick->rows;
    fclose(pf);
#endif // 0
/*---------------------------------------------------------------*/
/* debug ------------------------------------------------------------*/
#if 0
    //find(M)
    int d_count = 0;
    for(c = 0; c< pick->cols; c++){
        for(r = 0; r< pick->rows; r++){
            if(cvmGet(pick,r,c) != 0){
                d_count=d_count+1;
            }
        }
    }

#endif // 0
/*---------------------------------------------------------------*/
   /* matlab : range(I(:)); */
    double min, max;
    cvMinMaxLoc(x_mat,&min,&max, 0,0,0);

    int range = max -min; //difference between the maximun and the minimum of a I


    RecPF_constraint(m,n,aTV, aL1,pick,B,2, opts_para, WT_mat, W_mat,range,x_mat,constraint,wav);

//-------------------------------------------------------------------------------------------------------

    cvReleaseMat(&picks);
    cvReleaseMat(&FB);
    cvReleaseMat(&wav);
    //cvReleaseMat(&W);
   // cvReleaseMat(&WT);
    cvReleaseImage(&FB);
    //cvReleaseImage(&F);

  }
