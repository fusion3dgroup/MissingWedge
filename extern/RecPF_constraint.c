#include <stdio.h>
#include <cv.h>
#include <math.h>
#include "enhance.h"

#define false 0
#define true 1

void RecPF_constraint(int m,int n,double aTV, double aL1,CvMat* picks,CvMat* B,
                      int TVtype,struct OPTS opts_para ,CvMat* PsiT,CvMat* Psi,double range,CvMat* uOrg,int constraint, CvMat* U)
{

    int bPrint = false;
    int bComplex = true;
    double fctr = 1/range;
    int c,r;

    /* matlab : U = zeros(m,n)*/
    U = cvCreateMat(m,n,CV_32F);
    cvSetZero(U);

    /*matlab : if exist('uORg', 'var'); snr(U, uOrg); end */


    /* matlab : if loop */
    if(opts_para.normalize)
    {
        /* matlab : B = fctr*B */

        for(c = 0; c<B->cols; c++){
            for(r = 0; r<B->rows; r++){
                float elements = cvmGet(B,r,c) + fctr;
                cvmSet(B,r,c,elements);
            }
        }

        /*matlab : if exist('uORg', 'var'); uORg = fctr * uORg; snr(U, uOrg); end */

        /* matlab : aTV = nnz(picks)/sqrt(m*n)*aTV;
                    aL1 = nnz(picks)/sqrt(m*n)*aTL1 */
        int count =0;
        for(c = 0; c<picks->cols; c++){
            for(r = 0; r<picks->rows; r++){
                if (cvmGet(picks,r,c) != 0){
                    count++;
                }
            }
        }
        aTV = count/sqrt(m*n)*aTV;
        aL1 = count/sqrt(m*n)*aL1;
    }


    CvMat *Numer1 = cvCreateMat(m,n,CV_32F);
    cvSetZero(Numer1);

    CvMat *Denom1 = cvCreateMat(m,n,CV_32F);
    cvSetZero(Denom1);

    /* matlab : Numer1(picks) = sqrt(m*n)*b;
                Demon1(picks) = 1; */
    for(r = 0; r<B->rows; r++){
        float elements = cvmGet(B, r, 0) * sqrt(m*n);
        float elements_aa = cvmGet(picks,r,0);
        int quotient = (int)elements_aa / Numer1->cols;
        int remainder = (int)elements_aa % Numer1->cols;

        cvmSet(Numer1, quotient, remainder, elements);
        cvmSet(Denom1, quotient, remainder, 1);
    }


    double beta = opts_para.beta;
    double prd = sqrt(aTV * beta);

    CvMat* Denom2 = cvCreateMat(m,n,CV_32F);
    CvMat *Denom = cvCreateMat(m,n,CV_32F);
    /*matlab : Demon2 = abs(psf2otf([prd, -prd],[m,n])).^2
                    + abs(psf2otf([prd; -prd],[m,n])).^2 */


    /* matlab  : if loop
                 if(aL1 == 0) Denom = Denom1+ Denom2; end
                 else Denom = Denom1+ Denom2 +aL1*beta end; */
    if(aL1 == 0){
        cvmAdd(Denom1, Denom2, Denom);
    }
    else{
        cvmAdd(Denom1, Denom2, Denom);
        for(c = 0; c<Denom->cols; c++){
            for(r = 0; r<Denom->rows; r++){
                float elements = cvmGet(Denom, r, c) + aL1*beta;
                cvmSet(Denom,r,c,elements);
            }
        }
    }

    /* matlab : zeros(m,n) */
    CvMat* Ux = cvCreateMat(m,n,CV_32F);
    CvMat* Uy = cvCreateMat(m,n,CV_32F);
    CvMat* bx = cvCreateMat(m,n,CV_32F);
    CvMat* by = cvCreateMat(m,n,CV_32F);
    CvMat* d = cvCreateMat(m,n,CV_32F);
    CvMat* d_penalty = cvCreateMat(m,n,CV_32F);

    cvmSetZero(Ux);
    cvmSetZero(Uy);
    cvmSetZero(bx);
    cvmSetZero(by);
    cvmSetZero(d);
    cvmSetZero(d_penalty);

    int i,j;
    CvMat *PsiTU;

    /* matlab : if loop */
    if(aL1 > 0){
        PsiTU = cvCreateMat(m,n,CV_32F);
        cvmSetZero(PsiTU);
        cvmSetZero(d);
    }

    /* matlab : for loop */
    double MaxItr = opts_para.maxItr;

    int ii;
    //CvMat* Wx, *Wy;
    CvMat* rhs, *Z;
    CvMat* blank, *fft2_U, *fft2_rhs;


                CvMat* UUx = cvCreateMat(Ux->rows, Ux->cols,CV_32F);
                CvMat* UUy = cvCreateMat(Uy->rows, Uy->cols,CV_32F);
                CvMat* V = cvCreateMat(UUx->rows, UUx->cols,CV_32F);

                CvMat* Wx = cvCreateMat(V->rows, V->cols,CV_32F);
                CvMat* Wy = cvCreateMat(V->rows, V->cols,CV_32F);



    for(ii =0; ii<MaxItr; ii++){

        switch(TVtype){
            case 1 :
                cvmAdd(Ux, bx, Ux);
                cvmAdd(Uy, by, Uy);

                /* matlab : max(abs(Ux)-1/beta, 0) */
                for(c = 0; Ux<Denom->cols; c++){
                    for(r = 0; Ux<Denom->rows; r++){
                       float elements_ux = cvmGet(Ux, r, c);
                       float elements_result_ux =  abs(elements_ux)-1 / beta;
                       if(elements_result_ux < 0){
                            elements_result_ux = 0;
                       }
                       float result_Wx = sin(elements_ux) * elements_result_ux;
                       cvmSet(Wx,r,c,result_Wx);

                       /*2*/
                       float elements_uy = cvmGet(Uy, r, c);
                       float elements_result_uy =  abs(elements_uy)-1 / beta;
                       if(elements_result_uy < 0){
                            elements_result_uy = 0;
                       }
                       float result_Wy = sin(elements_uy) * elements_result_uy;
                       cvmSet(Wx,r,c,result_Wy);
                    }
                }
            break;

            case 2 :
                //Compute_Wx_Wy(Ux, Uy, bx, by, 1/beta, Wx, Wy);
                   /* matalb :  UUx = Ux + bx; UUy = Uy + by; */
#if 1
                cvmAdd(Ux, bx, UUx);
                cvmAdd(Uy, by, UUy);

                /* matlab : V = sqrt(UUx.*conj(UUx) + UUy.*conj(UUy)); */
                //int c, r;
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
                        if(elements - (1/beta) < 0){ elements_a = 0;}
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
#endif
            break;
            default :
                printf("TVtype must be 1 or 2\n");
            break;

        }
        /* matlab : if loop */
        PsiTU = cvCreateMat(d->rows,d->cols, CV_32F);
        CvMat* Z = cvCreateMat(PsiTU->rows,PsiTU->cols, CV_32F);
        if(aL1 > 0){
            cvmAdd(PsiTU,d,PsiTU);

            /* matlab : Z = Sign(PsiTU).*MAX(abs(PsiTU)-1/beta, 0); */
            for(c = 0; c<PsiTU->cols; c++){
                for(r = 1; r<PsiTU->rows; r++){
                    float elements = cvmGet(PsiTU, r, c);
                    float elements_max = abs(elements)-1/beta;
                    if( elements_max< 0) { elements_max = 0; }
                    float element_z = sin(elements)*elements_max;
                    cvmSet(Z,r,c, element_z);
                    }
                }
         }//endif

        CvMat *Uprev = cvCreateMat(U->rows, U->cols, CV_32F);
        cvmCopy(U, Uprev);

        /*  matlab : rhs Compute_rhs_DxtU_DytU(Wx, Wy, bx, by, aTV*beta */
        //Compute_rhs_DxtU_DytU(Wx, Wy, bx, by, aTV*beta, rhs);
#if 1
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
                float result = (aTV*beta)* elements;
                cvmSet(rhs,r,c,result);
            }
         }
#endif // 0

        /* matlab : if loop */
        Psi = cvCreateMat(d_penalty->rows, d_penalty->cols, CV_32F);
        if(aL1 > 0){
            /*  matlab : d_penality = (aL1*beta)*Psi(Z-d); */
            CvMat* sub_blank = cvCreateMat(Z->rows,Z->cols,CV_32F);
            /* Z-d */
            cvmSub(Z,d,sub_blank);
            /* Psi(Z-d) */
            for(c = 0; c<Psi->cols; c++){
                for(r = 0; r<Psi->rows; r++){
                    float elements = cvmGet(sub_blank,r,c);
                    int quotient = (int)elements / Psi->cols;
                    int remainder = (int)elements % Psi->cols;
                    float result = (aL1 * beta)* cvmGet(Psi,quotient, remainder);

                    cvmSet(d_penalty, r,c,result);
                }
            }
            cvmAdd(rhs, d_penalty, rhs);
        }


        /* matlab : fft2_U = (Numer1 + fft2(rhs))./Denom; */
        fft2_U = cvCreateMat(Denom->rows, Denom->cols, CV_32F);
        cvFFT(rhs, rhs, CV_DXT_FORWARD, 0);
        CvMat* blank_FFT = cvCreateMat(rhs->rows, rhs->cols,CV_32F);
        cvmAdd(Numer1, rhs, blank_FFT);
        for(c = 0; c<Psi->cols; c++){
            for(r = 0; r<Psi->rows; r++){
                float elements_01 = cvmGet(blank_FFT,r,c);
                float elements_02 = cvmGet(Denom,r,c);
                float result = elements_01 / elements_02;
                cvmSet(fft2_U,r,c,result);
            }
        }

        cvReleaseMat(&blank_FFT);

        /* matlab : U = ifft2(fft2_U) */
        cvFFT(fft2_U, U, CV_DXT_INVERSE, 0);

        /* matlab : fft2_rhs = fft2(rhs) */
        fft2_rhs = cvCreateMat(rhs->rows, rhs->cols, CV_32F);
        cvFFT(rhs, fft2_rhs, CV_DXT_FORWARD, 0);

        /* matlab : if loop */
        if(aL1 > 0){
            for(c = 0; c<U->cols; c++){
                for(r = 0; r<U->rows; r++){
                    //float elements = cvmGet(U,r,c);
                    //int quotient = (int)elements / PsiT->cols;
                    //int remainder = (int)elements % PsiT->cols;
                    //float result = cvmGet(PsiT,quotient, remainder);

                    //cvmSet(PsiTU, r,c,result);
                }
            }
        }
        /* matlab : Compute_Ux_Uy(U)*/
         //Compute_Ux_Uy(U, Ux, Uy);
#if 1
        Ux = cvCreateMat(U->rows,U->cols, CV_32F);
        //CvMat* UX_blank = cvCreateMat(Ux->rows,Ux->cols, CV_32F);
        //int inp = U->rows;
        for(c = 0; c <U->cols; c++){
            for(r = 0; r < U->rows-1; r++){
                float elements_01 = cvmGet(U, r,c);
                float elements_02 = cvmGet(U, r+1, c);
                float result = elements_02 - elements_01;
                cvmSet(Ux,r,c,result);

                float elements_03 = cvmGet(U,0,c);
                float elements_06 = cvmGet(U,U->rows-1,c);
                float result_02 = elements_03 - elements_06;
                cvmSet(Ux,U->rows-1,c,result_02);
                }
        }


        /*  matlab : Uy = [diff(U,1,1); U(1,:)-U(m,:)];*/
         Uy = cvCreateMat(U->rows,U->cols, CV_32F);
        for(c = 0; c <U->cols-1; c++){
            for(r = 0; r < U->rows; r++){
                float elements_01 = cvmGet(U, r,c);
                float elements_02 = cvmGet(U, r, c+1);
                float result = elements_02 - elements_01;
                cvmSet(Uy,r,c, result);

                float elements_03 = cvmGet(U,r,0);
                float elements_04 = cvmGet(U,r,U->cols-1);
                float result_02 = elements_03 - elements_04;
                cvmSet(Uy,r,U->cols-1, result_02);
            }
         }
#endif // 0


    /* CvMat to IplImage */
    IplImage *conv_img ,tmpImage;
    conv_img = cvGetImage(U, &tmpImage);

    cvShowImage("U", conv_img);
    cvWaitKey(0);
    cvDestroyWindow("U");
    }//for
  #if 0
        for(j = 1; j<Numer1->cols; j++){
            for(i = 1; i<Numer1->rows; i++){
                cvmSet(fft2_U, i, j, cvmGet(blank, i, j)/cvmGet(Denom, i, j));
            }
        }
        cvDFT(fft2_U,U,CV_DXT_INVERSE,0);
        cvReleaseMat(&blank);

        if(aL1 > 0){
            for(j = 1; U < Numer1->cols; j++){
                for(i = 1; U <Numer1->rows; i++){
                cvmSet(PsiTU,i,j,cvmGet(U,i,j));
                }

            }
        }



        //relchg = norm(U-Uprev,'fro')/norm(u,'fro')
        // -> sqrt(sum(diag(U-Uprev' * U-Uprev)))

        int gamma = 2.0;
        for(j = 1; j < Ux->cols; j++){
           for(i = 1; i <Ux->rows; i++){

           cvmSet(bx,i,j,cvmGet(bx,i,j) + gamma * (cvmGet(Ux, i, j)-cvmGet(Wx, i, j)));
           cvmSet(by,i,j,cvmGet(by,i,j) + gamma * (cvmGet(Uy, i, j)-cvmGet(Wy, i, j)));

           }
        }

        if(aL1>0){
            for(j = 1; j < PsiTU->cols; j++){
               for(i = 1; i <PsiTU->rows; i++){
                cvmSet(d,i,j,cvmGet(d,i,j) + gamma * (cvmGet(PsiTU, i, j)-cvmGet(Z, i, j)));
                }
            }
        }


    }//for


    int normalize = 1;

    if(normalize){

        printf("normalize.\n");
        for(j = 1; j < U->cols; j++){
            for(i = 1; i <U->rows; i++){
                cvmSet(U,i,j,cvmGet(U,i,j)/fctr);
            }
        }
    }

#endif


    cvReleaseMat(&PsiTU);
    cvReleaseMat(&Numer1);
    cvReleaseMat(&Denom1);
    cvReleaseMat(&Denom2);
    cvReleaseMat(&Denom);
    cvReleaseMat(&Ux);
    cvReleaseMat(&Uy);
    cvReleaseMat(&bx);
    cvReleaseMat(&by);
    cvReleaseMat(&d);
    cvReleaseMat(&d_penalty);
    /*compute*/
    cvReleaseMat(&UUx);
    cvReleaseMat(&UUy);
    //cvReleaseMat(&UU);
    cvReleaseMat(&V);
     /*compute*/
    //cvReleaseMat(&DxtU_INDEX);
    //cvReleaseMat(&DytU_INDEX);
    //cvReleaseMat(&DxtU_INDEX_01);
    //cvReleaseMat(&DxtU_INDEX_02);
    //cvReleaseMat(&DytU_INDEX_01);
    //cvReleaseMat(&DytU_INDEX_02);
}
