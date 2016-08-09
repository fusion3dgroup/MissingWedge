#include <stdio.h>
#include <cv.h>
#include <math.h>
#include <complex.h>

#include "enhance.h"

#define false 0
#define true 1
void snr(CvMat *sig, CvMat *ref){

    int c,r;

    CvMat* ref_save = cvCreateMat(ref->rows, ref->cols, CV_32F);
    int n = ref_save->cols*ref_save->rows;

    cvmCopy(ref, ref_save);

    float mean;
    float ref_aver=0.;
    for(c = 0; c<ref_save->cols; c++){
        for(r = 0; r<ref_save->rows; r++){
            float ref_aver = cvmGet(ref_save,r,c);
            float elements_B = cvmGet(sig,r,c);

           mean = (ref_aver - elements_B) * (ref_aver - elements_B);
           mean = mean+mean;

            ref_aver = ref_aver + ref_aver;

        }
    }
    float mse = mean/n;

    float ref_avr = ref_aver / n;

    float elements_sum;
    for(c = 0; c<ref_save->cols; c++){
        for(r = 0; r<ref_save->rows; r++){

            float elements = cvmGet(ref_save,r,c);
            elements_sum = (elements - ref_avr) * (elements - ref_avr);
            elements_sum = elements_sum+elements_sum;
        }
    }
    float dv = elements_sum /(n-1);

    float x = 10*log10(dv/mse);

}
void RecPF_constraint(int m,int n,double aTV, double aL1,CvMat* picks,CvMat* B,
                      int TVtype,struct OPTS opts_para ,CvMat* PsiT,CvMat* Psi,int range,CvMat* uOrg,int constraint, CvMat *wav)
{

    int bPrint = false;
    int bComplex = true;
    double fctr;

    int c,r, p,o;
    double gamma = opts_para.gamma;

    /* matlab : U = zeros(m,n)*/
    CvMat* U = cvCreateMat(m,n,CV_32F);
    cvSetZero(U);

    /*matlab : if exist('uORg', 'var'); snr(U, uOrg); end */
    snr(U,uOrg);


    /* matlab : if loop */
    if(opts_para.normalize)
    {
        fctr = 1.0/range;
        /* matlab : B = fctr*B */
        for(c = 0; c<B->cols; c++){
            for(r = 0; r<B->rows; r++){
                float elements = cvmGet(B,r,c) * fctr;
                cvmSet(B,r,c,elements);
            }
        }

        /*matlab : if exist('uORg', 'var'); uORg = fctr * uORg; snr(U, uOrg); end */
        for(c = 0; c<uOrg->cols; c++){
            for(r = 0; r<uOrg->rows; r++){
                float elements = cvmGet(uOrg,r,c) * fctr;
                cvmSet(uOrg,r,c,elements);
            }
        }
        snr(U, uOrg);

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
    }//if
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

    /*Numer1 = zeros(m,m); Denom1 = zeros(m,n);*/
    CvMat *Numer1 = cvCreateMat(m,n,CV_32F);
    cvSetZero(Numer1);

    CvMat *Denom1 = cvCreateMat(m,n,CV_32F);
    cvSetZero(Denom1);

    /* matlab : Numer1(picks) = sqrt(m*n)*B;
                Demon1(picks) = 1; */
    //for(r = 0; r<B->rows; r++){
    //    int elements = cvmGet(B, r, 0) * sqrt(m*n);
    //    //float elements_aa = cvmGet(picks,r,0);
    //    int quotient = elements_aa / Numer1->rows;
    //    int remainder = (int)elements_aa % Numer1->rows;

    //    cvmSet(Numer1, quotient, remainder, elements);
    //    //cvmSet(Denom1, quotient, remainder, 1);
    //}
    int count = 0;
    for(c = 0; c<picks->cols; c++){
        for(r = 0; r<picks->rows; r++){

            if (cvmGet(picks, r, c) !=0){
                float elements = cvmGet(B,count,0) * sqrt(m*n);
                cvmSet(Numer1,r,c,elements);
                cvmSet(Denom1,r,c,1.0);
                count = count + 1;
            }

        }
    }
    count = 0;

/* debug ------------------------------------------------------------*/
#if 0
    FILE * pf;
    pf = fopen("out_dst.txt", "w");
    for(r = 0; r<Denom1->rows; r++){
        for(c = 0; c<Denom1->cols; c++){
            float elements = cvmGet(Denom1,r,c);
            //printf("%f ",elements);
            fprintf(pf, "%f ", elements,0);

            if(c == Denom1->cols-1){
                fprintf(pf,"\n ",0);
            }
        }
    }
            int col_size = Denom1->cols;
            int row_size = Denom1->rows;
    fclose(pf);
#endif // 0
/*---------------------------------------------------------------*/

    double beta = opts_para.beta;
    double prd = sqrt(aTV * beta);

    CvMat* Denom2 = cvCreateMat(m,n,CV_32F);
    CvMat *Denom = cvCreateMat(m,n,CV_32F);

    /*matlab : Demon2 = abs(psf2otf([prd, -prd],[m,n])).^2
                    + abs(psf2otf([prd; -prd],[m,n])).^2 */

    /* circshift([prd, -prd],-foor(size([prd, -prd]/2)))
                    ots_A
        ffte(ots_A) */


     CvMat *ots_A = cvCreateMat(m,n,CV_32F);
     cvmSetZero(ots_A);
     cvmSet(ots_A,0,0,-prd);
     cvmSet(ots_A,0,1, prd);




    CvMat *ots_A_FFT = cvCreateMat(ots_A->rows, ots_A->cols, CV_32FC2);
    CvMat *imag_A = cvCreateMat(ots_A->rows, ots_A->cols, CV_32F);

    cvmSetZero(imag_A);
    cvMerge(ots_A, imag_A,NULL, NULL, ots_A_FFT);

    cvDFT(ots_A_FFT, ots_A_FFT, CV_DXT_FORWARD ,0);
    cvSplit(ots_A_FFT,ots_A,imag_A, NULL, NULL);
    cvmAdd(ots_A,imag_A,ots_A);
     //fabs(ots_A) * fabs(ots_A);

    CvMat *ots_B = cvCreateMat(m,n,CV_32F);
     cvmSetZero(ots_B);
    cvmSet(ots_B,0,0,-prd);
    cvmSet(ots_B,1,0, prd);

    CvMat *ots_B_FFT = cvCreateMat(ots_B->rows, ots_B->cols, CV_32FC2);
    CvMat *imag_B = cvCreateMat(ots_B->rows, ots_B->cols, CV_32F);
    cvmSetZero(imag_B);
    cvMerge(ots_B, imag_B,NULL, NULL, ots_B_FFT);

    cvDFT(ots_B_FFT, ots_B_FFT, CV_DXT_FORWARD ,0);
    cvSplit(ots_B_FFT,ots_B,imag_B, NULL, NULL);
    cvmAdd(ots_B,imag_B,ots_B);

    for(c = 0; c<Denom2->cols; c++){
       for(r = 0; r<Denom2->rows; r++){

                        //double complex ots_a = conj(cvmGet(ots_A, r, c));
                        //double complex ots_b = conj(cvmGet(ots_B, r, c));

            float elements_A = fabs(cvmGet(ots_A, r,c));
            float elements_B = fabs(cvmGet(ots_B, r,c));

            //double result = fabs(elements_A);
            float result =  elements_A*elements_A+ elements_B*elements_B;
            cvmSet(Denom2, r,c,result);
        }
    }

 /* debug ------------------------------------------------------------*/
#if 0
    FILE * pf;
    pf = fopen("out_dst.txt", "w");
    for(r = 0; r<Denom2->rows; r++){
        for(c = 0; c<Denom2->cols; c++){
            float elements = cvmGet(Denom2,r,c);
            //printf("%f ",elements);
            fprintf(pf, "%f ", elements,0);

            if(c == Denom2->cols-1){
                fprintf(pf,"\n ",0);
            }
        }
    }
            int col_size = Denom2->cols;
            int row_size = Denom2->rows;
    fclose(pf);
#endif // 0
/*---------------------------------------------------------------*/

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

    //int i,j;
    /**/
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



    for(ii =0; ii<MaxItr; ii++){

    CvMat *fft2_U = cvCreateMat(Denom->rows, Denom->cols, CV_32F);

    CvMat* UUx = cvCreateMat(Ux->rows, Ux->cols,CV_32F);
    CvMat* UUy = cvCreateMat(Uy->rows, Uy->cols,CV_32F);
    CvMat* V = cvCreateMat(UUx->rows, UUx->cols,CV_32F);

    CvMat* Wx = cvCreateMat(V->rows, V->cols,CV_32F);
    CvMat* Wy = cvCreateMat(V->rows, V->cols,CV_32F);
    CvMat* DxtU_INDEX = cvCreateMat(Wx->rows, Wx->cols, CV_32F);
    CvMat* DytU_INDEX = cvCreateMat(Wy->rows, Wy->cols, CV_32F);

    CvMat *DxtU_INDEX_01 = cvCreateMat(DxtU_INDEX->rows,1,CV_32F);
    CvMat* DxtU_INDEX_02 = cvCreateMat(DxtU_INDEX->rows,DxtU_INDEX->cols-1, CV_32F);
    CvMat* DxtU = cvCreateMat(DxtU_INDEX->rows,DxtU_INDEX->cols, CV_32F);
    CvMat *DytU_INDEX_01 = cvCreateMat(1,DytU_INDEX->cols,CV_32F);
    CvMat* DytU_INDEX_02 = cvCreateMat(DytU_INDEX->rows-1,DytU_INDEX->cols, CV_32F);
    CvMat* DytU = cvCreateMat(DytU_INDEX->rows,DytU_INDEX->cols, CV_32F);

    CvMat *rhs = cvCreateMat(DxtU->rows, DxtU->cols, CV_32F);
    CvMat *fft2_rhs = cvCreateMat(rhs->rows, rhs->cols, CV_32F);
    CvMat* Z = cvCreateMat(PsiTU->rows,PsiTU->cols, CV_32F);
    CvMat *Uprev = cvCreateMat(U->rows, U->cols, CV_32F);

    PsiTU = cvCreateMat(d->rows,d->cols, CV_32F);
    Psi = cvCreateMat(d_penalty->rows, d_penalty->cols, CV_32F);

    Ux = cvCreateMat(U->rows,U->cols, CV_32F);
    Uy = cvCreateMat(U->rows,U->cols, CV_32F);

        switch(TVtype){
            case 1 :
                cvmAdd(Ux, bx, Ux);
                cvmAdd(Uy, by, Uy);

                /* matlab : max(abs(Ux)-1/beta, 0) */
                for(c = 0; Ux<Denom->cols; c++){
                    for(r = 0; Ux<Denom->rows; r++){
                        int flag_elements = 0;
                       float elements_ux = cvmGet(Ux, r, c);
                       float elements_result_ux =  abs(elements_ux)-1 / beta;
                       if(elements_result_ux < 0){
                            elements_result_ux = 0;
                       }
                        if(elements_ux < 0){
                        flag_elements = -1;
                        }
                        else if(elements_ux == 0){
                        flag_elements = 0;
                        }
                        else{
                        flag_elements = 1;
                        }
                       float result_Wx = flag_elements * elements_result_ux;
                       cvmSet(Wx,r,c,result_Wx);

                       /*2*/
                       float elements_uy = cvmGet(Uy, r, c);
                       float elements_result_uy =  abs(elements_uy)-1 / beta;
                       if(elements_result_uy < 0){
                            elements_result_uy = 0;
                       }

                        if(elements_uy < 0){
                        flag_elements = -1;
                        }
                        else if(elements_uy == 0){
                        flag_elements = 0;
                        }
                        else{
                        flag_elements = 1;
                        }
                       float result_Wy = flag_elements * elements_result_uy;
                       cvmSet(Wx,r,c,result_Wy);
                    }
                }
            break;
//------------------------------------------------------------------------------------------------------
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

                        double complex uux_c = conj(cvmGet(UUx, r, c));
                        double complex uuy_c = conj(cvmGet(UUy, r, c));

                        float elements_uux = cvmGet(UUx, r, c);
                        float elements_uuy = cvmGet(UUy, r, c);
                        float result = sqrt((elements_uux * creal(uux_c) - elements_uux * cimag(uux_c) )+ (elements_uuy * creal(uuy_c) - elements_uuy * cimag(uuy_c)));
                        //float result = sqrt((elements_uux * (creal(uux_c) + cimag(uux_c)))+ (elements_uuy * (creal(uuy_c)+ cimag(uuy_c))));
                        cvmSet(V,r,c,result);
                    }
                }
                /* matlab : V = max(V - tau, 0) ./ max(V,eps); */
                float eps = 2.2204e-16;
                for(c = 0; c<V->cols; c++){
                    for(r = 0; r<V->rows; r++){
                        float elements = cvmGet(V,r,c);
                        float elements_a = 0, elements_b = 0;
                        if((elements - (1.0/beta)) <= 0.0){ elements_a = 0.0;}
                        else{elements_a = elements;}
                        if(elements <= eps) { elements_b = eps;}
                        else{elements_b = elements;}
                        float result = elements_a / elements_b;
                        cvmSet(V,r,c, result);
                    }
                }
                /* matlab  :   Wx = V.*UUx; Wy = V.*UUy;*/

                for(c = 0; c<V->cols; c++){
                    for(r = 0; r<V->rows; r++){
                        float elements = cvmGet(V,r,c);
                        float elements_a = cvmGet(UUx,r,c);
                        float result = elements * elements_a;
                        cvmSet(Wx,r,c,result);

                        float elements_b = cvmGet(UUy,r,c);
                        float result_b = elements * elements_b;
                        cvmSet(Wy,r,c,result_b);
                    }
                }
                //cvmMul(V, UUx, Wx);
                //cvmMul(V, UUy, Wy);
#endif
            break;
            default :
                printf("TVtype must be 1 or 2\n");
            break;

        } //switch
 /* debug ------------------------------------------------------------*/
#if 0
    FILE * pf;
    pf = fopen("out_dst.txt", "w");
    for(r = 0; r<Wy->rows; r++){
        for(c = 0; c<Wy->cols; c++){
            float elements = cvmGet(Wy,r,c);
            //printf("%f ",elements);
            fprintf(pf, "%f ", elements,0);

            if(c == Wy->cols-1){
                fprintf(pf,"\n ",0);
            }
        }
    }
            int col_size = Wy->cols;
            int row_size = Wy->rows;
    fclose(pf);
#endif // 0
/*---------------------------------------------------------------*/
        /* matlab : if loop */
       if(aL1 > 0){
            cvmAdd(PsiTU,d,PsiTU);

            /* matlab : Z = Sign(PsiTU).*MAX(abs(PsiTU)-1/beta, 0); */
            for(c = 0; c<PsiTU->cols; c++){
                for(r = 1; r<PsiTU->rows; r++){
                    float elements = cvmGet(PsiTU, r, c);
                    float elements_max = fabs(elements)-(1/beta);
                    //float elements_a =
                    if( elements_max<= 0) { elements_max = 0; }
                    else{elements_max = elements_max; }

                    float flag_elements = 0.0;
                    //sign
                    if(elements < 0){ flag_elements = -1.0;}
                    else if(elements == 0){ flag_elements = 0.0;}
                    else{flag_elements = 1.0;}

                    float element_z = flag_elements*elements_max;
                    cvmSet(Z,r,c, element_z);
                 }
             }
         }//endif
 /* debug ------------------------------------------------------------*/
#if 0
    FILE * pf;
    pf = fopen("out_dst.txt", "w");
    for(r = 0; r<PsiTU->rows; r++){
        for(c = 0; c<PsiTU->cols; c++){
            float elements = cvmGet(PsiTU,r,c);
            //printf("%f ",elements);
            fprintf(pf, "%f ", elements,0);

            if(c == PsiTU->cols-1){
                fprintf(pf,"\n ",0);
            }
        }
    }
            int col_size = PsiTU->cols;
            int row_size = PsiTU->rows;
    fclose(pf);
#endif // 0
/*---------------------------------------------------------------*/
        cvmCopy(U, Uprev);

        /*  matlab : rhs Compute_rhs_DxtU_DytU(Wx, Wy, bx, by, aTV*beta */
        //Compute_rhs_DxtU_DytU(Wx, Wy, bx, by, aTV*beta, rhs);
#if 1
        /* matlab : RHS = tau*(DxtU(Wx-bx)+DytU(Wy-by));
                                  DxtU         DyTU        */

        cvmSub(Wx, bx, DxtU_INDEX);
        cvmSub(Wy, by, DytU_INDEX);

        double complex tau = conj(aTV*beta);

         /* matlab  : dxtu = [U(:,end)-U(:, 1) U(:,1:end-1)-U(:,2:end)];
                                DxtU_INDEX_01         DxtU_INDEX_02     */

        for(r = 0; r<DxtU_INDEX->rows; r++){
            /* U(:,end)*/
            float elements_01 = cvmGet(DxtU_INDEX,r,DxtU_INDEX->cols-1);
            /* U(:, 1) */
            float elements_02 = cvmGet(DxtU_INDEX,r,0);
            float result = elements_01 - elements_02;
            cvmSet(DxtU_INDEX_01, r,0,result);
        }

        for(c = 0; c<DxtU_INDEX->cols-1; c++){
            for(r = 0; r<DxtU_INDEX->rows ; r++){
             /* U(:,1:end-1) */
            float elements_01 = cvmGet(DxtU_INDEX,r,c);
             /* U(:,2:end) */
            float elements_02 = cvmGet(DxtU_INDEX, r,c+1);
            float result = elements_01 - elements_02;
            cvmSet(DxtU_INDEX_02, r,c,result);

            }
        }

        for(c = 0; c<DxtU_INDEX->cols-1; c++){
            for(r = 0; r<DxtU_INDEX->rows  ; r++){
                if(c == 0){
                    float elements = cvmGet(DxtU_INDEX_01, r, 0);
                    cvmSet(DxtU, r,0,elements);
                }
                else{
                    float elements_02 = cvmGet(DxtU_INDEX_02, r, c-1);
                    cvmSet(DxtU, r,c,elements_02);
                }
            }
        }

        /* matlab :  dytu = [U(end,:)-U(1, :); U(1:end-1,:)-U(2:end,:)];
                                DytU_INDEX_01       DytU_INDEX_02      */


        for(c = 0; c<DytU_INDEX->cols; c++){
            /* U(end,:)*/
            float elements_01 = cvmGet(DytU_INDEX,DytU_INDEX->rows-1,c);
            /* U(1, :) */
            float elements_02 = cvmGet(DytU_INDEX,0,c);
            float result = elements_01 - elements_02;
            cvmSet(DytU_INDEX_01, 0,c,result);
        }

        for(c = 0; c<DxtU_INDEX->cols; c++){
            for(r = 0; r<DytU_INDEX->rows-1 ; r++){
             /* U(1:end-1,:) */
            float elements_01 = cvmGet(DytU_INDEX,r,c);
             /* U(2:end,:) */
            float elements_02 = cvmGet(DytU_INDEX, r+1,c);
            float result = elements_01 - elements_02;
            cvmSet(DytU_INDEX_02, r,c,result);

            }
        }

        for(c = 0; c<DytU_INDEX->cols; c++){
            for(r = 0; r<DytU_INDEX->rows-1; r++){
                if(r == 0){
                    float elements = cvmGet(DytU_INDEX_01, 0, c);
                    cvmSet(DytU, 0,c,elements);
                }
                else{
                    float elements_02 = cvmGet(DytU_INDEX_02, r-1, c);
                    cvmSet(DytU, r,c,elements_02);
                }
            }
        }

        cvmAdd(DxtU,DytU,rhs);

        for(c = 0; c<rhs->cols-1; c++){
            for(r = 0; r<rhs->rows ; r++){
                float elements = cvmGet(rhs,r,c);
                float result = creal(tau)* elements;
                cvmSet(rhs,r,c,result);
            }
         }
#endif // 0
 /* debug ------------------------------------------------------------*/
#if 0
    FILE * pf;
    pf = fopen("out_dst.txt", "w");
    for(r = 0; r<rhs->rows; r++){
        for(c = 0; c<rhs->cols; c++){
            float elements = cvmGet(rhs,r,c);
            //printf("%f ",elements);
            fprintf(pf, "%f ", elements,0);

            if(c == rhs->cols-1){
                fprintf(pf,"\n ",0);
            }
        }
    }
            int col_size = rhs->cols;
            int row_size = rhs->rows;
    fclose(pf);
#endif // 0
/*---------------------------------------------------------------*/
        /* matlab : if loop */
        CvMat* sub_blank = cvCreateMat(Z->rows,Z->cols,CV_32F);
        double d_x[Z->cols*Z->rows], d_wav[wav->rows*wav->cols];
        double W[Z->cols*Z->rows];
        if(aL1 > 0){
            /*  matlab : d_penality = (aL1*beta)*Psi(Z-d); */

            /* Z-d */
            cvmSub(Z,d,Z);

            /* Psi(Z-d) */
            for(c = 0; c< Z->cols; c++){
                for(r = 0; r< Z->rows; r++){
                 d_x[r+c*256] = cvmGet(Z,r,c);
                }
            }
            for(c = 0; c<wav->cols; c++){
                 d_wav[c] = cvmGet(wav,0,c);
            }

            //y :input
            midwt(W,d_wav,d_x);

            for(c = 0; c< Z->cols; c++){
                for(r = 0; r< Z->rows; r++){
                    cvmSet(sub_blank,r,c,W[r+c*256]);
                }
             }
            /* (aL1*beta)*Psi(Z-d)*/
            for(c = 0; c<d_penalty->cols; c++){
                for(r = 0; r<d_penalty->rows; r++){
                    float result = (aL1 * beta)* cvmGet(sub_blank,r, c);
                    cvmSet(d_penalty, r,c,result);
                }
            }
            cvmAdd(rhs, d_penalty, rhs);
        }
 /* debug ------------------------------------------------------------*/
#if 0
    FILE * pf;
    pf = fopen("out_dst.txt", "w");
    for(r = 0; r<Z->rows; r++){
        for(c = 0; c<Z->cols; c++){
            float elements = cvmGet(Z,r,c);
            //printf("%f ",elements);
            fprintf(pf, "%f ", elements,0);

            if(c == sub_blank->cols-1){
                fprintf(pf,"\n ",0);
            }
        }
    }
            int col_size = sub_blank->cols;
            int row_size = sub_blank->rows;
    fclose(pf);
#endif // 0
/*---------------------------------------------------------------*/

        /* matlab : fft2_U = (Numer1 + fft2(rhs))./Denom;
                                       F_Mat_blank_Rec               */

    CvMat *imaginaryInput_Rec = cvCreateMat(rhs->rows,rhs->cols,CV_32F);
    CvMat *imageinaryOutput_Rec = cvCreateMat(rhs->rows,rhs->cols,CV_32F);
    CvMat *realOutput_Rec = cvCreateMat(rhs->rows,rhs->cols,CV_32F);
    CvMat *F_Mat_blank_Rec = cvCreateMat(rhs->rows,rhs->cols,CV_32F);
    CvMat *F_Mat_FFT_Rec = cvCreateMat(rhs->rows,rhs->cols,CV_32FC2);
    CvMat *F_Mat_blank_out_Rec = cvCreateMat(rhs->rows,rhs->cols,CV_32FC2);


    cvmSetZero(imaginaryInput_Rec);
    cvMerge(rhs, imaginaryInput_Rec,NULL, NULL, F_Mat_FFT_Rec);

    cvDFT(F_Mat_FFT_Rec, F_Mat_blank_out_Rec, CV_DXT_FORWARD ,0);
    cvSplit(F_Mat_blank_out_Rec,F_Mat_blank_Rec,imageinaryOutput_Rec, NULL, NULL);
    cvmAdd(F_Mat_blank_Rec,imageinaryOutput_Rec,F_Mat_blank_Rec);

    cvReleaseMat(&realOutput_Rec);
    cvReleaseMat(&imaginaryInput_Rec);

    cvReleaseMat(&F_Mat_FFT_Rec);
    cvReleaseMat(&F_Mat_blank_out_Rec);

 /* debug ------------------------------------------------------------*/
#if 0
    FILE * pf;
    pf = fopen("out_dst.txt", "w");
    for(r = 0; r<F_Mat_blank_Rec->rows; r++){
        for(c = 0; c<F_Mat_blank_Rec->cols; c++){
            float elements = cvmGet(F_Mat_blank_Rec,r,c);
            //printf("%f ",elements);
            fprintf(pf, "%f ", elements,0);

            if(c == F_Mat_blank_Rec->cols-1){
                fprintf(pf,"\n ",0);
            }
        }
    }
            int col_size = F_Mat_blank_Rec->cols;
            int row_size = F_Mat_blank_Rec->rows;
    fclose(pf);
#endif // 0
/*---------------------------------------------------------------*/


        CvMat* blank_FFT = cvCreateMat(F_Mat_blank_Rec->rows, F_Mat_blank_Rec->cols,CV_32F);

        cvmAdd(Numer1, F_Mat_blank_Rec, blank_FFT);
        for(c = 0; c<fft2_U->cols; c++){
            for(r = 0; r<fft2_U->rows; r++){
                float elements_01 = cvmGet(blank_FFT,r,c);
                float elements_02 = cvmGet(Denom,r,c);
                float result = elements_01 / elements_02;
                cvmSet(fft2_U,r,c,result);
            }
        }

        //cvReleaseMat(&blank_FFT);

        //cvDFT(fft2_U, U, CV_DXT_INVERSE_SCALE, 0);

        /* matlab : U = ifft2(fft2_U) */
        CvMat *fft2_U_inv_blank = cvCreateMat(fft2_U->rows, fft2_U->cols, CV_32FC2);
        CvMat *U_blank = cvCreateMat(fft2_U->rows, Denom->cols, CV_32FC2);


        cvMerge(fft2_U, imageinaryOutput_Rec,NULL, NULL, fft2_U_inv_blank);
        cvDFT(fft2_U_inv_blank, U_blank, CV_DXT_INVERSE_SCALE, 0);
        cvSplit(U_blank,U,imageinaryOutput_Rec, NULL, NULL);


        cvReleaseMat(&imageinaryOutput_Rec);

 /* debug ------------------------------------------------------------*/
#if 0
    FILE * pf;
    pf = fopen("out_dst.txt", "w");
    for(r = 0; r<Numer1->rows; r++){
        for(c = 0; c<Numer1->cols; c++){
            float elements = cvmGet(Numer1,r,c);
            //printf("%f ",elements);
            fprintf(pf, "%f ", elements,0);

            if(c == Numer1->cols-1){
                fprintf(pf,"\n ",0);
            }
        }
    }
            int col_size = Numer1->cols;
            int row_size = Numer1->rows;
    fclose(pf);
#endif // 0
/*---------------------------------------------------------------*/
        /* matlab : fft2_rhs = fft2(rhs) */
        cvmCopy(F_Mat_blank_Rec, fft2_rhs);
        /* if ~bComplex; U=real(u); end*/


        //ERROR
        /* matlab : if loop */
        /* matlab : PsiTU = PsiT(U); */

        /* PsiT(U) */
        double d_x_dwt[U->cols * U->rows], d_wav_dwt[wav->rows * wav->cols];
        double WT[U->cols * U->rows];

        if(aL1 > 0){
            for(c = 0; c< U->cols; c++){
                for(r = 0; r< U->rows; r++){
                 d_x_dwt[r+c*256] = cvmGet(U,r,c);
                }
            }
            for(c = 0; c<wav->cols; c++){
                 d_wav_dwt[c] = cvmGet(wav,0,c);
            }

            //x : input
            mdwt(d_x_dwt,d_wav_dwt,WT);


            for(c = 0; c< U->cols; c++){
                for(r = 0; r< U->rows; r++){
                    cvmSet(PsiTU,r,c,WT[r+c*256]);
                }
             }
         }

 /* debug ------------------------------------------------------------*/
#if 0
    FILE * pf;
    pf = fopen("out_dst.txt", "w");
    for(r = 0; r<PsiTU->rows; r++){
        for(c = 0; c<PsiTU->cols; c++){
            float elements = cvmGet(PsiTU,r,c);
            //printf("%f ",elements);
            fprintf(pf, "%f ", elements,0);

            if(c == PsiTU->cols-1){
                fprintf(pf,"\n ",0);
            }
        }
    }
            int col_size = PsiTU->cols;
            int row_size = PsiTU->rows;
    fclose(pf);
#endif // 0
/*---------------------------------------------------------------*/


        /* matlab : Compute_Ux_Uy(U)*/
         //Compute_Ux_Uy(U, Ux, Uy);
#if 1

        //CvMat* UX_blank = cvCreateMat(Ux->rows,Ux->cols, CV_32F);
        //int inp = U->rows;
         /*  matlab : Ux = [diff(U,1,2), U(:,1)-U(:,n)];*/
        for(c = 0; c <U->cols-1; c++){
            for(r = 0; r < U->rows; r++){
                float elements_01 = cvmGet(U, r,c);
                float elements_02 = cvmGet(U, r, c+1);
                float result = elements_02 - elements_01;
                cvmSet(Ux,r,c,result);

                float elements_03 = cvmGet(U,r,0);
                float elements_06 = cvmGet(U,r,U->cols-1);
                float result_02 = elements_03 - elements_06;
                cvmSet(Ux,r,U->cols-1,result_02);
                }
        }


        /*  matlab : Uy = [diff(U,1,1); U(1,:)-U(m,:)];*/
        for(c = 0; c <U->cols; c++){
            for(r = 0; r < U->rows-1; r++){
                float elements_01 = cvmGet(U, r,c);
                float elements_02 = cvmGet(U, r+1, c);
                float result = elements_02 - elements_01;
                cvmSet(Uy,r,c, result);

                float elements_03 = cvmGet(U,0,c);
                float elements_04 = cvmGet(U, U->rows-1,c);
                float result_02 = elements_03 - elements_04;
                cvmSet(Uy,U->rows-1,c, result_02);
            }
         }
#endif // 0
/* debug ------------------------------------------------------------*/

#if 0
    FILE * pf;
    pf = fopen("out_dst.txt", "w");
    for(r = 0; r<Uy->rows; r++){
        for(c = 0; c<Uy->cols; c++){
            float elements = cvmGet(Uy,r,c);
            //printf("%f ",elements);
            fprintf(pf, "%f ", elements,0);

            if(c == Uy->cols-1){
                fprintf(pf,"\n ",0);
            }
        }
    }
            int col_size = Uy->cols;
            int row_size = Uy->rows;
    fclose(pf);
#endif // 0
/*---------------------------------------------------------------*/
    /* matlab : bx = bx + gamma*(Ux-Wx)
                by = by + gamma*(Uy-Wy)  */

    for(c = 0; c <Ux->cols; c++){
        for(r = 0; r < Ux->rows; r++){
            float elements =cvmGet(bx,r,c) + gamma*(cvmGet(Ux,r,c) - cvmGet(Wx,r,c));
            cvmSet(bx,r,c,elements);
            float elements_aa = cvmGet(by,r,c) + gamma*(cvmGet(Uy,r,c) - cvmGet(Wy,r,c));
            cvmSet(by,r,c,elements_aa);
        }
    }
    /* matlab  : if loop */
    if(aL1 > 0){
        for(c = 0; c <d->cols; c++){
            for(r = 0; r < d->rows; r++){
                float elements =cvmGet(d,r,c) + gamma*(cvmGet(PsiTU,r,c) - cvmGet(Z,r,c));
                cvmSet(d,r,c,elements);
            }
        }
    }

    }//for

    /* matlab : U = U/fctr*/
  //  int p,o;
    if(opts_para.normalize){
        for(c = 0; c <U->cols; c++){
            for(r = 0; r < U->rows; r++){
                float elements = cvmGet(U,r,c) / fctr;
                cvmSet(U, r, c, elements);
            }
        }

    }

    if(opts_para.real_sol){
        for(c = 0; c <U->cols; c++){
            for(r = 0; r < U->rows; r++){
                double complex uu = conj(cvmGet(U, r, c));
                //float elements = cvmGet(U,r,c);
                cvmSet(U,r,c,creal(uu));
            }
        }
    }



/* debug ------------------------------------------------------------*/
#if 0
    FILE * pf;
    pf = fopen("out_dst.txt", "w");
    for(r = 0; r<U->rows; r++){
        for(c = 0; c<U->cols; c++){
            float elements = cvmGet(U,r,c);
            //printf("%f ",elements);
            fprintf(pf, "%f ", elements,0);

            if(c == U->cols-1){
                fprintf(pf,"\n ",0);
            }
        }
    }
            int col_size = U->cols;
            int row_size = U->rows;
    fclose(pf);
#endif // 0
/*---------------------------------------------------------------*/

    IplImage *conv_img = cvCreateImage(cvSize(U->cols,U->rows), IPL_DEPTH_8U, 1);
    /* CvMat to IplImage */
    for(c = 0; c <U->cols; c++){
        for(r = 0; r < U->rows; r++){
            double elements =cvmGet(U,r,c)*255.0;
            cvSetReal2D(conv_img,c,r,elements);

        }
    }
    //double a_min = 0.0;
    //double a_max = 255.0;
    //cvNormalize(conv_img,conv_img,255,0,CV_C,0);


    cvShowImage("U", conv_img);
    cvWaitKey(0);
    cvDestroyWindow("U");


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
    //cvReleaseMat(&UUx);
   // cvReleaseMat(&UUy);
    //cvReleaseMat(&UU);
    //cvReleaseMat(&V);
     /*compute*/
    //cvReleaseMat(&fft2_U );

    //cvReleaseMat(&Wx);
    //cvReleaseMat(&Wy);
    //cvReleaseMat(&DxtU_INDEX );
    //cvReleaseMat(&DytU_INDEX);

   // cvReleaseMat(&DxtU_INDEX_01);
    //cvReleaseMat(&DxtU_INDEX_02 );
   // cvReleaseMat(& DxtU );
   //cvReleaseMat(&DytU_INDEX_01 );
    //cvReleaseMat(&DytU_INDEX_02 );
    //cvReleaseMat(&DytU );

   //cvReleaseMat(&rhs);
   //cvReleaseMat(&fft2_rhs);

   //cvReleaseMat(&Z);

   // cvReleaseMat(&Uprev);

}
