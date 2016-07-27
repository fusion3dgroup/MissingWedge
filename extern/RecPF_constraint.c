#include <stdio.h>
#include <cv.h>
#include <math.h>
//#include "RecPF_constraint.h"
#include "enhance.h"

void RecPF_constraint(int m,int n,double aTV, double aL1,CvMat* picks,CvMat* B,
                      int TVtype,struct OPTS opts_para ,CvMat* PsiT,CvMat* Psi,double range,IplImage* I,int constraint, CvMat* U)
{

    int bPrint = 0;
    int bComplex = 1;
    double fctr = 1/range;
    cvSetZero(U);

 //   CvMat* Ux, *Uy;


    if(/*OPTS.normalize == 1*/1) //true
    {

        //B = fctr*B
        int j,i;
        for(j = 1; j<B->cols; j++){
            for(i = 1; i<B->rows; i++){
                B->data.db[i*B->step+j] = fctr * B->data.db[i*B->step+j];
            }
        }

        int count =0;
        for(j = 1; j<picks->cols; j++){
            for(i = 1; i<picks->rows; i++){
                int element = cvmGet(picks,i,j);
                if (element == 0){
                    count++;
                }

            }
        }
        aTV = count/sqrt(m*n)*aTV;
        aL1 = count/sqrt(m*n)*aL1;
    }

    CvMat *Numer1, *Denom1, *Denom2, *Denom;
    cvZero(Numer1);

    int j, i;
    for(j = 1; j<B->cols; j++){
        for(i = 1; i<B->rows; i++){
            Numer1->data.db[i*Numer1->step+j] = sqrt(m*n) * Numer1->data.db[i*Numer1->step+j];
            Denom1->data.db[i*Denom1->step+j] = 1;
        }
    }

    double prd = sqrt(aTV * /*beta*/2.0);
    //denom2


    if(aL1 == 0){
        cvmAdd(Denom1, Denom2, Denom);
    }
    else{
        cvmAdd(Denom1, Denom2, Denom);
        for(j = 1; j<Denom->cols; j++){
            for(i = 1; i<Denom->rows; i++){
                double element = cvmGet(Denom, i, j) + aL1*/*beta*/2.0;
                cvmSet(Denom,i, j, element);
            }
        }
    }

    CvMat* Ux, *Uy, *bx, *by, *d, *d_penalty, *PsiTU;
    cvSetZero(Ux);
    cvSetZero(Uy);
    cvSetZero(bx);
    cvSetZero(by);
    cvSetZero(d);
    cvSetZero(d_penalty);

    if(aL1 > 0){
        cvSetZero(PsiTU);
        cvSetZero(d);
    }


    double MaxItr;
    double beta = 2.0;
    int ii;
    CvMat* Wx, *Wy;
    CvMat* rhs, *Z;
    CvMat* blank, *fft2_U, *fft2_rhs;

        for(ii =0; ii<MaxItr; ii++){

            switch(TVtype){
                case 1 :
                    cvmAdd(Ux, bx, Ux);
                    cvmAdd(Uy, by, Uy);
                break;

                case 2 :
                    Compute_Wx_Wy(Ux, Uy, bx, by, 1/beta, &Wx, &Wy);
                break;
                default :
                    printf("TVtype must be 1 or 2\n");
                break;

            }
            if(aL1 > 0){

                cvmAdd(PsiTU,d,PsiTU);
                int i,j;
                for(j = 1; j<PsiTU->cols; j++){
                    for(i = 1; i<PsiTU->rows; i++){

                        double element_sin = sin(cvmGet(PsiTU, i, j));
                        double element_abs = abs(cvmGet(PsiTU, i, j));

                        if(element_abs-1/beta > 0){
                            cvmSet(PsiTU,i,j,element_abs*element_sin);
                        }
                        else{
                            cvmSet(PsiTU,i,j,0.);
                        }

                    }
                }

            }

            Compute_rhs_DxtU_DytU(Wx, Wy, bx, by, aTV*beta, &rhs);

        if(aL1 > 0){

            for(j = 1; j<Psi->cols; j++){
                for(i = 1; i<Psi->rows; i++){

                    cvmSet(Psi, i, j, cvmGet(Z, i, j) - cvmGet(d, i, j));
                    double element = cvmGet(Psi, i, j) * aL1* beta;
                    cvmSet(d_penalty, i, j, element);
                }
            }
            cvmAdd(rhs, d_penalty, rhs);
        }
        //fft2_U = (Numer1 + fft2(rhs))./Denom;    cvFFT(image,image, 0,0);

        cvFFT(rhs, fft2_rhs, 0, 0);
        cvFFT(rhs, rhs, 0, 0);
        cvmAdd(Numer1, rhs, blank);
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

        Compute_Ux_Uy(U, &Ux, &Uy);

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
}
