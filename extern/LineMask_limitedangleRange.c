/* Include files */
#include <cv.h>
#include <math.h>
#define pi 3.141592

#define round_ban(a) (floor(a)*pow(10,0)+0.5)/pow(10,0)

void LineMask_limitedandleTange(int L, int N,  CvMat* mhi, CvMat* M, CvMat* Mh, CvMat* mi)
{

    /* matlab : thc = linspace(0, pi-pi/L, L); */

    CvMat *thc = cvCreateMat(1,L,CV_32F);

    int i,j;
    cvmSet(thc,0,0,0.);

    for(i = 1; i< L; i++){
        float elements = i*(pi-pi/L)/(L-1);
        cvmSet(thc,0,i,elements);
    }

    /*matlab : M = zeros(N);*/
    M = cvCreateMat(N,N,CV_32F);
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
    Mh = cvCreateMat(N,N,CV_32FC1);
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

    cvFFT(M,M,1,1); //0:fft, 1:ifft

    /* matlab : mMH = find(M)*/
    int count = 0;
    for(c = 0; c< M->cols; c++){
        for(r = 0; r< M->rows; r++){
            if(cvmGet(M,r,c) != 0){
                count=count+1;
            }
        }
    }
    int a_blank= count+1;
    mi = cvCreateMat(a_blank,1,CV_32F);

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

    cvFFT(Mh,Mh,1,1);

    /* matlab : mMH = find(M)*/
    count = 0;
    for(c = 0; c< Mh->cols; c++){
        for(r = 0; r< Mh->rows; r++){
            if(cvmGet(Mh,r,c) != 0){
                count=count+1;
            }
        }
    }

    a_blank= count;
    mhi = cvCreateMat(a_blank+1,1,CV_32F);
    count = 0;
    for(c = 0; c< Mh->cols; c++){
        for(r = 0; r< Mh->rows; r++){
            if(cvmGet(Mh,r,c) != 0){
                float elements = cvmGet(Mh,r,c);
                cvmSet(mhi,count, 0, elements);    //find the nonzero elements in matrix
                count=count+1;
            }
        }
    }

    cvReleaseMat(&thc);
    cvReleaseMat(&xc);
    cvReleaseMat(&yr);
    cvReleaseMat(&blank_aa);
}

