/* Include files */
#include <cv.h>

#define pi 3.141592

int round_x(double a)
{
    int round_x;

    if(a < 0)  round_x = (int)(a * 100 + 0.5);
    else   round_x  = (int)(a * 100 + -0.5);

}
void LineMask_limitedandleTange(int L, int N, CvMat* M, CvMat* Mh, CvMat* mi, CvMat* mhi)
{
    /*
    matlab : thc = linspace(0, pi-pi/L, L);
    */
    CvMat* thc = cvCreateMat(1,L,CV_32FC1);
    CvMat *xc = cvCreateMat(1,L,CV_32FC1);
    CvMat *yr = cvCreateMat(1,L,CV_32FC1);

    int i,j;
    for(i = 1; i< thc->rows; i++){
        cvmSet(thc,0,1,0);
        cvmSet(thc,i,1,i*(pi-pi/L)/(L-1));
    }

    //matlab : M = zeros(N);
    M = cvCreateMat(N,N,CV_32FC1);
    cvmSetZero(M);

    int ll, nn;
    int blank = -N/2+1;
    CvMat *blank_aa = cvCreateMat(1, N,CV_32FC1 );

    for(ll = 0; ll < N; ll++){
        cvmSet(blank_aa, ll,1,blank+ll);
    }

    for(ll = 0; ll < L; ll++){
        if(cvmGet(thc,ll,1)>= pi/3 && cvmGet(thc,ll,1)<= (pi/2+pi/18) ||
            cvmGet(thc,ll,1) >=(2*pi/3-pi/18) && cvmGet(thc,ll,1)<= 2*pi/3){

            cvmSet(xc, ll,1,round_x(tan(cvmGet(thc,ll,1)) * cvmGet(blank_aa,ll,1)) + N/2+1);  //round function

            for(nn=0; nn<N-1; nn++){
                cvmSet(M,cvmGet(xc,1,nn),nn+1,0);
            }
        }//if
        else{
            if(cvmGet(thc,ll,1) <= pi/4 || cvmGet(thc,ll,1)>3*pi/4){
                cvmSet(yr, ll,1,round_x(tan(cvmGet(thc,ll,1)) * cvmGet(blank_aa,ll,1)) + N/2+1);
                for(nn = 1; nn<N-1; nn++){
                    cvmSet(M,nn+1,cvmGet(yr,1,nn),1);
                }
            }
            else{
                cvmSet(xc, ll,1,round_x(tan(cvmGet(thc,ll,1)) * cvmGet(blank_aa,ll,1)) + N/2+1);  //round function
                for(nn=1; nn<N-1; nn++){
                      cvmSet(M,cvmGet(xc,1,nn),nn+1,1);
                }
            }
        }//else
    }//for


    //upper half plane mask
    Mh = cvCreateMat(N,N,CV_32FC1);
    cvmSetZero(Mh);
    cvmCopy(M,Mh);

    for(i = N/2+2; i<N; i++){
        for(j = 0; j<N; j++){
            cvmSet(Mh, j, i, 0);
        }
    }
    for(j = N/2+1; j<N; j++){
        cvmSet(Mh, j, N/2+2, 0);
    }

    cvFFT(M,M,1,1); //0:fft, 1:ifft
    int count = 0;
    for(i = 0; i<M->cols; i++){
        for(j = 0; j<M->rows; j++){
            if(cvmGet(M,j,i) != 0){
                count=count+1;
                cvmSet(mi,1,count, j+i*M->cols);    //find the nonzero elements in matrix
            }
        }
    }
    count = 0;

    cvFFT(Mh,Mh,1,1);
    for(i = 0; i<Mh->cols; i++){
        for(j = 0; j<Mh->rows; j++){
            if(cvmGet(Mh,j,i) != 0){
                count=count+1;
                cvmSet(mhi,1,count, j+i*Mh->cols);    //find the nonzero elements in matrix
            }
        }
    }
}

