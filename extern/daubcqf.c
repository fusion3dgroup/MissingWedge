#include <cv.h>
//#include "daubcqf.h"

void daubcqf(int N, CvMat * hH){


    int K = N/2, a=1, p=1, i,j;

    /* matlab : for loop */

    CvMat * h_0;
    int c;
    for(c = 0; c< 2; c++){
        cvmSet(h_0, c, 0, 1.);
    }

    for(j=1; j< K-1; j++){
        a = -a*0.25*(j+K-1)/j;

    }

    /*   matlab : q = sort(roots(q));
                 qt = q(1:K-1); */

    CvMat* q, *q_roots, *q_sort, *qt;
    cvSolvePoly(q,q_roots,20,100);
    cvSort(q_roots,q_sort,0,0);

    for(i = 0; i<K-1; i++){
        cvmSet(qt, i,1,cvmGet(q_sort,i,1));
    }

    // poly(qt) find coefficient qt vector

    /*
       matlab : h_0 = conv(h_0,real(poly(qt));
                h_0 = sqrt(2)*h_0/sum(h_0);
    */
    cvFilter2D(hH, hH,qt,cvPoint(-1,-1));
    double sum_elements=0;
    for(i=0; i<hH->cols; i++){
        for(j = 0; j< hH->rows; j++){
            double elements = cvmGet(hH, j, i);
            sum_elements = sum_elements + elements;
        }
    }
    for(i=0; i<hH->cols; i++){
        for(j = 0; j< hH->rows; j++){
            cvmSet(hH,j,i,sqrt(2) * cvmGet(hH,j,i)/sum_elements);
        }
    }

    /*
        matlab : h_1 = rot90(h_0,2) //rotation 180
    */
    CvMat *hH2 = cvCreateMat(hH->rows,hH->cols, CV_32FC1);
    for(i=0; i<hH->cols; i++){
        for(j = 0; j< hH->rows; j++){
            cvmSet(hH2, hH->rows-j, hH->cols-i,cvmGet(hH, j, i));
        }
    }

    /*
        matlab : h_1(1:2:N) = -h_1(1:2:N)
    */
    int count = 0 ;
    for(i=0; i<hH2->rows; i++){
        for(j = 0; j< hH2->cols; j++){
            count = count+1;
            if(count %2 == 1){
                cvmSet(hH2,i,j,-cvmGet(hH2, i,j));
            }
            else{
                cvmSet(hH2,i,j,cvmGet(hH2, i,j));
            }

        }
    }

    cvReleaseMat(&qt);
    cvReleaseMat(&q_sort);
    cvReleaseMat(&q_roots);
    cvReleaseMat(&q);
    cvReleaseMat(&hH2);
    cvReleaseMat(&hH);
}
