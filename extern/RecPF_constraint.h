#include <stdio.h>
//#include "enhance.h"

void RecPF_constraint(int m,int n,double aTV, double aL1,CvMat* pick,CvMat* B,
                      int TVtype,struct opts_para,CvMat* WT,CvMat* W,int range,IplImage* I,int constraint, CvMat* U);
