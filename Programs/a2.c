#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "mex.h" 

void sub(int n,int d,double *f,double *eps,double *y);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int n,d;
	double *eps,*f,*y;
	n = (int)(mxGetScalar(prhs[0]));
	d = (int)(mxGetScalar(prhs[1]));
	f = mxGetPr(prhs[2]);
    eps = mxGetPr(prhs[3]);
	plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
	y = mxGetPr(plhs[0]);
	sub(n,d,f,eps,y);
}

void sub(int n,int d,double *f,double *eps,double *y)
{
 int i = 0,k = 0;
 for (i=0;i<n;i++)
 {
     y[i] = 0;
     for (k=0;k<d;k++)
     {
         y[i] = f[i]+eps[i];
     }
 }
}