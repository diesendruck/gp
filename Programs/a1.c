#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "mex.h" 

void sub(int n,int d,double *x,double *f);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int n,d;
	double *x,*f;
	n = (int)(mxGetScalar(prhs[0]));
	d = (int)(mxGetScalar(prhs[1]));
	x = mxGetPr(prhs[2]);
	plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
	f = mxGetPr(plhs[0]);
	sub(n,d,x,f);
}

void sub(int n,int d,double *x,double *f)
{
 int i = 0,k = 0;
 for (i=0;i<n;i++)
 {
     f[i] = 0;
     for (k=0;k<d;k++)
     {
         f[i] = f[i]+x[i*d+k]*x[i*d+k];
     }
 }
}