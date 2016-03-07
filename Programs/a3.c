#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "mex.h" 

void sub(int n,int d,double *x,double *del);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int n,d;
	double *x,*del;
	n = (int)(mxGetScalar(prhs[0]));
	d = (int)(mxGetScalar(prhs[1]));
	x = mxGetPr(prhs[2]);
	plhs[0] = mxCreateDoubleMatrix(n*n*d,1,mxREAL);
	del = mxGetPr(plhs[0]);
	sub(n,d,x,del);
}

void sub(int n,int d,double *x,double *del)
{
 int i = 0,j = 0,k = 0;
 for (i=0;i<n;i++)
 {
     for (j=0;j<n;j++)
     {
         for (k=0;k<d;k++)
         {
             del[i*n*d+j*d+k] = x[i*d+k]-x[j*d+k];
         }
     }
 }
}