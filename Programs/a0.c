#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "mex.h" 

void fun(int n,int d,double *x,double *a);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int n,d;
	double *x,*a;
	n = (int)(mxGetScalar(prhs[0]));
	d = (int)(mxGetScalar(prhs[1]));
    x = mxGetPr(prhs[2]);
	plhs[0] = mxCreateDoubleMatrix(n*n,n*d,mxREAL);
	a = mxGetPr(plhs[0]);
	fun(n,d,x,a);
}

void fun(int n,int d,double *x,double *a)
{
 int i=0,j=0,k=0;
 for (i=0;i<n;i++)
 {
    for (j=0;j<n;j++)
    {
        for (k=0;k<d;k++)
        {
            *(a+i*n+j*n*n*d+(k*n*n+j)) = *(x+j*d+k) - *(x+i*d+k);
        }
    }
 }
}