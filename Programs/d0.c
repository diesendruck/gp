#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "mex.h" 

void fun(int n,double *D);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int n;
	double *D;
	n = (int)(mxGetScalar(prhs[0]));
	plhs[0] = mxCreateDoubleMatrix(n*n,n,mxREAL);
	D = mxGetPr(plhs[0]);
	fun(n,D);
}

void fun(int n,double *D)
{
 int i=0,j=0;
 for (i=0;i<n;i++)
 {
    for (j=0;j<n;j++)
    {
        *(D+i*n*(n+1)+j) = 1;
        *(D+i*n*n+j*n+i) = *(D+i*n*n+j*n+i) - 1;
    }
 }
}