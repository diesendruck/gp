#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "mex.h" 

void sub(int n,double *eta,double *th,double *del,double *xi,double *beta,int d);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int n,d;
	double *eta,*th,*del,*xi,*beta;
	n = (int)(mxGetScalar(prhs[0]));
	eta = mxGetPr(prhs[1]);
	th = mxGetPr(prhs[2]);
	del = mxGetPr(prhs[3]);
	xi = mxGetPr(prhs[4]);
       d = (int)(mxGetScalar(prhs[5]));
	plhs[0] = mxCreateDoubleMatrix(n*n,1,mxREAL);
	beta = mxGetPr(plhs[0]);
	sub(n,eta,th,del,xi,beta,d);
}

void sub(int n,double *eta,double *th,double *del,double *xi,double *beta,int d)
{
 int i = 0,j = 0,k = 0;
 double val = 0;
 for (i=0;i<n;i++)
 {
    for (j=0;j<n;j++)
    {
	*(beta+i*n+j) = 0;
       val = 0;
	for (k=0;k<d;k++)
	 val = val+del[i*n*d+j*d+k]*xi[j*d+k];
	*(beta+i*n+j) = (*(eta+i*n+j)) - (*(th+j) + val - *(th+i)) ;
    }
 }
}