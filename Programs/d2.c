#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "mex.h" 

void fun(int n,double *eta,double *del,double *xi,double *eta_T,int d);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int n,d;
	double *eta,*del,*xi,*eta_T;
	n = (int)(mxGetScalar(prhs[0]));
	eta = mxGetPr(prhs[1]);
	del = mxGetPr(prhs[2]);
	xi = mxGetPr(prhs[3]);
	d = (int)(mxGetScalar(prhs[4]));
	plhs[0] = mxCreateDoubleMatrix(n*n,1,mxREAL);
	eta_T = mxGetPr(plhs[0]);
	fun(n,eta,del,xi,eta_T,d);
}

void fun(int n,double *eta,double *del,double *xi,double *eta_T,int d)
{
 int i = 0,j = 0,k = 0;
 double val=0;
 for (i=0;i<n;i++)
 {
    for (j=0;j<n;j++)
    {
 	 eta_T[i*n+j] = 0;
	 val = 0;
        for (k=0;k<d;k++)
          val = val+del[i*n*d+j*d+k]*xi[j*d+k];
        eta_T[i*n+j] = eta[i*n+j] - val;
    }
 }
}