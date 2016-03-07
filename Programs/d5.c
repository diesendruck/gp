#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "mex.h" 

void sub(int n,double rho,double *nu,double *beta,double *nu_1);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int n;
	double rho,*nu,*beta,*nu_1;
	n = (int)(mxGetScalar(prhs[0]));
	rho = mxGetScalar(prhs[1]);
	nu = mxGetPr(prhs[2]);
	beta = mxGetPr(prhs[3]);
	plhs[0] = mxCreateDoubleMatrix(n*n,1,mxREAL);
	nu_1 = mxGetPr(plhs[0]);
	sub(n,rho,nu,beta,nu_1);
}

void sub(int n,double rho,double *nu,double *beta,double *nu_1)
{
 int i = 0,j = 0;
 for (i=0;i<n;i++)
 {
    for (j=0;j<n;j++)
    {
	*(nu_1+i*n+j) = 0;
	*(nu_1+i*n+j) = *(nu+i*n+j) + rho * (*(beta+i*n+j));
    }
 }
}