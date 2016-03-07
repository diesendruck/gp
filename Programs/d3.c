#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "mex.h" 

void sub(int n,double rho,double *nu,double *th,double *del,double *xi,double *eta,int d,double tol);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int n,d;
	double rho,tol,*nu,*th,*del,*xi,*eta;
	n = (int)(mxGetScalar(prhs[0]));
	rho = mxGetScalar(prhs[1]);
	nu = mxGetPr(prhs[2]);
	th = mxGetPr(prhs[3]);
	del = mxGetPr(prhs[4]);
	xi = mxGetPr(prhs[5]);
	d = (int)(mxGetScalar(prhs[6]));
    tol = mxGetScalar(prhs[7]);
	plhs[0] = mxCreateDoubleMatrix(n*n,1,mxREAL);
	eta = mxGetPr(plhs[0]);
	sub(n,rho,nu,th,del,xi,eta,d,tol);
}

void sub(int n,double rho,double *nu,double *th,double *del,double *xi,double *eta,int d,double tol)
{
 int i = 0,j = 0,k = 0;
 double val = 0;
 for (i=0;i<n;i++)
 {
    for (j=0;j<n;j++)
    {
	*(eta+i*n+j) = 0;
       val = 0;
	for (k=0;k<d;k++)
         val = val+del[i*n*d+j*d+k]*xi[j*d+k];
	*(eta+i*n+j) = *(th+j) + val - *(th+i) - *(nu+i*n+j)/rho;
	if (*(eta+i*n+j)>-tol)
		*(eta+i*n+j)=-tol;
    }
 }
}