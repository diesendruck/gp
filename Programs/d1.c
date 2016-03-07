#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "mex.h" 

void sub(int n,int d,double rho,double *nu,double *eta,double *th,double *del,double *invdel,double *xi);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int n,d;
	double rho,*nu,*eta,*th,*del,*invdel,*xi;
	n = (int)(mxGetScalar(prhs[0]));
	d = (int)(mxGetScalar(prhs[1]));
	rho = mxGetScalar(prhs[2]);
	nu = mxGetPr(prhs[3]);
	eta = mxGetPr(prhs[4]);
	th = mxGetPr(prhs[5]);
	del = mxGetPr(prhs[6]);
	invdel = mxGetPr(prhs[7]);
	plhs[0] = mxCreateDoubleMatrix(n*d,1,mxREAL);
	xi = mxGetPr(plhs[0]);
	sub(n,d,rho,nu,eta,th,del,invdel,xi);
}

void sub(int n,int d,double rho,double *nu,double *eta,double *th,double *del,double *invdel,double *xi)
{
 int i = 0,j = 0,k = 0,l = 0;
 double val[500];
 for (j=0;j<n;j++)
 {
  for (k=0;k<d;k++)
   val[k]=0;
  for (i=0;i<n;i++)
  {
   for (k=0;k<d;k++)
   {
    val[k] = val[k] + ((1/rho) * nu[i*n+j] + eta[i*n+j] - (th[j] - th[i])) * del[i*n*d+j*d+k];
   }
  }
  for (k=0;k<d;k++)
  {
    xi[j*d+k] = 0;
   for (l=0;l<d;l++)
    xi[j*d+k] = xi[j*d+k] + val[l] * invdel[j*d*d+l*d+k];
  }
 }
}