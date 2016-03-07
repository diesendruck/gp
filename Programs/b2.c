#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "mex.h" 

void fun(int n,int d,double rho,double *nu,double *eta, double *del,double *th,double *Y,double *xi,double *t);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int n,d;
	double rho,*nu,*eta,*del,*th,*t,*xi,*Y;
	n = (int)(mxGetScalar(prhs[0]));
	d = (int)(mxGetScalar(prhs[1]));
	rho = mxGetScalar(prhs[2]);
	nu = mxGetPr(prhs[3]);
    eta = mxGetPr(prhs[4]);
    del = mxGetPr(prhs[5]);
	th = mxGetPr(prhs[6]);
    Y = mxGetPr(prhs[7]);
    xi = mxGetPr(prhs[8]);
	plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
	t = mxGetPr(plhs[0]);
	fun(n,d,rho,nu,eta,del,th,Y,xi,t);
}

void fun(int n,int d,double rho,double *nu,double *eta, double *del,double *th,double *Y,double *xi,double *t)
{
 int i = 0,j = 0,k = 0;
 double nu_row[n],nu_col[n],eta_row[n],eta_col[n],deldot_row[n],deldot_col[n],th_sum[n];
 for (i=0;i<n;i++)
 {
    nu_row[i] = 0;
    nu_col[i] = 0;
    eta_row[i] = 0;
    eta_col[i] = 0;
    deldot_row[i] = 0;
    deldot_col[i] = 0;
    for (j=0;j<n;j++)
    {
     nu_row[i] = nu_row[i] + nu[i*n+j];
     nu_col[i] = nu_col[i] + nu[j*n+i];
     eta_row[i] = eta_row[i] + eta[i*n+j];
     eta_col[i] = eta_col[i] + eta[j*n+i];
     for (k=0;k<d;k++)
     {
         deldot_row[i] = deldot_row[i] + del[i*n*d+j*d+k]*xi[j*d+k];
         deldot_col[i] = deldot_col[i] + del[j*n*d+i*d+k]*xi[i*d+k];
     }
    }
  }
 for (i=0;i<n;i++)
 {
    th_sum[i] = 0;
    for (j=0;j<n;j++)
    {
        if (i!=j)
        {
            th_sum[i] = th_sum[i]+th[j];
        }
    }
    t[i] = (1/(1+2*n*rho-2*rho)) * (Y[i]+nu_col[i]-nu_row[i]+rho*(2*th_sum[i]+eta_col[i]-eta_row[i]-deldot_col[i]+deldot_row[i]));
 }
}