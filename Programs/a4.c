#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "mex.h" 

void sub(int n,int d,double *del,double *inv);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int n,d;
	double *del,*inv;
	n = (int)(mxGetScalar(prhs[0]));
	d = (int)(mxGetScalar(prhs[1]));
	del = mxGetPr(prhs[2]);
	plhs[0] = mxCreateDoubleMatrix(n*d*d,1,mxREAL);
	inv = mxGetPr(plhs[0]);
	sub(n,d,del,inv);
}

void sub(int n,int d,double *del,double *inv)
{
 int i = 0,j = 0,k = 0,l = 0;
 double val[500][500];
 for (j=0;j<n;j++)
 {
     for (k=0;k<d;k++)
         for (l=0;l<d;l++)
             val[k][l]=0;
     for (i=0;i<n;i++)
         for (k=0;k<d;k++)
             for (l=0;l<d;l++)
                 val[k][l] = val[k][l]+del[i*n*d+j*d+k]*del[i*n*d+j*d+l];
     for (k=0;k<d;k++)
         for (l=0;l<d;l++)
             inv[j*d*d+k*d+l] = 1/val[k][l];
 }      
}