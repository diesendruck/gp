#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "mex.h" 

void sub(int n,int d,int r,double *del,double *sum);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int n,d,r;
	double *del,*sum;
	n = (int)(mxGetScalar(prhs[0]));
	d = (int)(mxGetScalar(prhs[1]));
    r = (int)(mxGetScalar(prhs[2]));
	del = mxGetPr(prhs[3]);
	plhs[0] = mxCreateDoubleMatrix(d,d,mxREAL);
	sum = mxGetPr(plhs[0]);
	sub(n,d,r,del,sum);
}

void sub(int n,int d,int r,double *del,double *sum)
{
 int i = 0,k = 0,l = 0;
 for (k=0;k<d;k++)
     for (l=0;l<d;l++)
         *(sum+k*d+l)=0;
 for (i=0;i<n;i++)
     for (k=0;k<d;k++)
         for (l=0;l<d;l++)
             *(sum+d*k+l) = *(sum+d*k+l)+del[i*n*d+r*d+k]*del[i*n*d+r*d+l];
}