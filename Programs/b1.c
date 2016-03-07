#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "mex.h" 

void sub(int n,double *nu,double *v);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int n;
	double *nu,*v;
	n = (int)(mxGetScalar(prhs[0]));
    nu = mxGetPr(prhs[1]);
	plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
	v = mxGetPr(plhs[0]);
	sub(n,nu,v);
}

void sub(int n,double *nu,double *v)
{
 int i = 0,j = 0;
 double val=0;
 for (i=0;i<n;i++)
 {
     val = 0;
     for (j=0;j<n;j++)
     {
         val = val+nu[i*n+j];
         val = val-nu[j*n+i];
     }
	 v[i] = val;
 }
}