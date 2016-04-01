#include <matrix.h>
#include <cmath>
#include "mex.h"
#include <omp.h>
//#include "/usr/llvm-gcc-4.2/lib/gcc/i686-apple-darwin11/4.2.1/include/omp.h"

/* Definitions to keep compatibility with earlier versions of ML */
#ifndef MWSIZE_MAX
typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex;

#if (defined(_LP64) || defined(_WIN64)) && !defined(MX_COMPAT_32)
/* Currently 2^48 based on hardware limitations */
# define MWSIZE_MAX    281474976710655UL
# define MWINDEX_MAX   281474976710655UL
# define MWSINDEX_MAX  281474976710655L
# define MWSINDEX_MIN -281474976710655L
#else
# define MWSIZE_MAX    2147483647UL
# define MWINDEX_MAX   2147483647UL
# define MWSINDEX_MAX  2147483647L
# define MWSINDEX_MIN -2147483647L
#endif
#define MWSIZE_MIN    0UL
#define MWINDEX_MIN   0UL
#endif
typedef unsigned int uint;
/* [U_r,U_l,] = WENO_mex(U,dx,k,B,C,D,Ct,Dt,N)*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    const mxArray *psi_p, *w_p,*b_p, *N_p;
    const mwSize *dims;
    mxArray *u_p;
    double *psi, *w,*b, *N;
    double *u;
    int nq,nx,nmom;
//     int i,j,k;
    uint nThreads, tid, k;
//     double IntY, IntX, IntFull;
    
    
    if ( nrhs != 4 ) {
        mexErrMsgIdAndTxt("MATLAB:projectBasis:rhs","This function takes 4 input arguments.");
    }
    
//associate inputs
    
    psi_p = prhs[0]; //nq x nx
    w_p = prhs[1];  // nq x 1
    b_p = prhs[2];  // nmom x nq
    N_p  = prhs[3];
    
//figure out dimensions
    nq = mxGetM(psi_p);
    nx = mxGetN(psi_p);
    nmom = mxGetM(b_p);
    
//     mexPrintf("nq = %d\n",nq);
//     mexPrintf("nx = %d\n",nx);
//     mexPrintf("nmom = %d\n",nmom);
    
    
//associate pointers
    psi = mxGetPr(psi_p);
    w = mxGetPr(w_p);
    b = mxGetPr(b_p);
    
    N = mxGetPr(N_p);
    
    //associate outputs
    u_p = plhs[0] = mxCreateDoubleMatrix(nmom,nx,mxREAL);
    u = mxGetPr(u_p);
    int NumberThreads = (int) N[0]+0.5;
//     
// #pragma omp parallel private(tid) shared(nThreads) num_threads(NumberThreads)
//     {
//         tid = omp_get_thread_num();
//         
//         if (tid==0) {
//             nThreads = omp_get_num_threads();
//             
//         }
//         
//     }
//     mexPrintf("nThreads = %i\n",nThreads);
//     return;

    #pragma omp parallel for shared(nx,nmom,nq,u,w,psi,b) default(none) schedule(static) if(nx>1000 && NumberThreads>1) num_threads(NumberThreads)
    for(int i=0;i<nx;i++)
    {
        for(int n=0;n<nmom;n++)
        {
            for(int j=0;j<nq;j++)
            {
                u[n+nmom*i]+=w[j]*psi[j+nq*i]*b[n+nmom*j];
            }
        }
        
    }
    
    return;
}