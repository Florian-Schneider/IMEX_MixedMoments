/*=========================================================
 * matrixDivide.c - Example for illustrating how to use
 * LAPACK within a C MEX-file.
 *
 * X = matrixDivide(A,B) computes the solution to a
 * system of linear equations A * X = B
 * using LAPACK routine DGESV, where
 * A is an N-by-N matrix
 * X and B are N-by-1 matrices.
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2009-2010 The MathWorks, Inc.
 *=======================================================*/
/* $Revision: 1.1.6.3 $ $Date: 2011/01/28 18:11:56 $ */

#include "mex.h"
#include "lapack.h"
#include "blas.h"
#include "matrix.h"
#include <string.h>
#include <stdlib.h>
/* #include <math.h>*/
#include <limits.h>
#include <e_exp.c>




double signum(double x)
{
    if(x>0.0)
        return 1.0;
    else
        return -1.0;
    
}


void PrintMultipliers(mxArray *mult[],ptrdiff_t *nq, ptrdiff_t *NN)
{
    mxArray *alpha_p = mult[0] ;
    mxArray *beta_p = mult[1];
    mxArray *T_p = mult[2];
    mxArray *p_p = mult[3];
    mxArray *wG_p = mult[4];
    double *alpha = mxGetPr(alpha_p);
    double *beta = mxGetPr(beta_p);
    double *T = mxGetPr(T_p);
    double *p = mxGetPr(p_p);
    double *wG = mxGetPr(wG_p);
    ptrdiff_t i1,i2;
    mexPrintf("alpha =\n\t\t");
    for(i1=0;i1< *NN;i1++)
    {
        mexPrintf("%f\t",alpha[i1]);
    }
    mexPrintf("\n\n");
    
    mexPrintf("beta =\n\t\t");
    for(i1=0;i1< *NN;i1++)
    {
        mexPrintf("%f\t",beta[i1]);
    }
    mexPrintf("\n\n");
    
    mexPrintf("T =\n\t\t");
    for(i1=0;i1< *NN;i1++)
    {
        for(i2=0;i2< *NN;i2++)
        {
            mexPrintf("%f\t",T[i2*(*NN)+i1]);
        }
        mexPrintf("\n\t\t");
    }
    mexPrintf("\n\n");
    
    mexPrintf("p =\n\t\t");
    for(i1=0;i1< *NN;i1++)
    {
        for(i2=0;i2< *nq;i2++)
        {
            mexPrintf("%f\t",p[i2*(*NN)+i1]);
        }
        mexPrintf("\n\t\t");
    }
    mexPrintf("\n\n");
    
    mexPrintf("wG =\n\t\t");
    for(i1=0;i1< *nq;i1++)
    {
        mexPrintf("%f\t",wG[i1]);
    }
    mexPrintf("\n\n");
}

double fobjP(mxArray *wG_p,mxArray *beta_p,mxArray *phi_p,ptrdiff_t *nq, ptrdiff_t *NN)
{
    double *beta = mxGetPr(beta_p);
    double *wG = mxGetPr(wG_p);
    double *phi = mxGetPr(phi_p);
    double f=0;
    ptrdiff_t i,j;
    
    for(j=0;j<*nq;j++)
    {
        f += wG[j];
    }
    for(i=0;i<*NN;i++)
    {
        f-=beta[i]*phi[i];
    }
    
    return f;
    
}

void gradP(mxArray *mult[],mxArray *phi_p,ptrdiff_t *nq, ptrdiff_t *NN,mxArray *g_p)
{
    mxArray *alpha_p = mult[0] ;
    mxArray *beta_p = mult[1];
    mxArray *T_p = mult[2];
    mxArray *p_p = mult[3];
    mxArray *wG_p = mult[4];
    double *alpha = mxGetPr(alpha_p);
    double *beta = mxGetPr(beta_p);
    double *T = mxGetPr(T_p);
    double *p = mxGetPr(p_p);
    double *wG = mxGetPr(wG_p);
    double *phi = mxGetPr(phi_p);
    double *g = mxGetPr(g_p);
    ptrdiff_t i,j;
    
/*     mexPrintf("g = \n\t\t");*/
    for(i=0;i<*NN;i++)
    {
        g[i] = 0;
        for(j=0;j<*nq;j++)
        {
            g[i]+=p[j*(*NN)+i]*wG[j];
        }
        g[i]-=phi[i];
/*        mexPrintf("%f\t",g[i]);*/
    }
/*     mexPrintf("\n\n");*/
}


void hessP(mxArray *mult[],ptrdiff_t *nq, ptrdiff_t *NN,mxArray *H_p,double* cond, ptrdiff_t *stat)
{
    mxArray *alpha_p = mult[0] ;
    mxArray *beta_p = mult[1];
    mxArray *T_p = mult[2];
    mxArray *p_p = mult[3];
    mxArray *wG_p = mult[4];
    double *alpha = mxGetPr(alpha_p);
    double *beta = mxGetPr(beta_p);
    double *T = mxGetPr(T_p);
    double *p = mxGetPr(p_p);
    double *wG = mxGetPr(wG_p);
    double *H = mxGetPr(H_p);
    ptrdiff_t i,i2,j;
    double anorm = 0; /*Calculate ||H||_1 = anorm for condition number*/
    double rcond = 1;
    mwSignedIndex info;
    double f=0;
    double *WORK;
    ptrdiff_t *IWORK;
    IWORK=(ptrdiff_t *) malloc((*NN)*sizeof(ptrdiff_t)); /*For condition number estimation*/
    WORK=(double *) malloc(3*(*NN)*sizeof(double));
    
    for(i=0;i<(*NN);i++)
    {
        for(i2=0;i2<=i;i2++) /*H is symmetric, and we only need the lower triangle for the cholesky decomp.*/
        {
            H[i+i2*(*NN)] = 0;
            for(j=0;j<(*nq);j++)
            {
                H[i+i2*(*NN)] += wG[j]*p[j*(*NN)+i]*p[j*(*NN)+i2];
            }
/*             mexPrintf("H[%d][%d] = %f\n",i,i2,H[i+i2*(*NN)]);*/
        }
    }
    
    
    for(i=0;i<(*NN);i++)
    {
        rcond = fabs(H[i+i*(*NN)]);
        for(i2=0;i2<i;i2++)
        {
            rcond+=fabs(H[i+i2*(*NN)]);
        }
        for(i2=i+1;i2<(*NN);i2++)
        {
            rcond+=fabs(H[i2+i*(*NN)]);
        }
        if(rcond>anorm)
            anorm = rcond;
    }
    rcond = 0;
    
    dpotrf("L", (ptrdiff_t *) NN, H, (ptrdiff_t *) NN, &info); /*cholesky decomposition, L is now saved in H*/
    
     /*for(i=0;i<(*NN);i++)
     {
         for(i2=0;i2<=i;i2++) H is symmetric, and we only need the lower triangle for the cholesky decomp.
         {
             mexPrintf("L[%d][%d] = %f\n",i,i2,H[i+i2*(*NN)]);
         }
     }*/
    
    if(info==0)
    {
         /*mexPrintf("Cholesky succeeded\n");*/
        *stat = 1;
        dpocon( "L", (ptrdiff_t *) NN, H,(ptrdiff_t *) NN, &anorm, &rcond, WORK, (ptrdiff_t *) IWORK, &info);
    
    }
    else
    {
/*         mexPrintf("Cholesky failed\n");*/
        *stat = 0;
        rcond = 1.0e-17;
    }
/*     mexPrintf("rcond = %f\n",rcond);
     mexPrintf("anorm = %f\n",anorm);*/
    *cond = 1.0/rcond;
    
    free(IWORK);
    free(WORK);
}

void SolveLinearSystem(int nrhs, mxArray *prhs[])
{
    double *A, *B;    /* pointers to input matrices */
    double *A2;  /* in/out arguments to DGESV*/
    size_t m,n,p;     /* matrix dimensions */
    mwSignedIndex *iPivot;   /* inputs to DGESV */
    mxArray *Awork, *mxPivot;
    mwSignedIndex info, dims[2];
    
    /* Check for proper number of arguments. */
    if ( nrhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:matrixDivide:rhs",
                "This function requires 2 input matrices.");
    }
    
    A = mxGetPr(prhs[0]); /* pointer to first input matrix */
    B = mxGetPr(prhs[1]); /* pointer to second input matrix */
    /* dimensions of input matrices */
    m = mxGetM(prhs[0]);
    p = mxGetN(prhs[0]);
    n = mxGetN(prhs[1]);
    
    /* Validate input arguments */
    if (p != mxGetM(prhs[1])) {
        mexErrMsgIdAndTxt("MATLAB:matrixDivide:matchdims",
                "Inner dimensions of matrices do not match.");
    }
    if (p != m) {
        mexErrMsgIdAndTxt("MATLAB:matrixDivide:square",
                "LAPACK function requires input matrix 1 must be square.");
    }
    
    /* DGESV works in-place, so we copy the inputs first. */
    Awork = mxCreateDoubleMatrix(m, p, mxREAL);
    A2 = mxGetPr(Awork);
    memcpy(A2, A, m*p*mxGetElementSize(prhs[0]));
    
    /* Create inputs for DGESV */
    dims[0] = m;
    dims[1] = p;
    mxPivot = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
    iPivot = (mwSignedIndex*)mxGetData(mxPivot);
    
    /* Call LAPACK */
    dgesv(&m,&n,A2,&m,iPivot,B,&p,&info);
    /* B now holds X */
    
    mxDestroyArray(Awork);
    mxDestroyArray(mxPivot);
}


void RegularizePhi(mxArray *mult[],ptrdiff_t *NN,const mxArray *phi_p,mxArray *phiR_p,mxArray *phiRP_p, mxArray *phiIso_p,double *r)
{
    mxArray *T_p = mult[2];
    mxArray *Matrices[2];
    double *T = mxGetPr(T_p);
    const double *phi = mxGetPr(phi_p);
    double *phiR = mxGetPr(phiR_p);
    double *phiRP = mxGetPr(phiRP_p);
    double *phiIso = mxGetPr(phiIso_p);
    ptrdiff_t i;
/*     mexPrintf("Regularizing with r = %f\n",*r);*/
/*     mexPrintf("phiR = \n\t\t");*/
    for(i=0;i<*NN;i++)
    {
        phiR[i] = (1.0-*r)*phi[i] + *r*phiIso[i];
/*         mexPrintf("%f\t",phiR[i]);*/
    }
/*     mexPrintf("\n\n");*/
    memcpy(phiRP, phiR, *NN*1*sizeof(double));
    Matrices[0] = T_p;
    Matrices[1]= phiRP_p;
    SolveLinearSystem(2,Matrices); /*Update phiRP = T\phiR;*/
/*     mexPrintf("phiRP = \n\t\t");
     for(i=0;i<*NN;i++)
     {
         mexPrintf("%f\t",phiRP[i]);
     }
     mexPrintf("\n\n");*/
    
}


void CopyMultipliers(mxArray *mult2[],mxArray *mult[],ptrdiff_t *nq,ptrdiff_t *NN)
{
    mxArray *alpha_p = mult[0] ;
    mxArray *beta_p = mult[1];
    mxArray *T_p = mult[2];
    mxArray *p_p = mult[3];
    mxArray *wG_p = mult[4];
    
    double *alpha = mxGetPr(alpha_p);
    double *beta = mxGetPr(beta_p);
    double *T = mxGetPr(T_p);
    double *p = mxGetPr(p_p);
    double *wG = mxGetPr(wG_p);
    
    mxArray *alpha2_p = mult2[0] ;
    mxArray *beta2_p = mult2[1];
    mxArray *T2_p = mult2[2];
    mxArray *p2_p = mult2[3];
    mxArray *wG2_p = mult2[4];
    
    double *alpha2 = mxGetPr(alpha2_p);
    double *beta2 = mxGetPr(beta2_p);
    double *T2 = mxGetPr(T2_p);
    double *p2 = mxGetPr(p2_p);
    double *wG2 = mxGetPr(wG2_p);
    
    memcpy(alpha2, alpha, *NN*1*sizeof(double));
    memcpy(beta2, beta, *NN*1*sizeof(double));
    memcpy(T2, T, (*NN)*(*NN)*sizeof(double));
    memcpy(p2, p, (*NN)*(*nq)*sizeof(double));
    memcpy(wG2, wG, (*nq)*sizeof(double));
    
}


void changeBasis(mxArray *mult[],mxArray *H_p,mxArray *phiRP_p,mxArray *g_p,ptrdiff_t *nq,ptrdiff_t *NN)
{
    mxArray *alpha_p = mult[0] ;
    mxArray *beta_p = mult[1];
    mxArray *T_p = mult[2];
    mxArray *p_p = mult[3];
    mxArray *wG_p = mult[4];
    mxArray *T_tmp_p = mxDuplicateArray(T_p);
    mxArray *beta_tmp_p = mxDuplicateArray(beta_p);
    double *p = mxGetPr(p_p);    
    double *phiRP = mxGetPr(phiRP_p);    
    double *g = mxGetPr(g_p);    
    double *H = mxGetPr(H_p);
    double *T = mxGetPr(T_p);
    double *beta = mxGetPr(beta_p);
    double *T_tmp = mxGetPr(T_tmp_p);
    double *beta_tmp = mxGetPr(beta_tmp_p);
    double one = 1.0, zero = 0.0;
    ptrdiff_t i,i2;
    mxArray *Matrices[2];
    mxArray *RHS_p = mxCreateDoubleMatrix(*NN,*nq+1+1,mxREAL);
    double *RHS = mxGetPr(RHS_p);
    for(i2=0;i2<*nq;i2++)
    {
        for(i=0;i<*NN;i++)
        {
            RHS[i2*(*NN)+i] = p[i2*(*NN)+i];
        }
    }
    for(i=0;i<*NN;i++)
    {    
        RHS[(*nq)*(*NN)+i] = phiRP[i];
        RHS[(*nq+1)*(*NN)+i] = g[i];
    }
    
    Matrices[0] = H_p;
    Matrices[1] = RHS_p;
    SolveLinearSystem(2,Matrices); /*Solve L\[p,phiRP,g] for optimal LR*/
    
    for(i2=0;i2<*nq;i2++)
    {
        for(i=0;i<*NN;i++)
        {
            p[i2*(*NN)+i] = RHS[i2*(*NN)+i];
        }
    }
    for(i=0;i<*NN;i++)
    {
        phiRP[i] =  RHS[(*nq)*(*NN)+i];
        g[i] = RHS[(*nq+1)*(*NN)+i];
    }
   
    
    dgemm("N", "N", NN, NN, NN, &one, T, NN, H, NN, &zero, T_tmp, NN); /*Multiplication T*L*/
    memcpy(T, T_tmp, *NN**NN*sizeof(double)); /*T = T_tmp = T*L;*/
    
    for(i=0;i<(*NN);i++)
    {
        for(i2=0;i2<(*NN);i2++)
        {
            T_tmp[i+i2*(*NN)] = T[i2+i*(*NN)]; /*T_tmp = T';*/
        }
    }
    
    for(i=0;i<(*NN);i++)
    {
        beta_tmp[i] = 0;
        for(i2=i;i2<(*NN);i2++)
        {
            beta_tmp[i]+=H[i2+i*(*NN)]*beta[i2];
        }
    }
    
    memcpy(beta, beta_tmp, (*NN)*1*sizeof(double)); /*beta = beta_tmp = L'*beta;*/
    
    mxDestroyArray(T_tmp_p);
    mxDestroyArray(beta_tmp_p);
    mxDestroyArray(RHS_p);
}



void backTrackLineSearch(mxArray *mult[],mxArray *phiRP_p,mxArray *g_p,double *f,ptrdiff_t *nq,ptrdiff_t *NN,double *xi,bool *status)
{
    mxArray *alpha_p = mult[0] ;
    mxArray *beta_p = mult[1];
    mxArray *T_p = mult[2];
    mxArray *p_p = mult[3];
    mxArray *wG_p = mult[4];
    mxArray *w0_p = mult[6];
    
    mxArray *Matrices[2];
    double *g = mxGetPr(g_p);
    double *beta = mxGetPr(beta_p);
    double *alpha = mxGetPr(alpha_p);
    double *wG = mxGetPr(wG_p);
    double *w0 = mxGetPr(w0_p);
    double *p = mxGetPr(p_p);
    double *T = mxGetPr(T_p);
    
    mxArray *beta_new_p = mxDuplicateArray(beta_p);
    mxArray *wG_new_p = mxDuplicateArray(wG_p);
    mxArray *T_tmp_p =  mxDuplicateArray(T_p);
    double *T_tmp = mxGetPr(T_tmp_p); /*T'*/
    
    double *beta_new = mxGetPr(beta_new_p);
    double *wG_new = mxGetPr(wG_new_p);
    
    ptrdiff_t i,i2,j;
    
    /*Problem parameters*/
    double delta = 1.0e-3;
    double chi = 0.5;
    double aeps = 0;
    double norm_tmp = 0;
    double eps = 0.0000000000000002220446049250313080847263336181640625; /*machine accuracy*/
    
    double f_new = 0;
    
    for(i=0;i<*NN;i++)
    {
        norm_tmp += g[i]*g[i];
        aeps += beta[i]*beta[i];
    }
    aeps = eps*sqrt(aeps)/sqrt(norm_tmp);
    
    *xi = 1;
/*     mexPrintf("xi = %f\n",*xi);*/
    while(*xi>aeps)
    {
        double fTest = 0;
/*         mexPrintf("xi = %f\n",*xi);*/
        
        for(i=0;i<(*NN);i++) /*Update along search direction d = -g*/
        {
            beta_new[i] = beta[i] - *xi*g[i];
/*             mexPrintf("beta_new[%d] = %f\n",i,beta_new[i]);*/
        }
        for(j=0;j<*nq;j++)
        {
            wG_new[j] = 0;
            for(i=0;i<(*NN);i++)
            {
                wG_new[j]+=p[j*(*NN)+i]*beta_new[i];
                
            }
            wG_new[j] = w0[j]*__ieee754_exp(wG_new[j]);
/*             mexPrintf("wG_new[%d] = %f\n",j,wG_new[j]);*/
        }
        
        
        f_new = fobjP(wG_new_p,beta_new_p,phiRP_p,nq,NN);
/*         mexPrintf("f_new = %3.16f\n",f_new);*/
        
        for(i=0;i<(*NN);i++)
        {
            fTest+=-g[i]*g[i];
        }
        fTest = *f + delta*(*xi)*fTest;
/*         mexPrintf("fTest = %3.16f\n",fTest);*/
/*         mexPrintf("(1.0- f_new/fabs(f_new)*eps)*f_new = %E\n",(1.0- f_new/fabs(f_new)*eps)*f_new);*/
/*         mexPrintf("((1.0+ fTest/fabs(fTest)*eps)*fTest) = %E\n",(1.0+ fTest/fabs(fTest)*eps)*fTest);*/
/*         mexPrintf("Crit Difference = %E\n",(1.0- f_new/fabs(f_new)*eps)*f_new-(1.0+ fTest/fabs(fTest)*eps)*fTest);*/
        if( (1.0- signum(f_new)*eps)*f_new <= (1.0+ signum(fTest)*eps)*fTest)
        {
/*             mexPrintf("Criterion fulfilled\n");*/
            /*Overwrite beta and wG in mult*/
           memcpy(beta, beta_new, *NN*1*sizeof(double)); /*beta = beta_new;*/
           memcpy(wG, wG_new, *nq*1*sizeof(double)); /*wG = wG_new;*/
            
           gradP(mult,phiRP_p,nq,NN,g_p);
           *f = f_new;
           break;
        }
        else
        {
            *xi = chi*(*xi);
        }
    }
/*             mexPrintf("xi = %f\n",*xi);*/
    
    if(*xi<=aeps) /*Linesearch failed*/
    {
        *status = false;
    }
    else
    {
        
        for(i=0;i<(*NN);i++)
        {
            for(i2=0;i2<(*NN);i2++)
            {
                T_tmp[i+i2*(*NN)] = T[i2+i*(*NN)]; /*T_tmp = T';*/
            }
        }
/*         mexPrintf("T =\n\t\t");
         for(i=0;i< *NN;i++)
         {
             for(i2=0;i2< *NN;i2++)
             {
                 mexPrintf("%f\t",T[i2*(*NN)+i]);
             }
             mexPrintf("\n\t\t");
         }
         mexPrintf("\n\n");
         mexPrintf("T' =\n\t\t");
         for(i=0;i< *NN;i++)
         {
             for(i2=0;i2< *NN;i2++)
             {
                 mexPrintf("%f\t",T_tmp[i2*(*NN)+i]);
             }
             mexPrintf("\n\t\t");
         }
         mexPrintf("\n\n");
         mexPrintf("alphav = \n\t\t");
         for(i=0;i<*NN;i++)
         {
             mexPrintf("%f\t",alpha[i]);
         }
         mexPrintf("\n\n");*/
        
        Matrices[0] = T_tmp_p;
        Matrices[1] = alpha_p;
        memcpy(alpha, beta, (*NN)*1*sizeof(double)); /*beta = beta_tmp = L'*beta;*/
        SolveLinearSystem(2,Matrices); /*Calculate alpha_0*/
/*         mexPrintf("alphan = \n\t\t");*/
/*         for(i=0;i<*NN;i++)*/
/*         {*/
/*             mexPrintf("%f\t",alpha[i]);*/
/*         }*/
/*         mexPrintf("\n\n");*/
    }
    
    mxDestroyArray(beta_new_p);
    mxDestroyArray(wG_new_p);
    mxDestroyArray(T_tmp_p);
    
/*     mexPrintf("f = %f\n",*f);*/
}



void NewtonLoop(mxArray *plhs[],mxArray *mult[],mxArray *MinimizationVals[],ptrdiff_t *k,ptrdiff_t *itermax, double *tol, ptrdiff_t *nq, ptrdiff_t *NN,mxLogical *status,mxArray *normg_p,double *maxHcond)
{
 
    /*Inputs and Outputs*/
    mxArray *alpha_p = mult[0] ;
    mxArray *beta_p = mult[1];
    mxArray *T_p = mult[2];
    mxArray *p_p = mult[3];
    mxArray *wG_p = mult[4];
    mxArray *p0_p = mult[5];
    mxArray *w0_p = mult[6];
    
    double *alpha = mxGetPr(alpha_p);
    double *beta = mxGetPr(beta_p);
    double *T = mxGetPr(T_p);
    double *p = mxGetPr(p_p);
    double *wG = mxGetPr(wG_p);
    double *p0 = mxGetPr(p0_p);
    double *w0 = mxGetPr(w0_p);
    
    mxArray *f_p = MinimizationVals[0];
    mxArray *g_p = MinimizationVals[1];
    mxArray *H_p = MinimizationVals[2];
    mxArray *rList2_p = MinimizationVals[3];
    mxArray *phiIso_p = MinimizationVals[4];
    mxArray *phiR_p = MinimizationVals[5];
    mxArray *phiRP_p = MinimizationVals[6];
    double *f = mxGetPr(f_p);
    double *g = mxGetPr(g_p);
    double *H = mxGetPr(H_p);
    double *rList2 = mxGetPr(rList2_p);
    double *phiIso = mxGetPr(phiIso_p);
    double *phiR = mxGetPr(phiR_p);
    double *phiRP = mxGetPr(phiRP_p);
    
    double *normg = mxGetPr(normg_p);
    
    /*Helpers*/
    double r = 0;
    const double xitol = .9;
    ptrdiff_t r2idx = -1;
    ptrdiff_t xi0Counter = 0;
    double xiPrev = 1;
    double xi=1;
    ptrdiff_t zeroStepLimit = (ptrdiff_t) 1e16;
    
    double minNormg = 1.0e16;
    double *distribution=(double *) malloc((*nq)*sizeof(double));
    double normGrad = 0;
    double cond = 0;
    mxArray *alphaPrev_p = mxDuplicateArray(alpha_p);
    mxArray *betaPrev_p = mxDuplicateArray(beta_p);
    mxArray *TPrev_p = mxDuplicateArray(T_p);
    mxArray *pPrev_p = mxDuplicateArray(p_p);
    mxArray *wGPrev_p = mxDuplicateArray(wG_p);
    mxArray *alphaPrev2_p = mxDuplicateArray(alpha_p);
    mxArray *betaPrev2_p = mxDuplicateArray(beta_p);
    mxArray *TPrev2_p = mxDuplicateArray(T_p);
    mxArray *pPrev2_p = mxDuplicateArray(p_p);
    mxArray *wGPrev2_p = mxDuplicateArray(wG_p);
    mxArray *alphaMin_p = mxDuplicateArray(alpha_p);
    mxArray *betaMin_p = mxDuplicateArray(beta_p);
    mxArray *TMin_p = mxDuplicateArray(T_p);
    mxArray *pMin_p = mxDuplicateArray(p_p);
    mxArray *wGMin_p = mxDuplicateArray(wG_p);
    mxArray *multPrev[] = {alphaPrev_p,betaPrev_p,TPrev_p,pPrev_p,wGPrev_p,p0_p,w0_p};
    mxArray *multPrev2[] = {alphaPrev2_p,betaPrev2_p,TPrev2_p,pPrev2_p,wGPrev2_p,p0_p,w0_p};
    mxArray *multMin[] = {alphaMin_p,betaMin_p,TMin_p,pMin_p,wGMin_p,p0_p,w0_p};
    ptrdiff_t i,j; /*Indices*/
    ptrdiff_t cholStatus;
    bool BadHessian = false;
    bool LineSearchStatus = true;
    
/*     mexPrintf("Starting NewtonLoop\n");*/
    /*Function*/
    *status = false;
    
    while(*k<*itermax && *status==false)
    {
        *k+=1;
/*         mexPrintf("Iteration = %d\n",*k);*/
        /*--------Objective function-----------*/
        for(j=0;j<*nq;j++)
        {
            distribution[j] = 0;
            for(i=0;i<*NN;i++)
            {
                distribution[j]+=p0[j*(*NN)+i]*alpha[i];
                
            }
            distribution[j] = w0[j]*__ieee754_exp(distribution[j]);
        }
        normg[*k] = 0;
        for(i=0;i<*NN;i++)
        {
            double norm_tmp = 0;
            for(j=0;j<*nq;j++)
            {
                norm_tmp +=p0[j*(*NN)+i]*distribution[j];
            }
            normg[*k] += (norm_tmp-phiR[i])*(norm_tmp-phiR[i]);
        }
        normg[*k] = sqrt(normg[*k]);
        /*--------Objective function-----------*/
        
/*         mexPrintf("norm[%d] = %2.16f\n",*k,normg[*k]);*/
/*         mexPrintf("tol-norm = %2.16f\n",*tol);*/
        /*Stopping criterion*/
        if(normg[*k]<*tol)
        {
/*             mexPrintf("Succeeded\n");*/
            *status = true;
            break;
        }
        else
        {
            CopyMultipliers(multPrev2,multPrev,nq,NN);
            CopyMultipliers(multPrev,mult,nq,NN);
            
/*             PrintMultipliers(mult,nq,NN);*/
            hessP(mult,nq,NN,H_p,&cond,&cholStatus);
            
            if(cond>*maxHcond)
                *maxHcond = cond;
            
            BadHessian = (cholStatus==0 || cond>1.0e16);
                        
            if(BadHessian)
            {
                if((r2idx>=0 || *k==0))
                {
                    /*Give up*/
                    break;
                }
            }
            else
            {
                changeBasis(mult,H_p,phiRP_p,g_p,nq,NN);
/*                 mexPrintf("Mults after changeBasis\n\n");*/
/*                 PrintMultipliers(mult,nq,NN);*/
                backTrackLineSearch(mult,phiRP_p,g_p,f,nq,NN,&xi,&LineSearchStatus);
/*                 mexPrintf("backTrack gives xi = %f\n\n",xi);*/
/*                 mexPrintf("Mults after backTrack\n\n");*/
/*                 PrintMultipliers(mult,nq,NN);*/
/*                 for(i=0;i<*NN;i++)*/
/*                 {*/
/*                     mexPrintf("g[%d] = %E\n",i,g[i]);*/
/*                 }*/
                
                if(normg[*k]<minNormg)
                {
                    minNormg = normg[*k];
                    CopyMultipliers(multMin,multPrev,nq,NN);
                }
            }
            
            
/*             if((BadHessian) && *k>0 && r2idx == 0)*/
            if((BadHessian||xi<xitol) && *k>1 && r2idx == -1)
            {
                CopyMultipliers(mult,multPrev2,nq,NN);
/*                 mexPrintf("r2idx = %d\n",r2idx);*/
/*                 mexPrintf("k = %d\n",*k);*/
                r2idx += 1;
                r = rList2[r2idx];
/*                 mexPrintf("XI: Regularizing with r = %E\n",r);*/
/*                 mexPrintf("xi = %f\n",xi);*/
                RegularizePhi(mult,NN,phiR_p,phiRP_p,phiRP_p,phiIso_p,&r); /*We only want to change phiRP*/
                
/*                 PrintMultipliers(mult,nq,NN);*/
                *f = fobjP(wG_p,beta_p,phiRP_p,nq,NN);
                gradP(mult,phiRP_p,nq,NN,g_p);
/*                 mexPrintf("f = %E\n",*f);*/
                
                
/*                 for(i=0;i<*NN;i++)
                 {
                     mexPrintf("g[%d] = %E\n",i,g[i]);
                 }
                 for(i=0;i<*NN;i++)
                 {
                     mexPrintf("phiRP[%d] = %f\n",i,phiRP[i]);
                 }*/
                
            }
            else
            {
/*                 mexPrintf("r = %f\n",r);*/
/*                 mexPrintf("xi = %f\n",xi);*/
                if(xi>=xitol && r>0)
                {
                    r2idx += 1;
                    r = rList2[r2idx];
/*                     mexPrintf("Regularizing with r = %f\n",r);*/
                    RegularizePhi(mult,NN,phiR_p,phiRP_p,phiRP_p,phiIso_p,&r); /*We only want to change phiRP*/
                    
                    *f = fobjP(wG_p,beta_p,phiRP_p,nq,NN);
                    gradP(mult,phiRP_p,nq,NN,g_p);   
                    

/*                     for(i=0;i<*NN;i++)
                     {
                         mexPrintf("phiRP[%d] = %f\n",i,phiRP[i]);
                     }*/
                    
                }   
            }
            
            
            if(xi<=0.0 && xiPrev<=0)
            {
                xi0Counter+=1;
                if(xi0Counter>zeroStepLimit)
                    break;
                
            }
            xiPrev = xi;
            
            
        }
    }
    
    
    if(!*status)
    {
/*         mexPrintf("No Status!!");*/
        CopyMultipliers(mult,multMin,nq,NN);
    }
    
    free(distribution);
    mxDestroyArray(alphaPrev_p);
    mxDestroyArray(betaPrev_p);
    mxDestroyArray(TPrev_p);
    mxDestroyArray(pPrev_p);
    mxDestroyArray(alphaPrev2_p);
    mxDestroyArray(betaPrev2_p);
    mxDestroyArray(TPrev2_p);
    mxDestroyArray(pPrev2_p);
    mxDestroyArray(alphaMin_p);
    mxDestroyArray(betaMin_p);
    mxDestroyArray(TMin_p);
    mxDestroyArray(pMin_p);
    
    
    
}





void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mxArray *phi_p, *p_p,*alpha_p, *beta_p, *T_p, *iter_p,*normg_p, *status_p, *r_out_p, *el1_out_p;
    mxArray *T_tmp_p, *wG_p, *H_p, *g_p, *phiR_p, *phiRP_p, *f_p, *maxHCond_p;
    mxArray *rList_p, *rList2_p, *p0_p, *w0_p,*k0_p, *phiIso_p, *tolin_p;
    double *phi, *itermax, *rList, *rList2, *phiIso, *p, *w0, *wG, *N, *alpha, *beta, *T, *iter, *k0, *tolin, *H, *g,*phiR, *phiRP, *p0;
    mxLogical status=0;
    
    ptrdiff_t nq, length_rList;
    ptrdiff_t NN;
    ptrdiff_t i,i2,j, el1;
    
    mxArray *Matrices[2];
    double tol;
    double psi0;
    mxArray *NewtonOuts[3], *mult[7],*MinimizationVals[7];
    double *T_tmp;
    double cond=0;
    double *f;
    ptrdiff_t k = -1;
    ptrdiff_t stat=0;
    ptrdiff_t maxitLocal;
    double r=0;
    ptrdiff_t isoSwitchState = 0;
    if ( nrhs != 12) {
        mexErrMsgIdAndTxt("MATLAB:dualAdaptivePoly_mex:rhs",
                "This function requires 12 inputs.");
    }
    
    phi_p = mxDuplicateArray(prhs[0]);
    phi = mxGetPr(phi_p); /*NNxn*/
    
    itermax  = mxGetPr(prhs[3]);
    rList_p = mxDuplicateArray(prhs[4]);
    rList  = mxGetPr(rList_p);
    rList2_p = mxDuplicateArray(prhs[5]);
    rList2  = mxGetPr(rList2_p);
    phiIso_p = mxDuplicateArray(prhs[6]);
    phiIso  = mxGetPr(phiIso_p);
    p_p = mxDuplicateArray(prhs[7]);
    p  = mxGetPr(p_p); /*angular basis NNxnq*/
    p0_p = mxDuplicateArray(prhs[7]);
    p0  = mxGetPr(p0_p); /*angular basis NNxnq*/
    w0_p = mxDuplicateArray(prhs[8]);
    w0 = mxGetPr(w0_p);
    k0_p = mxDuplicateArray(prhs[9]);
    k0 = mxGetPr(k0_p);
    tolin_p = mxDuplicateArray(prhs[10]);
    tolin = mxGetPr(tolin_p);
    tol = tolin[0];
    N  = mxGetPr(prhs[11]); 

    alpha_p = plhs[0] = mxDuplicateArray(prhs[1]); /*NNx1*/
    alpha = mxGetPr(alpha_p);
    beta_p = plhs[1] = mxDuplicateArray(prhs[1]); /*NNx1*/
    beta = mxGetPr(beta_p);
    T_p = plhs[2] = mxDuplicateArray(prhs[2]); /*NNxNN*/
    T = mxGetPr(T_p);
    status_p = plhs[3] = mxCreateLogicalScalar(status);
    iter_p = plhs[4] = mxCreateDoubleScalar(0.0);
    iter = mxGetPr(iter_p);
    normg_p = plhs[5] = mxCreateDoubleMatrix((mwSize) (itermax[0]+0.5),1,mxREAL);
    
    
    
/*figure out dimensions*/
    NN = (size_t)mxGetM(prhs[0]);/*y = number of moments, x = number of systems to solve*/
    if(mxGetM(rList_p)>mxGetN(rList_p))
        length_rList = (ptrdiff_t) mxGetM(rList_p);
    else
        length_rList = (ptrdiff_t) mxGetN(rList_p);
    
    nq = (ptrdiff_t) mxGetN(p_p);
    
    
    /*Generate helpers*/
    wG_p = mxCreateDoubleMatrix(nq,1,mxREAL);
    wG = mxGetPr(wG_p);
    g_p = mxCreateDoubleMatrix(NN,1,mxREAL);
    g = mxGetPr(g_p);
    H_p = mxCreateDoubleMatrix(NN,NN,mxREAL);
    H = mxGetPr(H_p);
    phiR_p = mxCreateDoubleMatrix(NN,1,mxREAL);
    phiR = mxGetPr(phiR_p);
    phiRP_p = mxCreateDoubleMatrix(NN,1,mxREAL);
    phiRP = mxGetPr(phiRP_p);
    
    T_tmp_p =  mxCreateDoubleMatrix(NN,NN,mxREAL);
    T_tmp = mxGetPr(T_tmp_p); /*T0'*/
    for(i=0;i<NN;i++)
    {
        for(i2=0;i2<NN;i2++)
        {
            T_tmp[i+i2*NN] = T[i2+i*NN]; /*T_tmp = T';*/
        }
    }
    
    if(phi[0]<=0)
    {
        mexErrMsgIdAndTxt("MATLAB:dualAdaptivePoly_mex:NegativeDensity",
                "Unrealizable moments: psi[0] is negative.");
    }
    
    
    /*Normalize psi*/
    psi0 = phi[0];
    for(i=NN-1;i>=0;i--)
        phi[i] = phi[i]/psi0;
    
/*     mexPrintf("phi = \n\t\t");
     for(i=0;i<NN;i++)
     {
        mexPrintf("%f\t",phi[i]);
     }
     mexPrintf("\n\n");
     mexPrintf("phiISO = \n\t\t");
     for(i=0;i<NN;i++)
     {
         mexPrintf("%f\t",phiIso[i]);
     }
     mexPrintf("\n\n");*/
    tol = tol/psi0;
    
    /*Calculate multipliers*/
    
    Matrices[0] = T_tmp_p;
    Matrices[1] = alpha_p;
    SolveLinearSystem(2,Matrices); /*Calculate alpha_0*/
    Matrices[0] = T_p;
    Matrices[1] = p_p; /*p_p = p0_p*/
    SolveLinearSystem(2,Matrices); /*Calculate pP*/
    
    for(j=0;j<nq;j++)
    {
        
        wG[j] = 0;
        for(i=0;i<NN;i++)
        {
            wG[j]+=p[j*NN+i]*beta[i];
            
        }
        wG[j] = w0[j]*__ieee754_exp(wG[j]);
    }
    mult[0] = alpha_p;
    mult[1] = beta_p;
    mult[2] = T_p;
    mult[3] = p_p;
    mult[4] = wG_p;
    mult[5]= p0_p;
    mult[6] = w0_p;
    
    r = rList[0];
    RegularizePhi(mult,&NN,phi_p,phiR_p,phiRP_p,phiIso_p,&r);
            
    Matrices[1] = phiR_p; 
    SolveLinearSystem(2,Matrices); /*Calculate phiP*/
    
/*     PrintMultipliers(mult,&nq,&NN);*/
    f_p = mxCreateDoubleScalar(fobjP(wG_p,beta_p,phiRP_p,&nq,&NN));
    f = mxGetPr(f_p);
/*     mexPrintf("fobj = %f\n",*f);*/
    gradP(mult,phiRP_p,&nq,&NN,g_p);
    hessP(mult,&nq,&NN,H_p,&cond,&stat);
    
    /*Check if initial conditions are okay*/
    if(stat==0 || cond>1.0e16 || fabs(*f)>1.0e16)
    {
        isoSwitchState = 1;
    }
    else
    {
       for(i2=0;i2<nq;i2++)
       {
           if(fabs(wG[i2])>1.0e16)
           {
               isoSwitchState=1;
               break;
           }
       }
       for(i=0;i<NN;i++)
       {
           if(fabs(g[i])>1.0e16)
           {
               isoSwitchState=1;
               break;
           }
       } 
    }
    
/*     for(i=0;i<NN;i++)
     {
         mexPrintf("g[%d] = %f\n",i,g[i]);
     }
     for(i=0;i<NN;i++)
     {
         mexPrintf("phiRP[%d] = %f\n",i,phiRP[i]);
     }
     mexPrintf("Iso = %d\n",isoSwitchState);*/
    
    
    MinimizationVals[0] = f_p;
    MinimizationVals[1] = g_p;
    MinimizationVals[2] = H_p;
    MinimizationVals[3] = rList2_p;
    MinimizationVals[4] = phiIso_p;
    MinimizationVals[5] = phiR_p;
    MinimizationVals[6] = phiRP_p;
    
    status = 0; /*Set initial status to false*/
    for(i=isoSwitchState;i<2;i++)
    {
/*         mexPrintf("i = %d\n",i);*/
        if(i==1)
        {
            /*Isotropic multipliers*/
            alpha[0] = log(0.5);
            for(i=1;i<NN;i++)
                alpha[i]=0;
            for(i=0;i<(NN)*(NN);i++)
                T[i] = 0;
            for(i=0;i<NN;i++)
                T[i*(NN)+i] = 1;
            memcpy(beta, alpha, NN*1*sizeof(double));
            memcpy(p, p0, NN*nq*sizeof(double));
            for(i=0;i<nq;i++)
                wG[i] = w0[i]*__ieee754_exp(alpha[0]*p[i]);
            
/*              PrintMultipliers(mult,&nq,&NN);*/
            
            /*Recalculate initial values for f,g,H*/
            *f = fobjP(wG_p,beta_p,phiRP_p,&nq,&NN);
            gradP(mult,phiRP_p,&nq,&NN,g_p);
            hessP(mult,&nq,&NN,H_p,&cond,&stat);
        }
/*         r = rList[0];*/
/*         mexPrintf("Diff = %f\n",rList[0]-rList[1]);*/
/*         RegularizePhi(mult,&NN,phi_p,phiR_p,phiRP_p,phiIso_p,&r);*/
        
        for(el1 = 0;el1<length_rList;el1++) /*Regularization       */         
        {
            r = rList[el1];
/*             mexPrintf("r = %f\n",r);*/
/*             mexPrintf("el1 = %d\n",el1);*/
/*             mexPrintf("length_rList = %d\n",length_rList);*/
            RegularizePhi(mult,&NN,phi_p,phiR_p,phiRP_p,phiIso_p,&r);
            *f = fobjP(wG_p,beta_p,phiRP_p,&nq,&NN);
            gradP(mult,phiRP_p,&nq,&NN,g_p);
            hessP(mult,&nq,&NN,H_p,&cond,&stat);
            if(k+*k0 < (ptrdiff_t) (itermax[0]+0.5))
                maxitLocal = k+*k0;
            else
                maxitLocal = (ptrdiff_t) (itermax[0]+0.5);
                
            NewtonLoop(NewtonOuts,mult,MinimizationVals, &k,&maxitLocal, &tol, &nq, &NN, &status,normg_p, &cond);
/*             mexPrintf("k = %d\n",k);*/
            if(status || maxitLocal >= (ptrdiff_t) (itermax[0]+0.5))
            {
                break;
            }
            
        }
        if(status || maxitLocal >= (ptrdiff_t) (itermax[0]+0.5))
        {
            break;
        }

    }
    /*mxSetM(normg_p,k+1); /*Cut off the zeros at the end*/
    *iter = (double) k+1; /*matlab indexing*/
    *mxGetLogicals(status_p) = status;
    
    mxDestroyArray(phi_p);
    mxDestroyArray(p_p);
    mxDestroyArray(g_p);
    mxDestroyArray(H_p);
    mxDestroyArray(phiR_p);
    mxDestroyArray(phiRP_p);
    mxDestroyArray(T_tmp_p);
    mxDestroyArray(f_p);
    mxDestroyArray(rList_p);
    mxDestroyArray(rList2_p);
    mxDestroyArray(phiIso_p);
    mxDestroyArray(p0_p);
    mxDestroyArray(w0_p);
    mxDestroyArray(k0_p);
    mxDestroyArray(tolin_p);
    
    r_out_p = plhs[6] = mxCreateDoubleScalar(r);
    el1_out_p =  plhs[7] = mxCreateDoubleScalar(el1+1.0);/*matlab indexing*/
    maxHCond_p =  plhs[8] = mxCreateDoubleScalar(cond);
    plhs[9] = wG_p;
    
}
