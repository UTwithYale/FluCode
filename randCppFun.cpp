/*=================================================================
 *
 * randBin.cpp	Sample .MEX file corresponding to YPRIME.M
 *	        Solves simple 3 body orbit problem
 *
 * The calling syntax is:
 *
 *		[Np] = randBin(N p)
 *
 *  You may also want to look at the corresponding M-code, yprime.m.
 *
 * This is a MEX-file for MATLAB.
 * Copyright 1984-2017 The MathWorks, Inc.
 *
 *=================================================================*/
#include <math.h>
#include "mex.h"
#include <random>


/* Input Arguments */

#define C_IN prhs[0]
#define N_IN prhs[1]
#define P_IN prhs[2]


/* Output Arguments */

#define YP_OUT plhs[0]

#if !defined(MAX)
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#endif

std::mt19937_64 generator;

static void binom(int n1, int n2, double yp[],  double N[], double* p) {
    int i ;
    if( n2 == 1 ) {
        for( i=0; i<n1; i++) {
            std::binomial_distribution<int> distribution( (unsigned long)N[i], p[0]);
            yp[i] = distribution(generator) ;
        }
    } else if( n1 == 1 ) {
        for( i=0; i<n2; i++) {
            std::binomial_distribution<int> distribution( N[0], p[i]);
            yp[i] = distribution(generator) ;
        }        
    } else {
        for( i=0; i<n1; i++) {
            std::binomial_distribution<int> distribution( N[i], p[i]);
            yp[i] = distribution(generator) ;
        }
    }
    return;
}

static void poisson( int n, double yp[], double mu[] )
{
    int i ;
    for( i=0; i<n; i++) {
        std::poisson_distribution<int> distribution( mu[i]);
        yp[i] = distribution(generator) ;
    }
}

void init_rand( unsigned long seed ) {
    generator.seed( seed ) ;
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])

{
#if MX_HAS_INTERLEAVED_COMPLEX
    mxDouble* yp;
    mxDouble *c, *N, *p;
#else
    double* yp;
    double *c, *N, *p;
#endif
    size_t m, n, mP, nP, nC;

    /* Check for proper number of arguments */

    if (nrhs != 3) {
        mexErrMsgIdAndTxt("MATLAB:yprime:invalidNumInputs", "Three input arguments required.");
    } else if (nlhs > 1) {
        mexErrMsgIdAndTxt("MATLAB:yprime:maxlhs", "Too many output arguments.");
    }

    /* check to make sure the first input argument is a real matrix */
    if (!mxIsDouble(P_IN) || mxIsComplex(N_IN)) {
        mexErrMsgIdAndTxt("MATLAB:yprime:invalidT", "First input argument must be a real matrix.");
    }

    /* check to make sure the second input argument is a real matrix */
    if (!mxIsDouble(P_IN) || mxIsComplex(N_IN)) {
        mexErrMsgIdAndTxt("MATLAB:yprime:invalidY", "Second input argument must be a real matrix.");
    }

    /* Check the dimensions of Y.  Y can be 4 X 1 or 1 X 4. */

    nC = mxGetM( C_IN ) ;
        
    
    m = mxGetM(N_IN);
    n = mxGetN(N_IN);

    mP = mxGetM(P_IN);
    nP = mxGetN(P_IN);
    
    if( !( ((mP*nP==m*n) )  || ((mP*nP==1) ) || ((n*m==1)) ) ) {
        mexPrintf("P: %d %d N:%d %d\n",mP,nP,m,n);
        mexErrMsgIdAndTxt("MATLAB:yprime:invalidY", "N and p need equal dim or be scalar");
    }
    
    
    /* Create a matrix for the return argument */
    if( m*n > mP*nP)
        YP_OUT = mxCreateDoubleMatrix((mwSize)m,(mwSize)n, mxREAL);
    else
        YP_OUT = mxCreateDoubleMatrix((mwSize)mP,(mwSize)nP, mxREAL);
    
    /* Assign pointers to the various parameters */
#if MX_HAS_INTERLEAVED_COMPLEX
    yp = mxGetDoubles(YP_OUT);
    p = mxGetDoubles(P_IN);
    N = mxGetDoubles(N_IN);
    c = mxGetDoubles(C_IN);
#else
    yp = mxGetPr(YP_OUT);
    p = mxGetPr(P_IN);
    N = mxGetPr(N_IN);
    c = mxGetPr(C_IN);
#endif
        if( nC > 0) {
            if( abs( c[0]-1.0)<1e-10 ) {
                if( n > 0 ) {
                    init_rand( (unsigned long)N[0] ) ;
                }
            } else if( abs( c[0]-2.0)<1e-10 ) {
                poisson(n*m, yp, N ) ;
            } else {
               binom(n*m,nP*mP,yp, N, p);
            }
        } else {
            binom(n*m,nP*mP,yp, N, p);
        }
    return;
}

/* LocalWords:  yp maxlhs
 */
