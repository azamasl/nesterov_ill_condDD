//  Created by Azam Asl on 12/8/16.
//  Copyright © 2016 Azam Asl. All rights reserved.

#include <iostream>
#include "lbfgsb.h"
#include "yurileshouchessm.h"
using std::cout;
using std::endl;

int main(int argc, const char * argv[]) {

    // ensure that 80-bit arithmetic is not in place
    // this call forces 64-bit arithmetic
    unsigned int old_cw;
    fpu_fix_start(&old_cw);
    cout.precision(32);

    //std::cout.precision(16);  //setting the decimal precision globally
    static dd_real startper,normtol;
    static dd_real xx[1024], x[1024],f, fLB[1024], g[1024], mu, muvalues[1024];
    static dd_real l[1024],u[1024]; /*l   specifies the lower bounds,  u   specifies the upper bounds*/
    static dd_real wa[43251],dsave[29],pgtol, factr;
    
    static double fLBdouble[1024], muvaluesdouble[1024];
    
    static integer i, m, n,iprint,iterLB[1024];
    static integer nbd[1024], iwa[3072];
    static integer taskValue;
    static integer *task=&taskValue; /* must initialize !! */
    static integer csaveValue;
    static integer *csave=&csaveValue;
    static integer isave[44];
    static logical lsave[4];
    
//------------------------Initialization--------------------------
    iprint = 0;
    factr = dd_real("1.0");//exteremly high precision
    pgtol = dd_real("1e-100");// I set this to be the same as options.normtol = 1e-100 in Michael's code.
    normtol = dd_real("1e-100");


    int muvaluessize = 20;
    //int maxit = 100000;
    m = 5;
    n = 4;
    startper = dd_real("0.01");

    /* Azam: Michael is using normal distribution (randn) for X0.
     * I am using uniform distribution. But I think it's not a big deal*/
    for (i = 0; i < n; ++i) {
        xx[i] = dd_real(1.0+startper*random(-1.0,1.0));
    }

    for (i = 0; i < n; i++) {//This is a an unbounded problem
        nbd[i] = 0;
    }
    
    for(int j=0; j<muvaluessize; ++j){
        factr = dd_real("0.0");
        pgtol = dd_real("0.0");
        mu = pow(dd_real("10.0"), -(j+1));
        muvalues[j] = mu;
        for (i = 0; i < n; ++i) {
            x[i]=xx[i];
            //cout << "x[" << i<< "] = " << x[i] <<endl;
        }

        cout << "\n\n\n AZAM______________________calling L-BFGS-B-" << m << " with mu = " << mu <<"_______________________"<< endl;

        *task = START;
    L111:
        setulb(&n, &m, x, l, u, nbd, &f, g, &factr, &pgtol, wa, iwa, task, &iprint, csave, lsave, isave, dsave);
        //evaluate f and g in x:
        if ( IS_FG(*task) ) {
            /* Trying to impose the criteria as in Michael's matlab code however it doesn't work

            if (isave[9] >= maxit) {//Stop when the number of the iterations exceeds the maxit.
                *task = STOP_MAXITER;
                cout << "STOP MAX ITERATONS LIMIT REACHED"<< endl;
            }
            */
            yurileshoushessmoothed(&n, &mu, x, &f, g);
            goto L111;
        }
        
        if ( *task== NEW_X ) {

//            if (isave[29] >= maxit) {//Stop when the number of the iterations exceeds the maxit.
//                *task = STOP_MAXITER;
//                cout << "STOP MAX ITERATONS LIMIT REACHED"<< endl;
//            }
            
            
//            /* 1) Terminate if the total number of f and g evaluations */
            if (isave[33] >= 500000) {
                *task = STOP_ITER;
            }
            
            goto L111;
        }
        
        
        fLB[j] = f;//saving the final f for the current mu for the plot.
        iterLB[j] = isave[29];//saving the final iteration numbers for the current mu for the plot.
        /* The found solution: */
        for (i = 0; i < n; ++i) {


        }
        

    }
    
    //lbfgsbplot(fLB, muvalues, muvaluessize);

    //writeoutput(fLB, iterLB, muvalues, muvaluessize);

    for(int j=0; j<muvaluessize; ++j){
        fLBdouble[j] = to_double(fLB[j]);
        muvaluesdouble[j] = to_double(muvalues[j]);
    }

    lbfgsbcreatemat(fLBdouble, iterLB, muvaluesdouble, muvaluessize);

    fpu_fix_end(&old_cw);
    return 0;
    
} /* MAIN__ */




