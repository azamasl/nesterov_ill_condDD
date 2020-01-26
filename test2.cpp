//
//// Created by Fatemeh Asl on 5/25/17.
////  Created by Azam Asl on 12/8/16.
//  Copyright © 2016 Azam Asl. All rights reserved.

#include <iostream>
#include "lbfgsb.h"
#include "yurileshouchessm.h"
#include "T10_dd.h"
using std::cout;
using std::endl;

int main(int argc, const char * argv[]) {

    // ensure that 80-bit arithmetic is not in place
    // this call forces 64-bit arithmetic
    unsigned int old_cw;
    fpu_fix_start(&old_cw);
    cout.precision(32);

    static dd_real l[1024],u[1024]; /*l   specifies the lower bounds,  u   specifies the upper bounds*/
    static dd_real wa[43251],dsave[29],pgtol, factr;
    static integer nbd[1024], iwa[3072];
    static integer taskValue,iprint;
    static integer *task=&taskValue; /* must initialize !! */
    static integer csaveValue;
    static integer *csave=&csaveValue;
    static integer isave[44];
    static logical lsave[4];

    iprint = 101;
    //TODO: double check this!
    factr = dd_real("1.0");//exetremely high precision
    pgtol = dd_real("1e-100");// I set this to be the same as options.normtol = 1e-100 in Michael's code.


    static integer i, m, n;
    static dd_real x0[1024];
    static dd_real xf[1024];
    static dd_real f;
    static dd_real g[1024];

    int maxit = 10e8;		//maximum iterations
    m = 100;               //Setting to a large number to simulate full BFGS
    n = 10;

    for (i = 0; i < n; i++) {nbd[i] = 0;} //This is a an unbounded problem

    for (i = 0; i < n; ++i) {
        x0[i] = dd_real(random(-1.0,1.0));
        cout << "x[" << i<< "] = " << x0[i] <<endl;
    }
    for (int j=0; j<n; j++) {xf[j] = x0[j];}


    *task = START;
    L111:
    setulb(&n, &m, x0, l, u, nbd, &f, g, &factr, &pgtol, wa, iwa, task, &iprint, csave, lsave, isave, dsave);

    if ( IS_FG(*task) ) {       //evaluate f and g in x:
        test29f02_dd(&f, g, x0, n);
        goto L111;
    }

    if ( *task== NEW_X ) {
        if (isave[29] >= maxit) {
            *task = STOP_ITER;
        }
        goto L111;
    }


    fpu_fix_end(&old_cw);
    return 0;

} /* MAIN__ */





