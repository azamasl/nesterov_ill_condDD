
//  Copyright © 2016 Fatemeh Asl. All rights reserved.

#include "lbfgsb.h"
#include <iomanip>
#include <iostream>     // std::cout
#include <qd/dd_real.h>
#include <qd/fpu.h>
#include <algorithm>

double random(double lb, double ub);
int softmaxabs(integer *n, dd_real *mu, dd_real *x, dd_real *f, dd_real *g);
int yurileshoushessmoothed(integer *n, dd_real *mu, dd_real *x, dd_real *f, dd_real *g);
void writeoutput(dd_real *fLB, integer *iterLB, dd_real *muvals, int muvalsize);
int lbfgsbcreatemat(double *fLB, integer *iterLB, double *muvals, int muvalsize);

