//  Created by Azam Asl on 12/9/16.
//  Copyright © 2016 Azam	 Asl. All rights reserved.

#include "lbfgsb.h"
#include "yurileshouchessm.h"
#include <random>
#include <ctime>
#include <algorithm>
#include "mat.h"
#include <qd/qd_real.h>
#include <qd/fpu.h>
#include <iostream>
#include <fstream>
#include <string.h>

#define BUFSIZE 256
using namespace std;


double random(double lb, double ub) {
    std::uniform_real_distribution<double> unif(lb, ub);//Type of random number distribution
    std::random_device rand_dev;          // Use random_device to get a random seed.
    std::mt19937 rand_engine(rand_dev()); // mt19937 is a good pseudo-random number generator.
    double x = unif(rand_engine);
    return x;
}

//soft max of abs of vector x, avoiding overflow
int softmaxabs(integer *n, dd_real *mu, dd_real *x, dd_real *f, dd_real *g){
    dd_real z[1024],expz[1024];
    dd_real zmax,mu1;
    dd_real s=dd_real("0.0");
    integer n1,n2;
    
    n1=*n;
    mu1=*mu;
    zmax = dd_real("0.0");
    //computing the zmax
    for (int j = 0; j < n1; ++j) {
        z[j] = x[j]/mu1;
        z[n1+j] = -z[j];
        zmax = std::max(zmax, fabs(z[j]));
    }
    n2=2*n1;
    //shifting z by zmax
    for (int j = 0; j < n2; ++j) {
        z[j] = z[j] - zmax;
        expz[j]  = exp(z[j]);
        s+=expz[j];
    }
    //NOTE: shifting z by zmax when mu ~< 10-5  makes all expz=0  except maximum, which becomes 1.
    //same goes for g
    *f = mu1*zmax+mu1*log(s);
    //std::cout << "f = "<< *f << std::endl;
    for (int j = 0; j < n1; ++j) {
        g[j]= (expz[j] - expz[n1+j])/s;
	    //std::cout << "g[" << j << "] = "  << g[j] << std::endl;
    }
    return 0;
}

int yurileshoushessmoothed(integer *n, dd_real *mu, dd_real *x, dd_real *f, dd_real *g){
    integer n1 = *n;
    dd_real y[1024], gy[1024];
    
    y[0] = x[0];
    for (int j = 0; j < n1 -1; ++j) {
        y[j+1] = x[j+1]-2*x[j];
    }
    softmaxabs(n, mu, y, f, gy);
    g[0] = gy[0];
    for (int j = 0; j < n1 -1; ++j) {
        g[j] -=2*gy[j+1];
        g[j+1] = gy[j+1];
    }
    return 0;
}

void writeoutput(dd_real *fLB, integer *iterLB, dd_real *muvals, int muvalsize){

    ofstream myfile1 ("/Users/fatemehasl/Dropbox/matlabCodes/LBFGS/yuri_leshouches/fLBDD.txt");
    ofstream myfile2 ("/Users/fatemehasl/Dropbox/matlabCodes/LBFGS/yuri_leshouches/itLBDD.txt");
    if (myfile1.is_open() && myfile2.is_open()){
        for(int i = 0; i<muvalsize; i++){
            //cout << "writing value " << value << "to fLBDD file after double cast" << endl;
            myfile1 << std::fixed << std::setprecision(32) << fLB[i] << endl;
            myfile2 << std::fixed << std::setprecision(32) << iterLB[i] << endl;
        }
        myfile1.close();
        myfile2.close();
    }


//    if (myfile2.is_open()){
//        for(int i = 0; i<muvalsize; i++){
//            myfile2 << std::fixed << std::setprecision(32) << iterLB[i] << endl;
//        }
//        myfile2.close();
//    }


}

int lbfgsbcreatemat(double *fLB, integer *iterLB, double *muvals, int muvalsize){


    MATFile *pmat;
    const char *file = "/Users/fatemehasl/Dropbox/matlabCodes/LBFGS/yuri_leshouches/mattestDD.mat";
    printf("Creating file %s...\n\n", file);
    pmat = matOpen(file, "w");
    if (pmat == NULL) {
        printf("Error creating file %s\n", file);
        printf("(Do you have write permission in this directory?)\n");
        return(EXIT_FAILURE);
    }

    mxArray *mval = mxCreateDoubleMatrix(muvalsize, 1, mxREAL);
    double *pmval = mxGetPr(mval);

    mxArray *fLBC = mxCreateDoubleMatrix(muvalsize, 1, mxREAL);
    double *pfLBC = mxGetPr(fLBC);

    mxArray *itLBC = mxCreateDoubleMatrix(muvalsize, 1, mxREAL);
    double *pitLBC = mxGetPr(itLBC);

    int status;

    if (mval == NULL) {
        printf("%s : Out of memory on line %d\n", __FILE__, __LINE__);
        printf("Unable to create mxArray.\n");
        return(EXIT_FAILURE);
    }
    //memcpy((void *)(mxGetPr(pa2)), (void *)data, sizeof(data));
    for(int i = 0; i<muvalsize; i++){
        pmval[i] = muvals[i];
        pfLBC[i] = fLB[i];
        pitLBC[i] = iterLB[i];
    }

    status = matPutVariableAsGlobal(pmat, "muvalDD", mval);
    if (status != 0) {
        printf("Error using matPutVariableAsGlobal\n");
        return(EXIT_FAILURE);
    }
    status = matPutVariableAsGlobal(pmat, "fLBCDD", fLBC);
    if (status != 0) {
        printf("Error using matPutVariableAsGlobal\n");
        return(EXIT_FAILURE);
    }

    status = matPutVariableAsGlobal(pmat, "itLBCDD", itLBC);
    if (status != 0) {
        printf("Error using matPutVariableAsGlobal\n");
        return(EXIT_FAILURE);
    }

    mxDestroyArray(mval);
    mxDestroyArray(fLBC);

    if (matClose(pmat) != 0) {
        printf("Error closing file %s\n",file);
        return(EXIT_FAILURE);
    }
    printf("Done\n");
    return(EXIT_SUCCESS);
}
/*
//Testing softmaxabs gradient to make sure we got it right. Trnaslated from matlab findif by Michael.
//assuming that x is an array.
int findif(integer *n, dd_real *mu, dd_real *x, dd_real *f, dd_real *g){
    
    integer n1 = *n;
    dd_real xpert[1024], fpert , gpert[1024];
    dd_real h ,d[1024], fdif[1024], fderr[1024],f1, gd=0.0;
    int p=16 ,i,j;
    
    softmaxabs(n, mu, x, f, g);
    f1 = *f;
    for (i = 0; i < p; ++i) {
        h =pow((dd_real)10, -(i+1));
        std::cout << "h = "<< h << std::endl;
        gd=0.0;
        
        for (j = 0; j < n1; ++j) {
            d[j] = dd_real(random(0.0,1.0));
            xpert[j] = x[j]+ h*d[j];
            //std::cout << "xpert = "<< xpert[j] << std::endl;
            gd += g[j]*d[j];//directional derivative in d direction
        }
        
        softmaxabs(n, mu, xpert, &fpert, gpert);
        fdif[i] = (fpert-f1)/h;//gradient
        fderr[i] = fabs(fdif[i]-gd);
    }
    findifplot(fderr, p);
    return 0;
}
 */
