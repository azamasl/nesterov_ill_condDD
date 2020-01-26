# nesterov_ill_condDD
In this project we investigate the failure of L-BFGS on the smoothed version of the following very ill-conditioned nonsmooth function:

f(x) = max{|x_1|, |x_{i+1} - x_{i}|, i=1,...,n-1}

by using double double precision.  The original L-BFGS-B-C implemented by Stephen Becker and is  available from 
 
 https://github.com/stephenbeckr/L-BFGS-B-C
 
 and the Quad Double computation package is available from 
 
 https://www.davidhbailey.com/dhbsoftware/

