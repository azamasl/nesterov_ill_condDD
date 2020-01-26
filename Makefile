#  Created by Azam Asl on 12/9/16.
#  Copyright © 2016 Azam	 Asl. All rights reserved.
CC = gcc

CFLAGS = -g -Wall -Wno-uninitialized
LDFLAGS = -lqd -leng -lmx -lm -lmat -lut -lstdc++

# matlab library:
LIB = -L/Applications/MATLAB_R2015a.app/bin/maci64 -L/usr/local/lib

# matlab and Double-Double includes:
INCLUDES = -I/Applications/MATLAB_R2015a.app/extern/include/ -I/Users/fatemehasl/projects/BFGS/qd/include

LBFGSB  = lbfgsb.cpp linesearch.cpp subalgorithms.cpp print.cpp linpack.cpp miniCBLAS.cpp timer.cpp yurileshouchessm.cpp
#TESTFUNCS = yurileshouchessm.cpp
TESTFUNCS = T10_dd.cpp
MAIN =test2.cpp

SRC = $(LBFGSB) $(TESTFUNCS)

default: all exetest
all :  test

test :
	$(CC) $(CFLAGS) $(MAIN) $(SRC) $(INCLUDES) $(LIB) $(LDFLAGS) -o x

exetest : x
	./x
