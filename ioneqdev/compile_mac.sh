#!/bin/bash

# For Mac OSX
rm *.dylib
sed -i.ori "s,local_dir = '.*',local_dir = \'`pwd`\'," ioneq.f90
echo "initpackage ioneq lmodel_ioneq.dat "`pwd`"\nquit\ny" | xspec

rm *~ *.o
rm *FunctionMap.* lpack_* 
rm -f *.mod Makefile
