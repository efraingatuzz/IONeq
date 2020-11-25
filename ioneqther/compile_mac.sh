#!/bin/bash

# For Mac OSX
rm *.dylib
sed -i.ori "s,local_dir = '.*',local_dir = \'`pwd`\'," ioneqther.f90
echo "initpackage ioneqther lmodel_ioneqther.dat "`pwd`"\nquit\ny" | xspec

rm *~ *.o
rm *FunctionMap.* lpack_* 
rm -f *.mod Makefile
