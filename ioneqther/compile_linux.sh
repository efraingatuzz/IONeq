#!/bin/bash

# For Linux
rm lib*
sed -i "s,local_dir = '.*',local_dir = \'`pwd`\'," ioneqther.f90
echo "initpackage ioneqther lmodel_ioneqther.dat "`pwd`"\nquit\ny" | xspec

rm *~ *.o
rm *FunctionMap.* lpack_* 
rm -f *.mod Makefile
