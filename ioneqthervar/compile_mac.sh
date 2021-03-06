#!/bin/bash

# For Mac OSX
rm *.dylib
sed -i.ori "s,local_dir = '.*',local_dir = \'`pwd`\'," ioneqthervar.f90
echo "initpackage ioneqthervar lmodel_ioneqthervar.dat "`pwd`"\nquit\ny" | xspec

rm *~ *.o
rm *FunctionMap.* lpack_* 
rm -f *.mod Makefile
