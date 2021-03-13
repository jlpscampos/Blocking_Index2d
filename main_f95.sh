#!/bin/bash

FC="gfortran"
fflags="-fdefault-real-8 -fconvert=big-endian -frecord-marker=4 -w  -O3"
LIBNETCDF="-I/usr/include/ -L/usr/lib -lnetcdf -lnetcdff"
pwdd=`pwd`

$FC $fflags -c src_fortran/Blocked_Flow_Index_new_v2.f95 $LIBNETCDF
$FC $fflags -c main.f95 $LIBNETCDF
$FC $fflags main.f95 src_fortran/Blocked_Flow_Index_new_v2.f95 -o a.out $LIBNETCDF -I/$pwdd/src_fortran/
./a.out

# removing some files created in running main.f95.
# If you want to keep them comment or erase the lines bellow
rm a.out
rm *.mod
rm *.o
