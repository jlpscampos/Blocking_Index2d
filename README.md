# Two Dimensional Blocking Index #
 Program with the intent to detect **Atmospheric Blocking structures**

This implementation is an addaptation of the algorith presented in [Barnes et al., (2012)](https://link.springer.com/article/10.1007/s00382-011-1243-6), to detect atmospheric blockings structure in two fields, **Potential Temperature at Dynamical Tropopause**, [Pelly & Hoskins, (2012)](https://journals.ametsoc.org/view/journals/atsc/60/5/1520-0469_2003_060_0743_anpob_2.0.co_2.xml), [Berrisford et al., (2007)](https://journals.ametsoc.org/view/journals/atsc/64/8/jas3984.1.xml) and **Geopotential Height in the Middle Troposphere**, [Tibaldi & Molteni, (1990)](https://onlinelibrary.wiley.com/doi/abs/10.1034/j.1600-0870.1990.t01-2-00003.x).

## How to run ##

* (1) Before runnning you need to set up with input/output directories and begin and end dates in the file **config.fortran.txt** in folder *pwd/config/*  
* (2) The fortran95 version of this implementation is intended to read monthly netcdf data issued 6 by 6 hours. For other configurations, you must to addapt the main program.
* (3) You need to install the netcdf libraries in your UNIX system before running. 
* (4) The compile options for fortran95 implementation is in file *main_f95.sh*
>> ```bash
>> #!/bin/bash
>>FC="gfortran"
>>fflags="-fdefault-real-8 -fconvert=big-endian -frecord-marker=4 -w  -O3"
>>LIBNETCDF="-I/usr/include/ -L/usr/lib -lnetcdf -lnetcdff"
>>pwdd=`pwd`
>>$FC $fflags -c src_fortran/Blocked_Flow_Index_new_v2.f95 $LIBNETCDF
>>$FC $fflags -c main.f95 $LIBNETCDF
>>$FC $fflags main.f95 src_fortran/Blocked_Flow_Index_new_v2.f95 -o a.out $LIBNETCDF -I/$pwdd/src_fortran/
>>./a.out
>>```
