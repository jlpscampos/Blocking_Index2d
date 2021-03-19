# Two Dimensional Atmospheric Blocking Index #
 Program with the intent to detect **Atmospheric Blocking structures**

This implementation is a 2 dimensional addaptation of the algorithm presented by [Barnes et al., (2012)](https://link.springer.com/article/10.1007/s00382-011-1243-6), to detect atmospheric blockings structures in two atmospheric fields, **Potential Temperature at Dynamical Tropopause**, [Pelly & Hoskins, (2012)](https://journals.ametsoc.org/view/journals/atsc/60/5/1520-0469_2003_060_0743_anpob_2.0.co_2.xml), [Berrisford et al., (2007)](https://journals.ametsoc.org/view/journals/atsc/64/8/jas3984.1.xml) and **Geopotential Height in the Middle Troposphere**, [Tibaldi & Molteni, (1990)](https://onlinelibrary.wiley.com/doi/abs/10.1034/j.1600-0870.1990.t01-2-00003.x).

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

## Description ##

In middle and high latitudes, the upper level winds generally flows from west to east, forming a jet stream due to strengthening of meridional temperature gradient [Holton et al. 2016](https://aapt.scitation.org/doi/pdf/10.1119/1.1987371?casa_token=_TKypuiKE3YAAAAA%3AvNTMbLFXQqVZCLEAu6vWbBb_pO-iFynQe8m4a8d3XXPqTYjvjjD2L9CNsHZCfDP4j7nBRqC8XP0u5g&). When this flow cross some topographical barrier (like the Rokies or the Himalayas) or are influenced by stationary Rossby Waves, it gains potential vorticiy and begins to meander, generating circulation anomalies in the forms of high pressure poleward and lower pressure equatorward, making the flux between the high and low pressures be opposite to the jet stream. Thus, blocking the jet stream, if these anomalies are persisntent by more than 4 or 5 days, it is called **atmospheric blocking** ([Tibaldi & Molteni, (1990)](https://onlinelibrary.wiley.com/doi/abs/10.1034/j.1600-0870.1990.t01-2-00003.x)).

Atmospheric Blocking events can be seen as Rossby Wave breaking events ([Pelly & Hoskins, 2012](https://journals.ametsoc.org/view/journals/atsc/60/5/1520-0469_2003_060_0743_anpob_2.0.co_2.xml)). Strong convection in tropics and in convergence zones can generate vorticity anomalies due to latent heat loss to the atmosphere, exciting Rossby Waves, which propagates in an arch like path from the tropics to extratropics, generating positive and negative vorticity anomalies along its trajectory. High pressure systems acts as Rossby Waves sink ([Takaya & Nakamura, 2000](https://journals.ametsoc.org/view/journals/atsc/58/6/1520-0469_2001_058_0608_afoapi_2.0.co_2.xml?tab_body=fulltext-display)). So, the Rossby Wave generally decays in the form of a high pressure system, where in its northern flank, the opposition to the eastward flux occurs, if this system persist for a period of some days or weeks, we can classify this system as atmospheric blocking as well. (([Berrisford et al.,2007](https://journals.ametsoc.org/view/journals/atsc/64/8/jas3984.1.xml))

In classical works, there are identified two types of atmospheric blockings, the **omega blocking type**, depicted in fig.1a-c and the **dipole blocking type**, depicted in fig.1d-h.

![fig1](https://github.com/jlpscampos/Blocking_Index2d/blob/main/figs/b_all.png)

In summer season over South America, when a blocking event is set over the South Atlantic coast, it is expected heat waves over the continent, drougth conditions northwestward and wetness condition southwestward the blocking high, due to the meridional displacement of the South American Low Level Jet and the impediment of frontal system to propagate equatorward, due to the blocking high and the southward displacement of the Jet Stream.  

Before formulating the algorithm to detect an atmospheric blocking, some quantities must be computed before, such as the meridional gradient of potential temperature or geopotential high, using the potential temperature at the dynamical tropopause (\theta) and geopotential high at 500hPa level (Z) respectivelly, defined bellow:

![eqn1](https://github.com/jlpscampos/Blocking_Index2d/blob/main/figs/eqn1.png)

Equations 1 and 2 represent the quantities to be computed with the geopotential height and potential temperature respectivelly. To a flow be considered a blocked flow, the following criteria must be verified, we call this of **Instantaneous Blocked Latitude** (IBL).

![eqn1](https://github.com/jlpscampos/Blocking_Index2d/blob/main/figs/eqn2.png)

![Blocking](https://github.com/jlpscampos/Blocking_Index2d/blob/main/figs/blocking_19830126-19830204_full.png)
