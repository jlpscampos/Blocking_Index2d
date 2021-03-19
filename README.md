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

## About the Algorithm and Atmospheric Blocking ##

In summer season over South America, when a blocking event is set over the South Atlantic coast, it is expected a poleward displacement of the climatological humidity corridor, the convergence zone, leading to droughts over the most populated South America megalopolis, the São Paulo - Rio de Janeiro region  and wetter weather over Buenos Aires/Mar del La Plata region. Beyond the impacts in water supply and agriculture, hydroeletric and renewable energy markets are strong impacted as well, leading business to bankrupt if the blocking is not forecasted. So, the automatic identification can help in decision making by energy traders and civil defense professionals. This repository aims to provide a tool capable of identify atmospheric blockings.

In middle and high latitudes, the upper level winds generally flows from west to east, forming a jet stream due to strengthening of meridional temperature gradient [Holton et al. 2016](https://aapt.scitation.org/doi/pdf/10.1119/1.1987371?casa_token=_TKypuiKE3YAAAAA%3AvNTMbLFXQqVZCLEAu6vWbBb_pO-iFynQe8m4a8d3XXPqTYjvjjD2L9CNsHZCfDP4j7nBRqC8XP0u5g&). When this flow cross some topographical barrier (like the Rokies or the Himalayas) or are influenced by stationary Rossby Waves, it gains potential vorticiy and starts to meander, generating circulation patterns in the form of high pressure and lower pressure anomalies. When the configuration of this anomalies leads to an eastward flux opposing to the jet stream and last more than 3 or 5 days, we call these events as an **atmospheric blocking** event ([Tibaldi & Molteni, 1990](https://onlinelibrary.wiley.com/doi/abs/10.1034/j.1600-0870.1990.t01-2-00003.x)).

Atmospheric Blocking events can be seen as Rossby Wave breaking events ([Pelly & Hoskins, 2012](https://journals.ametsoc.org/view/journals/atsc/60/5/1520-0469_2003_060_0743_anpob_2.0.co_2.xml)) as well. Strong convection in tropics and in convergence zones (both leading upper level divergence) can generate vorticity anomalies due to latent heat loss to the atmosphere, exciting Rossby Waves, which propagates in an arch like path from tropics to extratropics, generating positive and negative vorticity anomalies along its trajectory. High pressure systems acts as Rossby Waves sink ([Takaya & Nakamura, 2000](https://journals.ametsoc.org/view/journals/atsc/58/6/1520-0469_2001_058_0608_afoapi_2.0.co_2.xml?tab_body=fulltext-display)), where in each high pressure formed along the Rossby Wave path, the energy is vanished, up to it decays in the form of a dipole type blocking, where a high pressure is formed poleward of a lower pressure system.

In classical works, there are identified two types of atmospheric blockings, the **omega blocking type**, depicted in fig.1a-c and the **dipole blocking type**, depicted in fig.1d-h. So, the purpose of this repository is to provide a set of routines capable to identify these kinds of structures.

![fig1](https://github.com/jlpscampos/Blocking_Index2d/blob/main/figs/b_all.png)

Before formulating the algorithm to detect an atmospheric blocking, some quantities must be computed before, such as the meridional gradient of potential temperature or geopotential high, using the potential temperature at the dynamical tropopause (\theta) and geopotential high at 500hPa level (Z) respectivelly, defined bellow:

![eqn1](https://github.com/jlpscampos/Blocking_Index2d/blob/main/figs/eqn1.png)

Equations 1 and 2 represent the quantities to be computed with the geopotential height and potential temperature respectivelly. To a flow be considered a blocked flow, the following criteria must be verified, we call this of **Instantaneous Blocked Latitude** (IBL).

![eqn1](https://github.com/jlpscampos/Blocking_Index2d/blob/main/figs/eqn2.png)

![Blocking](https://github.com/jlpscampos/Blocking_Index2d/blob/main/figs/blocking_19830126-19830204_full.png)
