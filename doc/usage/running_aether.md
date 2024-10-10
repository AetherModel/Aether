# Running Aether

This document assumes you have already downloaded and built the Aether model. If
not, you should return to [one](../../README.md) of
[these](../installation/installation.md) pages before continuing.

## The first run

Once you have compiled you can run Aether. To remember which runs you're doing,
we recommend creating a specific directory for each run.  This can be done by
copying the default run directory to have a special name, like this:

```bash
cd ..
cp -R share/run ./run.first_run
```

This creates the directory where you will do your run.  In that directory is a
link to the aether executable and an input file called aether.json.

You can then run the executable from the directory you created.

```bash
cd run.first_run
./aether
```

You should see something like:

```bash
run.first_run% ./aether
> Need to NOT adjust F10.7, but that isn't included yet!!!
> Writing file : 3DALL_20110320_000000
> Writing file : 3DBFI_20110320_000000
> Wall Time : 4s (left : 1111h); Current Time : 2011 3 20 0 0 0 0 
> Wall Time : 4s (left : 23m); Current Time : 2011 3 20 0 0 10 0 
> Wall Time : 5s (left : 14m); Current Time : 2011 3 20 0 0 20 0 
> Wall Time : 5s (left : 9m); Current Time : 2011 3 20 0 0 30 0 
> Wall Time : 5s (left : 7m); Current Time : 2011 3 20 0 0 40 0 
> Wall Time : 6s (left : 7m); Current Time : 2011 3 20 0 0 50 0 
> Wall Time : 6s (left : 5m); Current Time : 2011 3 20 0 1 0 0 
> Wall Time : 6s (left : 5m); Current Time : 2011 3 20 0 1 10 0 
> Wall Time : 7s (left : 5m); Current Time : 2011 3 20 0 1 20 0 
> Wall Time : 7s (left : 4m); Current Time : 2011 3 20 0 1 30 0 
etc.
```

The successful end of this will show a timing summary, similar to:

```bash
>main>advance>Chemistry::calc_chemistry
    nTimes called : 720
    timing_total (s) : 1.52599
>main>advance>Ions::calc_ion_temperature
    nTimes called : 720
    timing_total (s) : 19.0279
>main>advance>Neutrals::exchange
    nTimes called : 720
    timing_total (s) : 2.94398
```

## Output Files

Aether outputs to a subdirectory called UA/output. At this time, all processors
within the Aether run output their own files and they can be combined using the
post processor.

### Blocks

Aether is a block-based code, so the domain is split up into different blocks
that hold a given volume (for spherical grids, this is some range of latitude,
longitude and altitude; for other grids it could differ). If you asked for (for
example) 4 processors, then the Earth would be broken up into 4 different blocks
and each output would consist of 4 different blocks.

### Ensembles

Aether can run the same simulation multiple times simultaneously with different
drivers and internal parameters.  The collection of these simulations is called
an ensemble, with each individual simulation called a member.  All of the
ensemble members output to the UA/output directory.  See the ensembles.md
document to learn more about ensembles.

### Post processing

Because there can be multiple blocks and multiple ensemble members, there is a
post processor that can pull all of the files together. By default, the post
processor combines all of the block files for one ensemble member for one time
into one file.  This means that if there are 25 hourly outputs requested, there
will be 25 total files for each ensemble member.

If there is more than one simulation, then the post processor will create a mean
file for each time.  If there are more than two ensemble members, a standard
deviation file is also created.  The mean and standard deviations are calculated
across ensemble members.

To run the post processor, do the following:

```bash
cd UA/output
../../../srcPython/post_process.py
```

By default, this will also leave the raw files and will produce plots just to
make sure that everything is working ok.  If you don't need the raw files (there
is not a good reason to keep these, unless you are debugging the post
processor), you can run the post-processor with the "-rm" argument.  If you
don't want any plots, you can run with the "-alt=-1" argument.  An example of
this is (this is how it is often run):

```bash
../../../srcPython/post_process.py -rm -alt=-1
```

If you are going to be using Aether a lot, you may want to copy the
post-processor into your bin directory.

## Input Files

Aether reads in a bunch of files, most of which are specified by the settings
file.  In order to minimize the number of places where new settings need to be
specified, Aether uses a defaults file that sets the generic defaults of the
model.  This file is in UA/inputs/defaults.json.  

## defaults.json file

This is a json file that sets all of the defaults within Aether.  This file
should never be modified!

### For Developers

Within Aether, the inputs.cpp file has a large handful of of get_ routines to
get the values of the settings that the user has set.

## aether.json file

The file aether.json is read in AFTER the defaults file and these settings
overwrite the defaults.  So, if you want to modify the default settings, simply
copy a setting out of defaults.json and paste it into aether.json.  Then, you
can modify the setting in aether.json.

Because the files are json files, you don't actually have to set all of the
subsettings within a specific setting.  For example, within the defaults.json
file, the EUV setting is:

```bash
    "Euv" : {
        "doUse" : true,
        "Model" : "euvac",
        "File" : "UA/inputs/euv.csv",
        "IncludePhotoElectrons" : true,
        "HeatingEfficiency" : 0.05,
        "dt" : 60.0},
```

If you simply want to turn off the EUV, you can insert this into the aether.json
file:

```bash
    "Euv" : {
  "doUse" : false},
```

or, if you wanted to use the 59 wavelength bins instead of the 37 EUVAC bins,
you could do:

```bash
    "Euv" : {
  "File" : "UA/inputs/euv_59.csv"},
```

## planet.in file

Within the input file, the planet file can be set:

```bash
    "Planet" : {
        "name" : "earth",
        "file": "UA/inputs/earth.in"},
```

This file defines the neutrals and ions to be modeled. For example:

```bash
#NEUTRALS
name, mass, vibration, thermal_cond, thermal_exp, advect, BC
O, 16, 5, 5.6e-4, 0.69, 1, 5.0e17
N, 14, 5, 5.6e-4, 0.69, 1, 2.5e11
O2, 32, 7, 3.6e-4, 0.69, 1, 4.0e18
N2, 28, 7, 3.6e-4, 0.69, 1, 1.7e19
NO, 30, 7, 3.6e-4, 0.69, 1, 5.4e14
He, 4, 5, 5.6e-4, 0.69, 1, 1.5e14
N_2D, 14, 5, 5.6e-4, 0.69, 0, 2.5e9
N_2P, 14, 5, 5.6e-4, 0.69, 0, 2.5e7
H, 1, 5, 5.6e-4, 0.69, 0, 3.0e13
O_1D, 16, 5, 5.6e-4, 0.69, 0, 1.0e10
CO2, 46, 7, 3.6e-4, 0.69, 0, 4.5e15
```

This is a comma seperated list that includes:

- name - this is the common name of the species
- mass - the atomic mass of the species
- vibration - the degrees of freedom for the species
- thermal_cond - the thermal conductivity coeficient of the species
- thermal_exp - the exponent that is put on the temperature (Lambda = A * T^B)
- advect - whether the species is advected or not
- BC - the lower boundary condition on the species density

```bash
#IONS
name, mass, charge, advect
O+, 16, 1, 1
O2+, 32, 1, 0
N2+, 28, 1, 0
NO+, 30, 1, 0
N+, 14, 1, 0
He+, 4, 1, 1
O+_2D, 16, 1, 0
O+_2P, 16, 1, 0
```

This is a comma separated list that includes:

- name - this is the common name of the species
- mass - the atomic mass of the species
- charge - the charge of the species
- advect - whether the species is advected or not

```bash
#TEMPERATURE
alt, temp
100.0, 200.0
200.0, 800.0
300.0, 800.0
```

Within the input file (i.e., aether.json), the following can be set:

```bash
    "InitialConditions" : {
        "type" : "Planet"},    
    
    "BoundaryConditions" : {
        "type" : "Planet"},
```

If these are set to "Planet", then the temperature profile is set as the initial
condition.  If the boundary condition is set to "Planet", then the lowest
altitude is used as a boundary condition on the temperature.

## orbits.csv file

This file contains all of the orbital, mass, rotation, and magnetic field
characteristics of the planets.  All of the planets within the solar system are
included.  Other planets (like exo-planets or artificial planets) can be
included by adding lines.  

## chemistry file

The chemistry file defines all of the chemical reactions within the system.

For example:

```bash
He+ + e- -> He, with  R = 4.8e-18 * (250 / Te) ^ 0.7
```

With an exothermal energy of xxx. Also, the uncertainty of the reaction can be
set too (set to 10%, or 0.1).

```bash
R11,He+,e-,,,He,,,4.80E-18,(250/Te)^0.7,,1,0,,0.1,,,,gitm,250,Te,0.7,,,,1,
```
