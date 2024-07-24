# Ensembles

Aether is capable of running with ensembles, meaning that you can run
multiples of the same simulation in one run. This is enabled by adding
to the aether.json file:

```json
    "Ensembles" : {
        "nMembers" : N},
```

where N is the number of members that you would like to run.  In order
to do this, the number of nodes that are requested have to be N times
the number of nodes that you need for one simulation.

The default is 1, so only one ensemble member.

In order to make the ensemble members unique from each other, you can
perturb the indices and internal parameterizations that Aether is run
with.

To perturb indices in Aether, you can use the "Perturb" command in the
aether.json file. As an example:

```json
"Perturb": {
    "f107" : { "Mean" : 1.0,
               "Std" : 0.10,
               "Add" : false,
               "Constant" : true},
    "f107a" : { "Mean" : 1.0,
                "Std" : 0.10,
                "Add" : false,
                "Constant" : true},
    "imfbz" : { "Mean" : 0.0,
                "Std" : 2.0,
                "Add" : true,
                "Constant" : false}
    }
```

This perturbs the F10.7, F10.7a, and IMF Bz indices.  The sub-parts of
the command include:

The random value (perturbation) is given as the mean + value, where
value is the value that comes out of the random number generator, with
the std fed in as the standard deviation of the random numbers
desired.

Add: The tells Aether whether you want to add perturbations or
multiply perturbations. Adding perturbations (add = true) is good for
indices that have both positive and negative values. Multiplying (add
= false) is good for indices that have to be positive.

Mean: This is the bias that is either added or multiplied.  For an
unbiased distribution, use a mean of 0 for addition and 1.0 for
multiplication.

Std: Standard deviation of the perturbations.  For addition, this
number could, in theory, be whatever the user would like.  For
multiplication, if this number becomes too close to one (1.0), there
could be negative values, which may cause issues (i.e., if the std is
0.33, then there is a 1% chance that the output from the random number
generator will be less than -1.0, which would make the mean + value
less than zero.)

Constant: perturb all of the times by the same value.  If this is used
for indices, then they will be biased by the perturbed value,
resulting in a biased ensemble member. The entire ensemble will be
unbiased if there are enough members and the means are set as
unbiased.  Having non-constant values (Constant = false), adds
(non-biased) noise to the index.

When Aether is run, it will create and perturb all of the ensemble
members, outputting member files with the addition of "_m0NML", where
NML is the member number.  These files can be viewed on their own, or
the post processing code can create mean and std files.

You can also perturb chemical reaction rates, but doing something like:

```json
    "Perturb": {
        "Chemistry" : ["R2", "R10", "R20"] }
```

This will perturb the reaction rates in the chemistry file.  The Rxx
names have to be in the "name" column.  The percentage of uncertainty
is given in the "uncertainty" column and is expressed as a percentage.
This is used as the standard deviation in the perturbation.  In the
language used above, this perturbation is multiplicative with a mean
of 1.0 and a standard deviation set as the uncertainty.
