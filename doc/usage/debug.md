# Debug

The Debug command in the input file sets the type and amount of
information that is fed to the user. An example of the Debug command
is:

```json
    "Debug" : {
        "iVerbose" : 0,
        "dt" : 60.0,
        "TimingPercent" : 1.0,
        "iTimingDepth" : 5,
        "iProc" : 0},
```

The "iVerbose" command can be used under "Debug" in aether.json to set
the overall verbose level in the code.  This sets the amount of
debugging information that is output.  In addition, many functions
report when they are entered (and existed).  These functions
automatically report when the depth of the function is less than
iVerbose.

Because Aether is a parallel code, users typically don't want to see
output from all processors, so the user can specify which processor
outputs using "iProc".

"dt" specifies how often to report timing information and is in
seconds in simulation time (not walltime).

The "Timing" variables specify how much information is passed to the
user for code profiling at the end of the simulation.  These variables
limit the amount of information.  For example, "iTimingDepth" limits
the reporting of timing for functions by depth of call.
"TimingPercent" specifies to report the timing only if the time that
the function took was larger than the TimingPercent of the total
run-time.

Users can assign specific verbose levels to be used for
certain functions by specifying the function names and the
corresponding verbose levels using the "iFunctionVerbose" command as
follows:

```json
    "Debug" : {
    "iVerbose" : 0,
    "iFunctionVerbose" : {
        "func1" : 1,
        "func2" : 2},
    "dt" : 10.0
    }
```

When a sub-function is entered, the verbose level stays the same if
the verbose level for that sub-function is not specified. When that
sub-function exits, the verbose level will accordingly be set back to
that of the function it is returning to.
