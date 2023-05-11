
The "iVerbose" command can be used under "Debug" in aether.json to set
the overall verbose level in the code. Optionally, users can also
assign specific verbose levels to be used for certain functions by
specifying the function names and the corresponding verbose levels
using the "iFunctionVerbose" command as follows:

    "Debug" : {
    "iVerbose" : 0,
    "iFunctionVerbose" : {
        "func1" : 1,
        "func2" : 2},
    "dt" : 10.0
    }

When a sub-function is entered, the verbose level stays the same if
the verbose level for that sub-function is not specified. When that
sub-function exits, the verbose level will accordingly be set back to
that of the function it is returning to.