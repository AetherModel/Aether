# FISM-2

FISM-2 can be used as a EUV model.

Needs the FISM-2 files automatically made by `srcPython/fism.py`.

FISM-2 contains the binned flux already, however Aether must still be provided with an
EUV csv file containing the cross-sections.

The different fism models available in fism.py contain different numbers of bins, so
the euv file must be different.

|  Model  |  number of bins  | euv file |
|  :---  |  :-------------: | -------: |
| HFG    | 23  | euv_solomon.csv |
| Solomon | 23  | euv_solomon.csv |
| NEUVAC  | 37/59   | euv.csv / euv_59.csv |
| EUVAC  |  37  | euv.csv |

The input format, when using fism data:

    "Euv" : {
        "doUse" : true,
        "Model" : "fism",
        "File" : "UA/inputs/euv_59.csv",
        "fismFile": "fism2_file_59.txt",
        "IncludePhotoElectrons" : true,
        "HeatingEfficiency" : 0.05,
        "dt" : 60.0
    },

To generate FISM-2 irradiances between two dates, one may call fism.py like so:

srcPython/fism.py 20110319 20110321 -b neuvac

Note that the optional argument '-b' defaults to the binning scheme of the NEUVAC model,
which employs 59 bins. To use the 37 bins of the EUVAC model, the argument should be
'euvac', while to use the 23 bins used by the HFG model, the argument should be 'solomon'.

fism.py should always be run before Aether is run using the FISM-2 model. Even though the
entire FISM-2 irradiances are stored in the repository, Aether will need a separate .csv
file containing the temporal subset of FISM-2 irradiances in the desired binning scheme.
