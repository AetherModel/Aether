# FISM

FISM can be used as a EUV model.

Needs the FISM files automatically made by `srcPython/fism.py`.

FISM contains the binned flux already, however Aether must still be provided with an
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
