#!/bin/sh

# Documentation taken from:
# https://omniweb.gsfc.nasa.gov/html/HROdocum.html
# res=min and spacecraft=omnbi_min were taken from looking at the html
#   (i.e., documentation didn't cover this!)
#
# These are the variables needed.  They are off by 1 from
# the above the documentation:
# vars=14   -> BX, nT (GSE, GSM)
# vars=17   -> BY, nT (GSM)
# vars=18   -> BZ, nT (GSM)
# vars=22   -> Vx Velocity,km/s
# vars=23   -> Vy Velocity, km/s
# vars=24   -> Vz Velocity, km/s
# vars=25   -> Proton Density, n/cc
# vars=26   -> Temperature, K
# vars=37   -> AE-index, nT
# vars=38   -> AL-index, nT
# vars=39   -> AU-index, nT
#


curl -d  "activity=retrieve&res=min&spacecraft=omni_min&start_date=20110620&end_date=20110623&vars=14&vars=17&vars=18&vars=22&vars=23&vars=24&vars=25&vars=26&vars=37&vars=38&vars=39&back=" https://omniweb.gsfc.nasa.gov/cgi/nx1.cgi > omni_20110622.txt
#curl -d  "activity=retrieve&res=min&spacecraft=omni_min&start_date=20050101&end_date=20050105&vars=14&vars=17&vars=18&vars=22&vars=23&vars=24&vars=25&vars=26&vars=37&vars=38&vars=39&back=" https://omniweb.gsfc.nasa.gov/cgi/nx1.cgi > test_curl_min.txt
