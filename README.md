# floods

This repository contains data and R code for robust local likelihood estimation for nonstationary flood frequency analysis. The data and code presented here are used in the 2024 paper by John Grego and Philip Yates, "Robust Local Likelihood Estimation for Non-stationary Flood Frequency Analysis," appearing in the _Journal of Agricultural, Biological, and Environmental Statistics_.

To perform the robust local likelihood estimation:

  - For the Congaree River: Read the file `peak.txt` into R. This file contains the annual maximum discharge for the Congaree River at Columbia, SC (USGS Gage 02169500) for the flood years 1892 to 2022. The data can also be found on the [USGS Water website](https://nwis.waterdata.usgs.gov/sc/nwis/peak?site_no=02169500&agency_cd=USGS&format=html).
  - For the Illinois River: Read the file `illinois.csv` into R. This file contains the annual maximum discharge for the Illinois River at Marseilles, IL (USGS Gage 05543500) for the flood years 1892 to 2022. The data can also be found on the [USGS Water Website](https://nwis.waterdata.usgs.gov/il/nwis/peak?site_no=05543500&agency_cd=USGS&format=html).
  - For the Winooski River: Read the file `winooski.csv` into R. This file contains the annual maximum discharge for the Winooski River at Montpelier, VT (USGS Gage 04286000) for the flood years 1912 to 2023. The data can also be found on the [USGS Water Website](https://nwis.waterdata.usgs.gov/vt/nwis/peak?site_no=04286000&agency_cd=USGS&format=html).
  - The program `RLLgev.R` produces two components:
    - `parm.est` which for each flood year produces parameter estimates for the robust local likelihood model as well as estimates for the 99<sup>th</sup> percentile, mean, and median
    - `wts` which contains the robust weights for each observation used in a flood year's robust local likelihood model; each **column** represents the robust weights associated with a given flood year while each **row** represents the weights used in a flood year's smoothing window.
  - The program `RLLqsd.R` produces the standard error estimate for the 99<sup>th</sup> percentile for a specific flood year. This program uses the program `RLL.sandwich.mat.R`.

 See the program `steps for analysis.R` on the main page to see how to perform the robust local likelihood estimation.

 The program `l.cv.R` can be used to find the cross-validated log-likelihood for a given bandwidth, h.

 Results are found in the CSV file `congareeRLL.csv`, `illinoisRLL.csv`, and `winooskiRLL.csv`, respectively.
    
   
      
