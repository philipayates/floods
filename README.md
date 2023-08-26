# floods

This repository contains data and R code for robust local likelihood estimation for nonstationary flood frequency analysis. The data and code presented here are used in the (submitted 2023) paper by John Grego and Philip Yates.

To perform the robust local likelihood estimation:

  - Read the file `peak.txt` into R. This file contains the annual maximum discharge for the Congaree River at Columbia, SC (USGS Gage 02169500) for the flood years 1892 to 2022. The data can also be found on the [USGS Water website](https://nwis.waterdata.usgs.gov/sc/nwis/peak?site_no=02169500&agency_cd=USGS&format=html).
  - The program `RLLgev.R` produces two components:
    - `parm.est` which for each flood year produces parameter estimates for the robust local likelihood model as well as estimates for the 99<sup>th</sup> percentile, mean, and median
    - `wts` which contains the robust weights for each observation used in a flood year's robust local likelihood model; each **column** represents a flood year
  - The program `RLLqsd.R` produces the standard error estimate for the 99<sup>th</sup> percentile for a specific flood year. This program uses the program `RLL.sandwich.mat.R`.

 See the program `steps for analysis.R` on the main page to see how to perform the robust local likelihood estimation.
    
   
      
