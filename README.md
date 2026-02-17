# Rfast

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/Rfast)](https://cran.r-project.org/package=Rfast) [![metacran downloads](https://cranlogs.r-pkg.org/badges/grand-total/Rfast)](https://cran.r-project.org/package=Rfast) [![metacran downloads](https://cranlogs.r-pkg.org/badges/Rfast)](https://cran.r-project.org/package=Rfast) [![CRAN_latest_release_date](https://www.r-pkg.org/badges/last-release/Rfast)](https://cran.r-project.org/package=Rfast) [![Stack Overflow](https://img.shields.io/badge/stackoverflow-Rfast-blue.svg)](https://stackoverflow.com/search?q=Rfast)


## 1. About
A collection of fast (utility) functions for data analysis. Column- and row- wise means, medians, variances, minimums, maximums, many t, F and G-square tests, many regressions (normal, logistic, Poisson), are some of the many fast functions. 

## 2. Note
* The vast majority of the functions accept matrices only, not data.frames.
* Do not have matrices or vectors with have missing data (i.e NAs). We do no check about them and C++ internally transforms them into zeros (0), so you may get wrong results.
* In general, make sure you give the correct input, in order to get the correct output. We do no checks and this is one of the many reasons we are fast.

## 3. References
*  Tsagris M., Papadakis M. (2018). Taking R to its limits: 70+ tips. PeerJ Preprints 6:e26605v1 <\doi:10.7287/peerj.preprints.26605v1>.
*  Tsagris M. and Papadakis M. (2018). Forward regression in R: from the extreme slow to the extreme fast. Journal of Data Science, 16(4): 771–780. <\doi:10.6339/JDS.201810_16(4).00006>.

## Donate

[![paypal](https://www.paypalobjects.com/en_US/i/btn/btn_donateCC_LG.gif)](https://www.paypal.com/paypalme/rfastofficial)
