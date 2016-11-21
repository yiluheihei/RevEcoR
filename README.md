# Reverse Ecology analysis in  R

[![Travis-CI Build Status](https://travis-ci.org/yiluheihei/RevEcoR.svg?branch=master)](https://travis-ci.org/yiluheihei/RevEcoR)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/RevEcoR)](https://cran.r-project.org/package=RevEcoR)

This package implements the reverse ecology framework. 
Reverse ecology refers to the use of genomics to study ecology with no a priori
assumptions about the organism(s) under consideration, linking the organism and
their environment and the interactions among species.

**Prerequisites**

**RevEcoR** is free available on CRAN. You can install the latest released 
version from CRAN as following:

```{r,eval=FALSE} 
install.packages("RevEcoR") 
```

or the latest development version from github. To install packages from GitHub,
you first need install the **devtools** package on your system with 
`install.packages("devtools")`. Note that devtools sometimes needs some 
extra non-R software on your system -- more specifically, an Rtools download for
Windows or Xcode for OS X. There's more information about devtools
[here](https://github.com/hadley/devtools).
  
```{r,eval=FALSE} 
if (!require(devtools) 
  install.packages("devtools") 
devtools::install_github("yiluheihei/RevEcoR") 
```


After installation, you can load **RevEcoR** into current workspace by typing or pasting the following codes:

 ```R
library("RevEcoR")
 ```
## Contributing

For very simple changes such as fixing typos, you can just edit the file by clicking the button `Edit`. 
For more complicated changes, you will have to manually create a pull request after forking this repository.
 
##License

`RevEcoR` is a free and open source software, licensed under GPL 2.0.
