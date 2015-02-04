# Reverse Ecology analysis in  R

This package implementation the applications of reverse ecology. 
Reverse ecology refers to the use of genomics to study ecology with no a priori
assumptions about the organism(s) under consideration, linking the organism and
their environment and the interaction among species.

**Prerequisites**

Until **RevEcoR** is ready for CRAN or Bioconductor, you can install it directly from GitHub using devtools.

However, `install_github` can only install the dependencies that have been released on CRAN. Three dependencencies for `RevEcoR` were released on **Bioconductor**, including 
`Biobase`, `KEGGREST` and `mmnet`. Moreover, the vignette of `RevEcoR` 
was built using package `knitr`. Users should install these four packages manually first. 

```r
source("http://bioconductor.org/biocLite.R")
biocLite(c("mmnet","KEGGREST","Biobase")) 
if(!require(knitr))
 install.packages("knitr") 
```
Then install `RevEcoR` using devtools.

```r
if (!require(devtools)
  install.packages("devtools")
devtools::install_github("yiluheihei/RevEcoR")
```

You'll also need to make sure your machine is able to build packages from source. Rtools

After installation, you can load **RevEcoR** into current workspace by typing or pasting the following codes:

 ```R
library("RevEcoR")
 ```
## Contributing

For very simple changes such as fixing typos, you can just edit the file by clicking the button `Edit`. 
For more complicated changes, you will have to manually create a pull request after forking this repository.
 
##License

`RevEcoR` is a free and open source software, licensed under GPL.
