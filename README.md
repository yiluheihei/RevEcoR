# Reverse Ecology analysis in  R

This package implementation the applications of reverse ecology. 
Reverse ecology refers to the use of genomics to study ecology with no a priori
assumptions about the organism(s) under consideration, linking the organism and
their environment and the interaction among species.

**Prerequisites**
Until **ReverseEcologyR** is ready for CRAN or Bioconductor, you can install it directly from GitHub using devtools.

However, `install_github` can only install the dependencies that have been released on CRAN. Three dependencencies for 
`ReveseEcologyR were released on **Bioconductor**£º`Biobase`, `KEGGREST` and `mmnet`, and vignette of `ReverseEcologyR` 
was built with `knitr`. Users should install these four packages manually first. 

```r
source("http://bioconductor.org/biocLite.R")
biocLite(c("mmnet","KEGGREST","Biobase")) 
if(!require(knitr))
 install.packages("knitr") 
```
Then install `ReverseEcologyR` using devtools.
```{r}
if (!require(devtools)£©
  install.packages("devtools")
devtools::install_github("yiluheihei/ReverseEcologyR")
```

You'll also need to make sure your machine is able to build packages from source. See [Package Development Prerequisites]
(https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites) for the tools needed for your operating system.

After installation, you can load **ReverseEcologyR** into current workspace by typing or pasting the following codes:

 ```R
library("ReverseEcologyR")
 ```
## Contributing

For very simple changes such as fixing typos, you can just edit the file by clicking the button `Edit`. 
For more complicated changes, you will have to manually create a pull request after forking this repository.
 
##License

`ReverseEcologyR` is a free and open source software, licensed under GPL.
