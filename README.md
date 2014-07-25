# Reverse Ecology analysis in  R

This package implementation the applications of reverse ecology. 
Reverse ecology refers to the use of genomics to study ecology with no a priori
assumptions about the organism(s) under consideration, linking the organism and
their environment and the interaction among species.

Until **ReverseEcologyR** is ready for CRAN or Bioconductor, you can install it directly from GitHub using devtools.

```{r}
install.packages("devtools")
devtools::install_github("yiluheihei/ReverseEcologyR")
```

You'll also need to make sure your machine is able to build packages from source. See [Package Development Prerequisites](https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites) for the tools needed for your operating system.

After installation, you can load **ReverseEcologyR** into current workspace by typing or pasting the following codes:

 ```R
library("ReverseEcologyR")
 ```
**Note
For insall our package using devtools correctly from source, make sure packages `ReverseEcologyR` import in: `KEGGREST`, `mmnet` and `Biobase` which is in Bioconductor have been installed in your computer.
. 

