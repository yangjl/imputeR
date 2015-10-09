# imputeR

GBS data normally have high error rate for heterozygote sites.
**ImputeR** is a package to infer the most likely genotypes using raw GBS data.

## Install

Install [devtools](https://github.com/hadley/devtools) first, and then use `devtools` to install `imputeR` from github.

```R
devtools::install_github("hadley/devtools")
library(devtools)
install_github("yangjl/imputeR")
library(imputeR)
```

## Documentation

A vignette can be found [here](https://github.com/yangjl/imputeR/blob/master/vignettes/imputeR-vignette.pdf).
Documented functions are listed as below. Their usage information can by found by typing `?function_name` or `help(function_name)`.

 - impute_mom

