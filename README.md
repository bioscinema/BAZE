# BAZE Package

BAZE is an R package for Bayesian tree-based selection of microbiome features with enhanced visualizations.

## Installation

You can install the released version of BAZE from GitHub with:

```r
if (!require(devtools)) {
  install.packages("devtools")
}
devtools::install_github("bioscinema/BAZE")
```

## basic usage

Here's a basic example of how to use the functions from the BAZE package:

```r
## Load package
library(BAZE)

## Perform main analysis within BAZE algorithm (when you install package, there will be a phyloseq named "ps" in BAZE package)
result <- result(ps, nburnin=10000, niter=5000, a=-9, level="Genus", response="bmi")

## Begin from the selection plot and effect size plot
### Extract taxonomy table first
mytax <- as.data.frame(tax_table(ps))
### Generate selection probability plot
aggregate_plot <- (nburnin = 10000, niter = 5000, result, taxa_table = mytax, level="Genus")
```
![image](https://github.com/bioscinema/BAZE/assets/90227639/4284d743-5db9-48ee-9d83-f01ea087f105)

```r
## Generate effect size plot
p <- effect_plot(result,ps,nburnin = 10000,niter = 5000,mode = "meidan", level = "Genus")
p
```
![image](https://github.com/bioscinema/BAZE/assets/90227639/8862813d-2589-486a-b5ae-b3a824590013)

```r
## Generate annotation files for GraPhlAn
### Generate taxonomy tree file 
create_tax(ps,"test_tax.txt")

### Generate overall annotation file
create_tax_annot(ps,"test_annot.txt")

![test](https://github.com/bioscinema/BAZE/assets/90227639/d15e1c3e-b202-43e1-b782-ca70c9b8415e)


### Generate annotation file with selected taxa
create_tax_selected(ps,nburnin = 10000,niter = 5000,result = result,annotation_file = "test_select.txt")





