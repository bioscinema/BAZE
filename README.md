# BAZE Package

BAZE is an R package for Bayesian tree-based selection of microbiome features with enhanced visualizations providing Graphlan annotation file or using ggtree to generate figures directly in R.

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
## Generate overall taxonomy tree plot at user-specified level based on ggtree package
### First, delete duplicate taxa
ps1 <- fix_duplicate_tax(ps)

### Convert phyloseq to treeio data and generate plot at phylum level
tr1 <- phy_to_tax(ps1)
p <- taxview(ps1, tr1, branch_thickness=0.5, levele="Phylum")
```

![ibd_select](https://github.com/bioscinema/BAZE/assets/90227639/9b225374-df76-466a-8d64-73af1cf900f6)


## Perform main analysis within BAZE algorithm (when you install package, there will be a phyloseq named "ps" in BAZE package)
```r
result <- result(ps, nburnin=10000, niter=5000, a=-9, level="Genus", response="bmi")
```

## Begin from the selection plot and effect size plot
```r
## Extract taxonomy table first
mytax <- as.data.frame(tax_table(ps))
## Generate selection probability plot
aggregate_plot <- (nburnin = 10000, niter = 5000, result, taxa_table = mytax, level="Genus")
```

![image](https://github.com/bioscinema/BAZE/assets/90227639/4284d743-5db9-48ee-9d83-f01ea087f105)

```r
## Generate effect size plot
p <- effect_plot(result,ps,nburnin = 10000,niter = 5000,mode = "meidan", level = "Genus")
p
```
![image](https://github.com/bioscinema/BAZE/assets/90227639/8862813d-2589-486a-b5ae-b3a824590013)


### Generate selection tree plot
We provide two kinds of taxonomy tree plot with showing relative abundance in different ways: node size and bar plot.

First, you need to run BAZE main algorithm at each level you want

```r
## Run BAZE at any level you want
result_genus <- BAZE::result(ps,nburnin = 15000,niter = 10000,a=-5,response = "metabolic_disease_score",level = "Genus")
result_class <- result(ps,nburnin = 15000,niter = 10000,a=-3,response = "metabolic_disease_score",level = "Class")
result_order <- result(ps,nburnin = 15000,niter = 10000,a=-4,response = "metabolic_disease_score",level = "Order")
result_family <- result(ps,nburnin = 15000,niter = 10000,a=-5,response = "metabolic_disease_score",level = "Family")

## Merge result into one data frame
eff_genus <- BAZE::effect_sign(result_genus,ps,nburnin = 15000, niter = 10000, level = "Genus")
eff_family <- BAZE::effect_sign(result_family,ps,nburnin = 15000, niter = 10000, level = "Family")
eff_order <- BAZE::effect_sign(result_order,ps,nburnin = 15000, niter = 10000, level = "Order")
eff_class <- BAZE::effect_sign(result_class,ps,nburnin = 15000, niter = 10000, level = "Class")
effect <- rbind(eff_genus,eff_family,eff_class,eff_order,eff_species)

## Generate tree plot with different function
p2 <- tax_rel(tr1,anno.data = effect,anno.depth = 6,scale_size = 0.2)
p2
```

![hiv_select](https://github.com/bioscinema/BAZE/assets/90227639/cc52c178-d7d2-4a18-a798-aa892f7343b1)

```r
## fix duplicate taxa
ps1 <- fix_duplicate_tax(ps)

##create treeio subject
tr1 <- phy_to_tax(ps1)

##create annotation data frame
effect <- data.frame(node=c("g__Streptococcus","g__Dialister"),color=c("green","blue"))

##create annotated figure
p2 <- plottax(tr1,anno.data = effect,anno.depth = 6)
p2
```
![hiv_select](https://github.com/bioscinema/BAZE/assets/90227639/5e4a514e-4ee7-4b7b-b7a7-890962ea5a9e)





