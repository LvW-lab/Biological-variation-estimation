# Biological-variation-estimation README

## Introduction
An R script for the estimation of biological variation. The formulas and general script are also published in van Winden et al.(Clin Chim Acta, 2021).

## Instructions

### Data structure

An example dataset is provided in the "Synthetic Dataset" folder.
We recommended using this data structure when calculating BV using our R scripts. 

The first column  “pat” describes the patients and is numbered in ascending order.
The second column “time.point” indicates the number (ascending order) of within-subject samples. 
It is not required to balance the number of samples per subject thus allowing for varying number of within-subject samples. 
The third column “rep” represents the replicates that are measured. 
Each within-subject sample contains two replicate measures.
The final column “meas” lists the values of measurement by a specific assay.
All values, except for the column headings, should be numeric.
All information in the dataset should be inserted as numeric data and should be saved as a tab separated text file.

### Custom functions

Please note that before running BV calculation and BV plot R script the Functions R script should be run as they won't function without the custom functions.

### BV calculation protocol scripts

It is recommended to first check the dataset on normality.
If within-subject values or subject means are not normally distributed, transformation scripts can be considered.
Various protocols are provided generating BV estimates from untransformed, log transformed and CV transformed data. 
CV transformation is an especially powerful method for calculation of analytical and within person variation. 
However, all subject means are transformed to a value of 1. 
Therefore, untransformed or log transformed data voor between-subject variation calculation should be used. 
Protocols for both combinations are provided (CV- and untransformed or CV- and log transformed).
After running the BV calculations script, a BV plot of the data can be generated using the BV Plot Generation script.

### Example of BV Calculations output (Untransformed Data)

```js
# Loading external dataset 
> dm.df <- read.table(file = "dataset_example.txt", 
+                     header = TRUE,
+                     sep = "\t",
+                     stringsAsFactors = FALSE)
> 
> # Assign pat as factor variable - necessary for c test
> dm.df$pat <- as.character(dm.df$pat)
> 
> # Duplicate measures outlier identification and removal
> # Cocnran C test
> dm.ex.df <- dm.ctest.func(dm.df)
[1] "Duplicate measures 23.2 has outlying variance, subject is removed from dataset"
[1] "No outlier detected"
> 
> # Create ws.df by calculating means of duplicate measures
> ws.df <- dm.ex.df %>%
+   group_by(pat, time.point)%>%
+   summarise(
+     meas.mean = mean(meas))
`summarise()` has grouped output by 'pat'. You can override using the `.groups` argument.
> 
> # Within-subject outlier identification and removal
> # Cochran C test
> ws.ex.df <- ws.ctest.func(ws.df)
[1] "Subject 5 has outlying variance, subject is removed from dataframe"
[1] "No outlier detected"
> 
> # Create bs.df by calculating means of within-subject values
> bs.df <- ws.ex.df %>%
+   group_by(pat) %>%
+   summarise(
+     patmeas.mean = mean(meas.mean)) 
> 
> # Between-subject outlier identification and removal
> # Reeds Criterion Test
> bs.ex.df <- bs.reed.func(bs.df)
[1] "The highest subject mean is outlying. The subject is removed from the dataset"
[1] "No outliers detected"
> 
> # Create final dfs by removing all outlying patients that were identified
> dm.finaldf <- subset(dm.ex.df, pat %in% bs.ex.df$pat)
> ws.finaldf <- subset(ws.ex.df, pat %in% bs.ex.df$pat)
> bs.finaldf <- bs.ex.df
> 
> # Perform steady state analysis 
> # Linear regression on all values
> steadystate.func(ws.finaldf)
[1] "Group is in steady state"
> 
> # Perform normality test
> # Shapiro-Wilk test on within- and between-subject values
> normality.func(ws.finaldf, bs.finaldf)
[1] "Data OK! >50% of within-subject values normally distributed"
[1] "Means of patients are normally distributed"
> 
> # Calculation of all BV estimates
> # Results are saved in external results csv file
> results <- parameter.calc.func(dm.finaldf, ws.finaldf, bs.finaldf)
$n
[1] 27

$mean
[1] 2.931196

$minimum
[1] 2.163081

$maximum
[1] 4.698221

$cva
[1] 4.870302

$cva.lower
[1] 4.514552

$cva.upper
[1] 5.226053

$cvi
[1] 10.73668

$cvi.lower
[1] 9.627566

$cvi.upper
[1] 11.84579

$cvg
[1] 19.65791

$cvg.lower
[1] 14.41471

$cvg.upper
[1] 24.90111

$ii
[1] 0.5997413

$rcv
[1] 32.65735

$imprecision
[1] 5.368338

$bias
[1] 5.599718

$total.error
[1] 13.63572

> 
> results <- as.data.frame(results)
> 
> write.table(results, file = "results.csv", sep = ";", row.names = FALSE)

```

The code was written in R (Version 4.2.1) and RStudio (Version 2022.02.3). 

The packages and versions used for our scripts:

- GAD (version 1.1.1)
- outliers (version 0.15)
- tidyverse (version 1.3.1)
-	VCA (version 1.4.3)
-	lme4 (version 1.1-29)
-	data.table (version 1.14.2)
