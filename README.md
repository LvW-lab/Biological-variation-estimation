# Biological-variation-estimation
An R script for the estimation of biological variation.

Description: The scripts allow for calculation of BV parameters. The formulas and general script are also published in van Winden et al.(Clin Chim Acta, 2021).

Instructions: The formulas should be loaded first, then the general BV script. The BV plot is optional.

Data structure: 
The first column  “pat” describes the patients and is numbered in ascending order, whereas the second column “time.point” indicates 
the within-subject samples that are analyzed. 
It is not required to balance the number of samples per subject thus allowing for varying number of within-subject samples. 
The third column “rep” represents the two replicates that are measured to generate the analytical BV in the data. 
The last column “meas” lists the concentrations that have been measured by a specific assay. 
All information in the dataset should be inserted as numeric data and should be saved as a tab separated text file.
