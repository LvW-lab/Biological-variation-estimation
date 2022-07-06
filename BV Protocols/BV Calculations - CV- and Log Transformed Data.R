# Packages
library(GAD)
library(outliers) 
library(tidyverse)
library(VCA)
library(lme4)
library(data.table)

# Loading external dataset 
dm.df <- read.table(file = "dataset_example.txt",
                    header = TRUE,
                    sep = "\t",
                    stringsAsFactors = FALSE)

# Convert dataframe to data.table for CV transformation
dm.df <- data.table(dm.df, key = "pat")

# Transform values into CV values
# All values are divided by the corresponding subject mean 
dm.df <- dm.df[ , cv.meas := meas /  
                  mean(meas), "pat"]

#Log Transformation of all values
dm.df <- dm.df %>%
  mutate(log.meas = log(meas))

# Assign pat as factor variable - necessary for c test
dm.df$pat <- as.character(dm.df$pat)

# Duplicate measures outlier identification and removal
# Cocnran C test
dm.ex.df <- cv.dm.ctest.func(dm.df)

# Create ws.df by calculating means of duplicate measures
ws.df <- dm.ex.df %>%
  group_by(pat, time.point)%>%
  summarise(
    meas.mean = mean(meas),
    log.meas.mean = mean(log.meas),
    cv.meas.mean = mean(cv.meas),
  )

# Within-subject outlier identification and removal
# Cochran C test
ws.ex.df <- cv.ws.ctest.func(ws.df)

# Create bs.df by calculating means of within-subject values
bs.df <- ws.ex.df %>%
  group_by(pat) %>%
  summarise(
    patmeas.mean = mean(meas.mean),
    log.patmeas.mean = mean(log.meas.mean),
    cv.patmeas.mean = mean(cv.meas.mean)) 

# Between-subject outlier identification and removal
# Reeds Criterion Test
bs.ex.df <- bs.reed.func(bs.df)

# Create final dfs by removing all outlying patients that were identified
dm.finaldf <- subset(dm.ex.df, pat %in% bs.ex.df$pat)
ws.finaldf <- subset(ws.ex.df, pat %in% bs.ex.df$pat)
bs.finaldf <- bs.ex.df

# Perform steady state analysis 
# Linear regression on all values
cv.steadystate.func(ws.finaldf)

# Perform normality test for CV values
# Shapiro-Wilk test on within- and between-subject values
cv.normality.func(ws.finaldf)

# Perform normality test for normal values
# Shapiro-Wilk test on within- and between-subject values
normality.func(ws.finaldf, bs.finaldf)

# Calculation of all BV estimates
# Results are saved in external results csv file
results <- cv.norm.parameter.calc.func(dm.finaldf, ws.finaldf, bs.finaldf)

results <- as.data.frame(results)

write.table(results, file = "results.csv", sep = ";", row.names = FALSE)
