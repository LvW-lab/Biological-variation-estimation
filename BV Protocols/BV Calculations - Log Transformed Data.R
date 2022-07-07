#Packages

library(GAD)
library(outliers) 
library(tidyverse)
library(VCA)
library(lme4)

# Loading external dataset
dm.df <- read.table(file = "dataset_example.txt",
                    header = TRUE,
                    sep = "\t",
                    stringsAsFactors = FALSE)

# Log Transformation of all values
dm.df <- dm.df %>%
  mutate(log.meas = log(meas))

# Assign pat as character variable - necessary for c test
dm.df$pat <- as.character(dm.df$pat)

# Duplicate measures outlier identification and removal
# Cocnran C test
dm.ex.df <- log.dm.ctest.func(dm.df)

# Create ws.df by calculating means of duplicate measures
ws.df <- dm.ex.df %>%
  group_by(pat, time.point)%>%
  summarise(
    meas.mean = mean(meas),
    log.meas.mean = mean(log.meas))

# Within-subject outlier identification and removal
# Cochran C test
ws.ex.df <- log.ws.ctest.func(ws.df)

# Create bs.df by calculating means of within-subject values
bs.df <- ws.ex.df %>%
  group_by(pat) %>%
  summarise(
    patmeas.mean = mean(meas.mean),
    log.patmeas.mean = mean(log.meas.mean)) 

# Between-subject outlier identification and removal
# Reeds Criterion Test
# Note that at least two values are required to perform this test!
bs.ex.df <- log.bs.reed.func(bs.df)

# Create final dfs by removing all outlying patients that were identified
dm.finaldf <- subset(dm.ex.df, pat %in% bs.ex.df$pat)
ws.finaldf <- subset(ws.ex.df, pat %in% bs.ex.df$pat)
bs.finaldf <- bs.ex.df

# Perform steady state analysis 
# Linear regression on all values
log.steadystate.func(ws.finaldf)

# Perform normality test
# Shapiro-Wilk test on within- and between-subject values
# Note that at least three values are required to perform this test!
log.normality.func(ws.finaldf, bs.finaldf)

# Calculation of all BV estimates
# Results are saved in external results csv file
results <- log.parameter.calc.func(dm.finaldf, ws.finaldf, bs.finaldf)

results <- as.data.frame(results)

write.table(results, file = "results.csv", sep = ";", row.names = FALSE)
