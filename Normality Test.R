# Packages
library(tidyverse)

# Loading external dataset 
dm.df <- read.table(file = "dataset_example.txt",
                    header = TRUE,
                    sep = "\t",
                    stringsAsFactors = FALSE)

# Create ws.df by calculating means of duplicate measures
ws.df <- dm.df %>%
  group_by(pat, time.point)%>%
  summarise(
    meas.mean = mean(meas))

# Create bs.df by calculating means of within-subject values
bs.df <- ws.df %>%
  group_by(pat) %>%
  summarise(
    patmeas.mean = mean(meas.mean)) 

# Perform normality test for normal values
# Shapiro-Wilk test on within- and between-subject values
normality.func(ws.df, bs.df)

#If normality cannot be assumed consider log cv and/or transformation
