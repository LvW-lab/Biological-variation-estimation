#BIOLOGICAL VARIATION PLOT

#Packages
library(GAD)
library(outliers) 
library(tidyverse)
library(VCA)

# Create df with the means, min and max per patient per group
descrip.df <- ws.finaldf %>%
  group_by(pat) %>%
  summarise(max = max(meas.mean),
            min = min(meas.mean),
            mean = mean(meas.mean))

# Renumber pats from 1 to nrow
descrip.df$pat <- 1:nrow(descrip.df)

descrip.df <- as.data.frame(descrip.df)

# Turn the pat variable from a numeric into a factor variable 
# and reorder the factor according to the mean values with fct_reorder
descrip.df$pat <- as.factor(descrip.df$pat)

descrip.df <- descrip.df %>%
  mutate(pat = fct_reorder(pat, mean)) 

# ggplot the means with the min/max error bars 
# and facet_wrap for treatment (free_x to use different axes) 
bv.plot <- descrip.df %>%
  ggplot(aes(x = pat,y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = min, ymax = max))

# Create final plot 
bv.plot + 
  theme_classic() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.length.x= unit (-1,"mm"),
        axis.ticks.x.top = element_line(1),
        panel.spacing = unit(0, "lines"),
        panel.border = element_rect(size = 0.25, fill = NA)) +
  labs(
    #title= "[Insert Plot Title Here]",     #possible title 
    y = "[Analyte] (SI units)", 
    x = "Subjects")
