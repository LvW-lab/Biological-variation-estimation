#Packages

library(GAD)
library(outliers) 
library(tidyverse)
library(VCA)

###### Supporting functions#####################

C.test.alt <- function (object)  #Altered version of C.test applicable for the custom log.dm.ctest.func() function                     
{
  model <- deparse(substitute(object))
  by.factor <- as.factor(1:object$rank)
  n <- length(object$model[, 1])/object$rank
  k <- object$rank
  var <- tapply(object$model[, 1], rep(1:k, each = n), var)
  int <- interaction(object$model[, -1], lex.order = TRUE)
  f.int <- factor(int, levels = unique(int))
  names(var) <- levels(f.int)
  mean <- tapply(object$model[, 1], rep(1:k, each = n), mean)
  C <- max(var)/sum(var)
  group <- names(var)[which(var == max(var))]
  method <- "Cochran test of homogeneity of variances"
  alt <- paste(group)
  f <- (1/C - 1)/(k - 1)
  p <- 1 - pf(f, (n - 1) * (k - 1), (n - 1)) * k
  pval <- 1 - p
  result <- list(statistic = c(C = C), parameter = c(n = n, 
                                                     k = k), alternative = alt, p.value = pval, method = method, 
                 estimate = round(var, 4), mean = mean, var = var, data.names = model)
  class(result) <- "htest"
  return(result)
}

cochran.test.alt <- function (object, data, inlying = FALSE) #Altered version of cochran.test applicable for the custom log.dm.ctest.func() function 
{
  DNAME <- deparse(substitute(object))
  if (is.vector(object)) {
    by.factor <- as.factor(1:length(data))
    vars <- object
    names(vars) <- levels(by.factor)
    k <- length(data)
    df <- mean(data)
  }
  else {
    if (missing(data)) 
      data <- environment(object)
    bn <- as.character(attr(terms(object), "variables")[-1])
    by.factor <- as.factor(data[[bn[2]]])
    vars <- tapply(data[[bn[1]]], by.factor, var)
    names(vars) <- levels(by.factor)
    k <- nlevels(by.factor)
    df <- length(data[[bn[1]]])/k
  }
  if (inlying) {
    value <- min(vars)/sum(vars)
    group <- levels(by.factor)[which(vars == min(vars))]
    method <- "Cochran test for inlying variance"
    alt <- paste(group)
    pval <- pcochran(value, df, k)
  }
  else {
    value <- max(vars)/sum(vars)
    group <- levels(by.factor)[which(vars == max(vars))]
    method <- "Cochran test for outlying variance"
    alt <- paste(group)
    pval <- 1 - pcochran(value, df, k)
  }
  RVAL <- list(statistic = c(C = value), parameter = c(df = df, 
                                                       k = k), alternative = alt, p.value = pval, method = method, 
               estimate = vars, data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}

############ Reiterating function for stepwise exclusion of duplicate measures outliers - Normal ##############################

dm.ctest.func <- function(dm.df)
{
  # Store C test in object
  pat <- as.fixed(dm.df$patient)
  tp <- as.random(dm.df$time.point) 
  ctest.model <- C.test.alt(lm(meas ~ pat + tp%in%pat,data= dm.df))
  
  # Check for outlying variance and remove it from dm df
  if(ctest.model$p.value <= 0.05) {
    print(paste("Duplicate measures", ctest.model$alt, "has outlying variance, subject is removed from dataset"))
    excl <- strsplit(ctest.model$alt, "\\.")[[1]]
    newDF<-dm.df[!(dm.df$pat == excl[1]), ]
    dm.ctest.func(newDF)
  }
  else {
    print("No outlier detected")
    return(dm.df)
  }
}

############ Reiterating function for stepwise exclusion of duplicate measures outliers - Log ##############################

log.dm.ctest.func <- function(dm.df)
{
  # Store C test in object
  pat <- as.fixed(dm.df$patient)
  tp <- as.random(dm.df$time.point) 
  ctest.model <- C.test.alt(lm(log.meas ~ pat + tp%in%pat,data= dm.df))
  
  # Check for outlying variance and remove it from dm df
  if(ctest.model$p.value <= 0.05) {
    print(paste("Duplicate measures", ctest.model$alt, "has outlying variance, subject is removed from dataset"))
    excl <- strsplit(ctest.model$alt, "\\.")[[1]]
    newDF<-dm.df[!(dm.df$pat == excl[1]), ]
    log.dm.ctest.func(newDF)
  }
  else {
    print("No outlier detected")
    return(dm.df)
  }
}

############ Reiterating function for stepwise exclusion of duplicate measures outliers - CV ##############################

cv.dm.ctest.func <- function(dm.df)
{
  # Store C test in object
  pat <- as.fixed(dm.df$patient)
  tp <- as.random(dm.df$time.point) 
  ctest.model <- C.test.alt(lm(cv.meas ~ pat + tp%in%pat,data= dm.df))
  
  # Check for outlying variance and remove it from dm df
  if(ctest.model$p.value <= 0.05) {
    print(paste("Duplicate measures", ctest.model$alt, "has outlying variance, subject is removed from dataset"))
    excl <- strsplit(ctest.model$alt, "\\.")[[1]]
    newDF<-dm.df[!(dm.df$pat == excl[1]), ]
    cv.dm.ctest.func(newDF)
  }
  else {
    print("No outlier detected")
    return(dm.df)
  }
}

############ Reiterating function for stepwise exclusion of within-subject outliers - Normal ##############################

ws.ctest.func <- function(ws.df)
{
  # Store C test in object
  ctest.model <- cochran.test.alt(meas.mean ~ pat, ws.df)
  
  # Check for outlying variance and remove it from ws df
  if(ctest.model$p.value <= 0.05) {
    print(paste("Subject", ctest.model$alt, "has outlying variance, subject is removed from dataframe"))
    excl <- ctest.model$alt
    newDF <- ws.df[!(ws.df$pat == excl), ]
    ws.ctest.func(newDF)
  }
  else {
    print("No outlier detected")
    return(ws.df)
  }
}

############ Reiterating function for stepwise exclusion of within-subject outliers - Log ##############################

log.ws.ctest.func <- function(ws.df)
{
  # Store C test in object
  ctest.model <- cochran.test.alt(log.meas.mean ~ pat, ws.df)
  
  # Check for outlying variance and remove it from ws df
  if(ctest.model$p.value <= 0.05) {
    print(paste("Subject", ctest.model$alt, "has outlying variance, subject is removed from dataframe"))
    excl <- ctest.model$alt
    newDF <- ws.df[!(ws.df$pat == excl), ]
    log.ws.ctest.func(newDF)
  }
  else {
    print("No outlier detected")
    return(ws.df)
  }
}

############ Reiterating function for stepwise exclusion of within-subject outliers - CV ##############################

cv.ws.ctest.func <- function(ws.df)
{
  # Store C test in object
  ctest.model <- cochran.test.alt(cv.meas.mean ~ pat, ws.df)
  
  # Check for outlying variance and remove it from ws df
  if(ctest.model$p.value <= 0.05) {
    print(paste("Subject", ctest.model$alt, "has outlying variance, subject is removed from dataframe"))
    excl <- ctest.model$alt
    newDF <- ws.df[!(ws.df$pat == excl), ]
    cv.ws.ctest.func(newDF)
  }
  else {
    print("No outlier detected")
    return(ws.df)
  }
}

############ Reiterating function for stepwise exclusion of between-subject outliers - Normal ##############################

bs.reed.func <- function(bs.df)
{
  # Reorder mean subject values
  bs.df <- bs.df[order(bs.df$patmeas.mean), ]
 
  # Calculate Reeds Criterion parameters
  a <- bs.df$patmeas.mean[2] - bs.df$patmeas.mean[1]
  b <- bs.df$patmeas.mean[nrow(bs.df)] - bs.df$patmeas.mean[nrow(bs.df) - 1]
  c <- (bs.df$patmeas.mean[nrow(bs.df)] - bs.df$patmeas.mean[1])/3
  
  # Test whether the minimum or maximum values are outlying
  if(a > c){
    print("The lowest subject mean is outlying. The subject is removed from the dataset")
    newDF <- bs.df[-1, ]
    bs.reed.func(newDF)
  }
  else if(b > c) {
    print("The highest subject mean is outlying. The subject is removed from the dataset")
    newDF <- bs.df[-nrow(bs.df), ]
    bs.reed.func(newDF)
  }
  else {
    print("No outliers detected")
    return(bs.df)
  }
}

############ Reiterating function for stepwise exclusion of between-subject outliers - Log ##############################

log.bs.reed.func <- function(bs.df)
{
  # Reorder mean subject values
  bs.df <- bs.df[order(bs.df$log.patmeas.mean), ]
  
  # Calculate Reeds Criterion parameters
  a <- bs.df$log.patmeas.mean[2] - bs.df$log.patmeas.mean[1]
  b <- bs.df$log.patmeas.mean[nrow(bs.df)] - bs.df$log.patmeas.mean[nrow(bs.df) - 1]
  c <- (bs.df$log.patmeas.mean[nrow(bs.df)] - bs.df$log.patmeas.mean[1])/3
  
  # Test whether the minimum or maximum values are outlying
  if(a > c){
    print("The lowest subject mean is outlying. The subject is removed from the dataset")
    newDF <- bs.df[-1, ]
    log.bs.reed.func(newDF)
  }
  else if(b > c) {
    print("The highest subject mean is outlying. The subject is removed from the dataset")
    newDF <- bs.df[-nrow(bs.df), ]
    log.bs.reed.func(newDF)
  }
  else {
    print("No outliers detected")
    return(bs.df)
  }
}

############ Normality test for within- and between-subject values - Normal ##############################

normality.func <- function(ws.df, bs.df) 
{
  # Extract p values of Shapiro Wilk tests and list them with tapply
  norm.pat <- as.data.frame(do.call("rbind", with(ws.df, tapply(meas.mean, pat, function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))))
  norm.patmean <- shapiro.test(bs.df$patmeas.mean)
  # Check whether at least 50% of the within-subject values are normally distributed 
  if(sum(norm.pat$p.value < 0.05) > (0.5*nrow(norm.pat))) {
    return(print("Warning! <50% of within-subject values normally distributed"))
  }
  else {
    print("Data OK! >50% of within-subject values normally distributed")
    
    # Check whether subject mean values are normally distributed 
    if(norm.patmean$p.value < 0.05) {
      return(print("Warning! Means of patients are not normally distributed"))
    }
    else {
      return(print("Means of patients are normally distributed"))
    }
  }
}

############ Normality test for within- and between-subject values - Log ##############################

log.normality.func <- function(ws.df, bs.df) 
{
  # Extract p values of Shapiro Wilk tests and list them with tapply
  norm.pat <- as.data.frame(do.call("rbind", with(ws.df, tapply(log.meas.mean, pat, function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))))
  norm.patmean <- shapiro.test(bs.df$log.patmeas.mean)
  
  # Check whether at least 50% of the within-subject values are normally distributed 
  if(sum(norm.pat$p.value < 0.05) > (0.5*nrow(norm.pat))) {
    return(print("Warning! <50% of within-subject values normally distributed"))
  }
  else {
    print("Data OK! >50% of within-subject values normally distributed")
    
    # Check whether subject mean values are normally distributed 
    if(norm.patmean$p.value < 0.05) {
      return(print("Warning! Means of patients are not normally distributed"))
    }
    else {
      return(print("Means of patients are normally distributed"))
    }
  }
}

############ Normality test for within-subject values - CV ##############################

cv.normality.func <- function(ws.df) 
{
  # Extract p values of Shapiro Wilk tests and list them with tapply
  norm.pat <- as.data.frame(do.call("rbind", with(ws.df, tapply(cv.meas.mean, pat, function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))))
  
  # Check whether at least 50% of the within-subject values are normally distributed 
  if(sum(norm.pat$p.value < 0.05) > (0.5*nrow(norm.pat))) {
    return(print("Warning! <50% of within-subject values normally distributed"))
  }
  else {
    return("Data OK! >50% of within-subject values normally distributed")
    
  }
}


############ Steady-State analysis - Normal ##############################

steadystate.func <- function(ws.df)
{
  # linear regression and extraction of slope + CI
  summary.fit <- summary(lm(meas.mean ~ time.point, ws.df))
  slope <- summary.fit$coefficients[2,1]
  lower <- slope - (1.96 * summary.fit$coefficients[2,2])
  upper <- slope + (1.96 * summary.fit$coefficients[2,2])
  
  # Check whether CI of slope contains 0
  if(0 < upper & 0 > lower) {
    return(print("Group is in steady state"))
  }
  else {
    return(print("Group is not in steady state, apply inverse regression"))
  }
}


############ Steady-State analysis - Log ##############################

log.steadystate.func <- function(ws.df)
{
  # linear regression and extraction of slope + CI
  summary.fit <- summary(lm(log.meas.mean ~ time.point, ws.df))
  slope <- summary.fit$coefficients[2,1]
  lower <- slope - (1.96 * summary.fit$coefficients[2,2])
  upper <- slope + (1.96 * summary.fit$coefficients[2,2])
  
  # Check whether CI of slope contains 0
  if(0 < upper & 0 > lower) {
    return(print("Group is in steady state"))
  }
  else {
    return(print("Group is not in steady state, apply inverse regression"))
  }
}

############ Steady-State - CV ##############################

cv.steadystate.func <- function(ws.df)
{
  # linear regression and extraction of slope + CI
  summary.fit <- summary(lm(cv.meas.mean ~ time.point, ws.df))
  slope <- summary.fit$coefficients[2,1]
  lower <- slope - (1.96 * summary.fit$coefficients[2,2])
  upper <- slope + (1.96 * summary.fit$coefficients[2,2])
  
  # Check whether CI of slope contains 0
  if(0 < upper & 0 > lower) {
    return(print("Group is in steady state"))
  }
  else {
    return(print("Group is not in steady state, apply inverse regression"))
  }
}

############ Parameter calculation - Normal ##############################

parameter.calc.func <- function(dm.finaldf, ws.finaldf, bs.finaldf) 
{
  # Variance and CVA calculation from DM C test
  pat <- as.fixed(dm.finaldf$patient)
  tp <- as.random(dm.finaldf$time.point)
  ctest.dm.model <- C.test(lm(meas ~ pat + tp%in%pat,data = dm.finaldf))
  dmvar <-  data.frame("var" = ctest.dm.model$var)
  s2a <- mean(dmvar$var)
  cva <- (sqrt(exp(s2a)-1)) * 100
  cva.upper95 <- cva + (1.96 * (cva / (sqrt(2 * nrow(dm.finaldf)))))
  cva.lower95 <- cva - (1.96 * (cva / (sqrt(2 * nrow(dm.finaldf)))))
  
  # Variance and CVI calculation from WS C test
  ctest.ws.model <- cochran.test(meas.mean ~ pat, ws.finaldf)
  wsvar <- data.frame("var" = ctest.ws.model$estimate)
  s2i <- (mean(wsvar$var) - (0.5 * (s2a)))
  cvi <- (sqrt(exp(s2i)-1)) * 100
  cvi.upper95 <- cvi + (1.96 * (cvi / (sqrt(2 * nrow(ws.finaldf)))))
  cvi.lower95 <- cvi - (1.96 * (cvi / (sqrt(2 * nrow(ws.finaldf)))))
  
  # Variance and CVG calculation with nested ANOVA
  lmer.fit <- summary(lmer(meas ~ (1|pat/time.point), dm.finaldf))
  lmer.df <- as.data.frame(lmer.fit$varcor)
  s2g <- lmer.df$sdcor[2]^2
  cvg <- (sqrt(exp(s2g)-1)) * 100
  cvg.upper95 <- cvg + (1.96 * (cvg / (sqrt(2 * nrow(bs.finaldf)))))
  cvg.lower95 <- cvg - (1.96 * (cvg / (sqrt(2 * nrow(bs.finaldf)))))
  
  # index of individuality and RCV calculation
  ii <- (sqrt((cva^2) + (cvi^2)))/cvg
  rcv <- 2.77 * (cva^2 + cvi^2)^0.5
  
  # Descriptives calculation
  mean.conc <- exp(mean(bs.finaldf$patmeas.mean))
  min.conc <- exp(min(bs.finaldf$patmeas.mean))
  max.conc <- exp(max(bs.finaldf$patmeas.mean))
  n <- nrow(bs.finaldf)
  
  # Desirable APS calculation
  imprecision <- 0.5 * cvi
  bias <- 0.25 * (sqrt((cvi^2) + (cvg ^ 2)))
  toterr <- (1.65 * cva) + bias
 
  # Store all results in object and print results
  result <- list(n = n,
                 mean = mean.conc, 
                 minimum = min.conc, 
                 maximum = max.conc, 
                 cva = cva, 
                 cva.lower = cva.lower95, 
                 cva.upper = cva.upper95,
                 cvi = cvi, 
                 cvi.lower = cvi.lower95, 
                 cvi.upper = cvi.upper95,
                 cvg = cvg, 
                 cvg.lower = cvg.lower95, 
                 cvg.upper = cvg.upper95,
                 ii = ii,
                 rcv = rcv,
                 imprecision = imprecision,
                 bias = bias,
                 total.error = toterr)
  as.data.frame(result)
  print(result)
  return(result)
}

############ Parameter calculation - Log ##############################

log.parameter.calc.func <- function(dm.finaldf, ws.finaldf, bs.finaldf) 
{
  # Variance and CVA calculation from DM C test
  pat <- as.fixed(dm.finaldf$patient)
  tp <- as.random(dm.finaldf$time.point)
  ctest.dm.model <- C.test(lm(log.meas ~ pat + tp%in%pat,data = dm.finaldf))
  dmvar <-  data.frame("var" = ctest.dm.model$var)
  s2a <- mean(dmvar$var)
  cva <- (sqrt(exp(s2a)-1)) * 100
  cva.upper95 <- cva + (1.96 * (cva / (sqrt(2 * nrow(dm.finaldf)))))
  cva.lower95 <- cva - (1.96 * (cva / (sqrt(2 * nrow(dm.finaldf)))))
  
  # Variance and CVI calculation from WS C test
  ctest.ws.model <- cochran.test(log.meas.mean ~ pat, ws.finaldf)
  wsvar <- data.frame("var" = ctest.ws.model$estimate)
  s2i <- (mean(wsvar$var) - (0.5 * (s2a)))
  cvi <- (sqrt(exp(s2i)-1)) * 100
  cvi.upper95 <- cvi + (1.96 * (cvi / (sqrt(2 * nrow(ws.finaldf)))))
  cvi.lower95 <- cvi - (1.96 * (cvi / (sqrt(2 * nrow(ws.finaldf)))))
  
  # Variance and CVG calculation with nested ANOVA
  lmer.fit <- summary(lmer(log.meas ~ (1|pat/time.point), dm.finaldf))
  lmer.df <- as.data.frame(lmer.fit$varcor)
  s2g <- lmer.df$sdcor[2]^2
  cvg <- (sqrt(exp(s2g)-1)) * 100
  cvg.upper95 <- cvg + (1.96 * (cvg / (sqrt(2 * nrow(bs.finaldf)))))
  cvg.lower95 <- cvg - (1.96 * (cvg / (sqrt(2 * nrow(bs.finaldf)))))
  
  # index of individuality and RCV calculation
  ii <- (sqrt((cva^2) + (cvi^2)))/cvg
  rcv <- 2.77 * (cva^2 + cvi^2)^0.5
  
  # Descriptives calculation
  mean.conc <- exp(mean(bs.finaldf$log.patmeas.mean))
  min.conc <- exp(min(bs.finaldf$log.patmeas.mean))
  max.conc <- exp(max(bs.finaldf$log.patmeas.mean))
  n <- nrow(bs.finaldf)
  
  # Desirable APS calculation
  imprecision <- 0.5 * cvi
  bias <- 0.25 * (sqrt((cvi^2) + (cvg ^ 2)))
  toterr <- (1.65 * cva) + bias
  
  # Store all results in object and print results
  result <- list(n = n,
                 mean = mean.conc, 
                 minimum = min.conc, 
                 maximum = max.conc, 
                 cva = cva, 
                 cva.lower = cva.lower95, 
                 cva.upper = cva.upper95,
                 cvi = cvi, 
                 cvi.lower = cvi.lower95, 
                 cvi.upper = cvi.upper95,
                 cvg = cvg, 
                 cvg.lower = cvg.lower95, 
                 cvg.upper = cvg.upper95,
                 ii = ii,
                 rcv = rcv,
                 imprecision = imprecision,
                 bias = bias,
                 total.error = toterr)
  as.data.frame(result)
  print(result)
  return(result)
}

############ Parameter calculation - CV norm ##############################

cv.norm.parameter.calc.func <- function(dm.finaldf, ws.finaldf, bs.finaldf) 
{
  # Variance and CVA calculation from DM C test
  pat <- as.fixed(dm.finaldf$patient)
  tp <- as.random(dm.finaldf$time.point)
  ctest.dm.model <- C.test(lm(cv.meas ~ pat + tp%in%pat,data = dm.finaldf))
  dmvar <-  data.frame("var" = ctest.dm.model$var)
  s2a <- mean(dmvar$var)
  cva <- (sqrt(exp(s2a)-1)) * 100
  cva.upper95 <- cva + (1.96 * (cva / (sqrt(2 * nrow(dm.finaldf)))))
  cva.lower95 <- cva - (1.96 * (cva / (sqrt(2 * nrow(dm.finaldf)))))
  
  # Variance and CVI calculation from WS C test
  ctest.ws.model <- cochran.test(cv.meas.mean ~ pat, ws.finaldf)
  wsvar <- data.frame("var" = ctest.ws.model$estimate)
  s2i <- (mean(wsvar$var) - (0.5 * (s2a)))
  cvi <- (sqrt(exp(s2i)-1)) * 100
  cvi.upper95 <- cvi + (1.96 * (cvi / (sqrt(2 * nrow(ws.finaldf)))))
  cvi.lower95 <- cvi - (1.96 * (cvi / (sqrt(2 * nrow(ws.finaldf)))))
  
  # Variance and CVG calculation with nested ANOVA
  lmer.fit <- summary(lmer(meas ~ (1|pat/time.point), dm.finaldf))
  lmer.df <- as.data.frame(lmer.fit$varcor)
  s2g <- lmer.df$sdcor[2]^2
  cvg <- (sqrt(exp(s2g)-1)) * 100
  cvg.upper95 <- cvg + (1.96 * (cvg / (sqrt(2 * nrow(bs.finaldf)))))
  cvg.lower95 <- cvg - (1.96 * (cvg / (sqrt(2 * nrow(bs.finaldf)))))
  
  # index of individuality and RCV calculation
  ii <- (sqrt((cva^2) + (cvi^2)))/cvg
  rcv <- 2.77 * (cva^2 + cvi^2)^0.5
  
  # Descriptives calculation
  mean.conc <- mean(bs.finaldf$patmeas.mean)
  min.conc <- min(bs.finaldf$patmeas.mean)
  max.conc <- max(bs.finaldf$patmeas.mean)
  n <- nrow(bs.finaldf)
  
  # Desirable APS calculation
  imprecision <- 0.5 * cvi
  bias <- 0.25 * (sqrt((cvi^2) + (cvg ^ 2)))
  toterr <- (1.65 * cva) + bias
  
  # Store all results in object and print results
  result <- list(n = n,
                 mean = mean.conc, 
                 minimum = min.conc, 
                 maximum = max.conc, 
                 cva = cva, 
                 cva.lower = cva.lower95, 
                 cva.upper = cva.upper95,
                 cvi = cvi, 
                 cvi.lower = cvi.lower95, 
                 cvi.upper = cvi.upper95,
                 cvg = cvg, 
                 cvg.lower = cvg.lower95, 
                 cvg.upper = cvg.upper95,
                 ii = ii,
                 rcv = rcv,
                 imprecision = imprecision,
                 bias = bias,
                 total.error = toterr)
  as.data.frame(result)
  print(result)
  return(result)
}


############ Parameter calculation - CV Log ##############################

cv.log.parameter.calc.func <- function(dm.finaldf, ws.finaldf, bs.finaldf) 
{
  # Variance and CVA calculation from DM C test
  pat <- as.fixed(dm.finaldf$patient)
  tp <- as.random(dm.finaldf$time.point)
  ctest.dm.model <- C.test(lm(cv.meas ~ pat + tp%in%pat,data = dm.finaldf))
  dmvar <-  data.frame("var" = ctest.dm.model$var)
  s2a <- mean(dmvar$var)
  cva <- (sqrt(exp(s2a)-1)) * 100
  cva.upper95 <- cva + (1.96 * (cva / (sqrt(2 * nrow(dm.finaldf)))))
  cva.lower95 <- cva - (1.96 * (cva / (sqrt(2 * nrow(dm.finaldf)))))
  
  # Variance and CVI calculation from WS C test
  ctest.ws.model <- cochran.test(cv.meas.mean ~ pat, ws.finaldf)
  wsvar <- data.frame("var" = ctest.ws.model$estimate)
  s2i <- (mean(wsvar$var) - (0.5 * (s2a)))
  cvi <- (sqrt(exp(s2i)-1)) * 100
  cvi.upper95 <- cvi + (1.96 * (cvi / (sqrt(2 * nrow(ws.finaldf)))))
  cvi.lower95 <- cvi - (1.96 * (cvi / (sqrt(2 * nrow(ws.finaldf)))))
  
  # Variance and CVG calculation with nested ANOVA
  lmer.fit <- summary(lmer(log.meas ~ (1|pat/time.point), dm.finaldf))
  lmer.df <- as.data.frame(lmer.fit$varcor)
  s2g <- lmer.df$sdcor[2]^2
  cvg <- (sqrt(exp(s2g)-1)) * 100
  cvg.upper95 <- cvg + (1.96 * (cvg / (sqrt(2 * nrow(bs.finaldf)))))
  cvg.lower95 <- cvg - (1.96 * (cvg / (sqrt(2 * nrow(bs.finaldf)))))
  
  # index of individuality and RCV calculation
  ii <- (sqrt((cva^2) + (cvi^2)))/cvg
  rcv <- 2.77 * (cva^2 + cvi^2)^0.5
  
  # Descriptives calculation
  mean.conc <- exp(mean(bs.finaldf$log.patmeas.mean))
  min.conc <- exp(min(bs.finaldf$log.patmeas.mean))
  max.conc <- exp(max(bs.finaldf$log.patmeas.mean))
  n <- nrow(bs.finaldf)
  
  #Desirable APS calculation
  imprecision <- 0.5 * cvi
  bias <- 0.25 * (sqrt((cvi^2) + (cvg ^ 2)))
  toterr <- (1.65 * cva) + bias
  
  # Store all results in object and print results
  result <- list(n = n,
                 mean = mean.conc, 
                 minimum = min.conc, 
                 maximum = max.conc, 
                 cva = cva, 
                 cva.lower = cva.lower95, 
                 cva.upper = cva.upper95,
                 cvi = cvi, 
                 cvi.lower = cvi.lower95, 
                 cvi.upper = cvi.upper95,
                 cvg = cvg, 
                 cvg.lower = cvg.lower95, 
                 cvg.upper = cvg.upper95,
                 ii = ii,
                 rcv = rcv,
                 imprecision = imprecision,
                 bias = bias,
                 total.error = toterr)
  as.data.frame(result)
  print(result)
  return(result)
}