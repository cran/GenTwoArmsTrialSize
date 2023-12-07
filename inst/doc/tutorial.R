## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(GenTwoArmsTrialSize)

## -----------------------------------------------------------------------------
getSizeMean(design="parallel", test="equivalence", alpha=0.05, 
            beta=0.20, sigma=0.10, k=1, delta=0.05, TTE=0.01, 
            rho=c(0.05, 0.07), r=0.1)

## -----------------------------------------------------------------------------
getSizeProp(design="crossover", test="superiority", alpha=0.05,
            beta=0.20, varsigma=c(0.5, 0.5), k=1, seqnumber=2,
            delta=0.10, TTE=0, rho=c(0.05, 0.07), r=0.1)

## -----------------------------------------------------------------------------
getSizeTTE(design="parallel", test="equality", alpha=0.05, beta=0.20,
           varlambda=c(1,2), k=1, ttotal=3, taccrual=1, 
           gamma=0.00001, delta=1, rho=c(0.05, 0.07), r=0.1)

## -----------------------------------------------------------------------------
getSizeOrd(design="parallel", test="equality", alpha=0.05, beta=0.10,
           varcatprob= list(c(0.2,0.5,0.2,0.1),
                            c(0.378,0.472,0.106,0.044)), 
           k=1, theta=0.877, delta=0, rho=c(0.05, 0.07), r=0.1)

