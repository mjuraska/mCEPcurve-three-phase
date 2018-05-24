# this code computes Monte-Carlo point estimates of the TE(s1) curve in Figure 1 using the MLE-TE estimator and saves them 
# in a file for evaluation of bias and MSE

rm(list=ls(all=TRUE))

library(RCurl)
library(rstudioapi)

# source 'TEcurveMC_myFunctions.R'
script <- getURL("https://raw.githubusercontent.com/mjuraska/mCEPcurve-three-phase/master/TEcurveMC_myFunctions.R", ssl.verifypeer = FALSE)
eval(parse(text=script))

# the total sample size
n <- 5000

# the indicator of whether marker distributions are truncated
truncateMarker <- TRUE

# grid of biomarker values on which the performance of the estimators is evaluated
if (truncateMarker){
  # calculate the values of the probit risk model coefficients
  beta2 <- -0.02
  # the marginal probability of infection in the placebo group set to 0.1
  beta0 <- getBeta0(beta2, meanS0=2, varS0=1, 0.1, scenario="A2", threshold=1.5)
  beta1start <- startBeta1(beta0, 0.75)
  beta3 <- getBeta3(beta0, beta1start, beta2, meanS1=3, varS1=1, 0.25)
  # the marginal probability of infection in the vaccine group set to 0.05
  beta1 <- getBeta1(beta0, beta2, beta3, meanS1=3, varS1=1, 0.05, scenario="A2", threshold=1.5)
  beta <- c(beta0, beta1, beta2, beta3)
  
  # the grid of marker values
  s1 <- seq(1.5, qnorm(0.95, mean=3, sd=1), by=0.05)
} else {
  # calculate the values of the probit risk model coefficients
  beta2 <- -0.02
  # the marginal probability of infection in the placebo group set to 0.1
  beta0 <- getBeta0(beta2, meanS0=2, varS0=1, 0.1)
  beta1start <- startBeta1(beta0, 0.75)
  beta3 <- getBeta3(beta0, beta1start, beta2, meanS1=3, varS1=1, 0.25)
  # the marginal probability of infection in the vaccine group set to 0.05
  beta1 <- getBeta1(beta0, beta2, beta3, meanS1=3, varS1=1, 0.05)
  beta <- c(beta0, beta1, beta2, beta3)
  
  # the grid of marker values
  q <- qnorm(c(0.05,0.95), mean=3, sd=1)
  s1 <- seq(q[1], q[2], by=0.05)
}

# the number of MC iterations
nMC <- 1000

# a sampling probability for sampling a Phase 2/3 subcohort with biomarker measurements
pi <- c(0.1, 0.25, 0.5)

# save the below point estimates of the TE(s1) curve obtained from each MC iteration in 'currentDir'
currentDir <- dirname(getSourceEditorContext()$path)

for (i in 1:length(pi)){
  TE.MCestimates <- sapply(1:nMC, function(seed, s1grid, n, beta, pi, truncateMarker){ suppressWarnings(getPestVE(s1grid, n, beta, pi, truncateMarker, seed)) }, s1grid=s1, n=n, beta=beta, pi=pi[i], truncateMarker=truncateMarker)
  # save the matrix 'TE.MCestimates' as an .RData file
  save(TE.MCestimates, file=file.path(currentDir, paste0("pTE_MCestimates_ns.df=3_nMC=",nMC,"_N=",n,"_truncateMarker=",truncateMarker,"_pi=",pi[i],".RData")))
}
