# for each Monte-Carlo iteration, this code computes indicators of CI coverage, and indicators of rejection of the null hypothesis 
# in tests of H^1_0, H^2_0, and H^4_0 under validity of both the null and an alternative hypothesis
# this code considers the MLE-TE estimator

rm(list=ls(all=TRUE))

library(RCurl)
library(rstudioapi)

currentDir <- dirname(getSourceEditorContext()$path)

# source 'TEcurveMC_myFunctions.R'
script <- getURL("https://raw.githubusercontent.com/mjuraska/mCEPcurve-three-phase/master/TEcurveMC_myFunctions.R", ssl.verifypeer = FALSE)
eval(parse(text=script))

# the total sample size
n <- 5000

# the indicator of whether marker distributions are truncated
truncateMarker <- TRUE

# the number of MC iterations
nMC <- 1000

# the number of bootstrap iterations within each MC iteration
nBoot <- 500

# a sampling probability for sampling a Phase 2/3 subcohort with biomarker measurements
pi <- c(0.1, 0.25, 0.5)

# the grid of marker values
s1 <- seq(1.5, qnorm(0.95, mean=3, sd=1), by=0.05)

# betas for coverage, power of (T1), and size of (T3)
# these betas are also used in one population for power of (T3)
beta <- getBeta(beta2=-0.02, meanS0=2, varS0=1, meanS1=3, varS1=1, probP=0.1, probV=0.05, rr=0.25, rr0=0.75, scenario="A2", threshold=1.5)

# VE(s1) curve evaluated on the grid
trueVEcurve <- trueTE(s1, beta[1], beta[2], beta[3], beta[4], 2, 1, 3, 1, 0.7, scenario="A2", threshold=1.5)

# betas for size of (T1)
# infection probability in either treatment group does not depend on S(0) or S(1), i.e., beta2 = beta3 = 0
betaSizeT1 <- getBeta(beta2=0, beta3=0, meanS0=2, varS0=1, meanS1=3, varS1=1, probP=0.1, probV=0.05, rr=1, rr0=0.5, scenario="A2", threshold=1.5)

# betas in the other population for power of (T3)
betaPowerT3 <- getBeta(beta2=-0.02, meanS0=2, varS0=1, meanS1=3, varS1=1, probP=0.1, probV=0.05, rr=0.5, rr0=0.75, scenario="A2", threshold=1.5)

for (i in 1:length(pi)){
  inferenceVE <- lapply(1:nMC, function(seed, s1grid, trueVEcurve, n, beta, betaSizeT1, betaPowerT3, pi, truncateMarker, nBoot){ 
    getPinferenceVE2(s1grid, trueVEcurve, n, beta, betaSizeT1, betaPowerT3, pi, truncateMarker, seed, nBoot) }, s1grid=s1, trueVEcurve=trueVEcurve, n=n, beta=beta, 
    betaSizeT1=betaSizeT1, betaPowerT3=betaPowerT3, pi=pi[i], truncateMarker=truncateMarker, nBoot=nBoot)
  
  save(inferenceVE, file=file.path(currentDir, paste0("pTE_coverSizePower_RoyBose_ns.df=3_nMC=",nMC,"_N=",n,"_nBoot=",nBoot,"_truncateMarker=",truncateMarker,"_pi=",pi[i],".RData")))
}
