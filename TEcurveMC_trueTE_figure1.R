# this code plots the true TE(s1) curves in Figure 1 assumed by the data generating mechanism in the simulation study

rm(list=ls(all=TRUE))

library(RCurl)
library(rstudioapi)

# source 'TEcurveMC_myFunctions.R'
script <- getURL("https://raw.githubusercontent.com/mjuraska/mCEPcurve-three-phase/master/TEcurveMC_myFunctions.R", ssl.verifypeer = FALSE)
eval(parse(text=script))

# regression coefficients in the probit model for the association of (S(0),S(1)) with the endpoint indicator in each treatment group
beta2 <- -0.02
# the marginal probability of infection in the placebo group set to 0.1
print(beta0 <- getBeta0(beta2, meanS0=2, varS0=1, 0.1, scenario="A2", threshold=1.5))
print(beta1start <- startBeta1(beta0, 0.75))
print(beta3 <- getBeta3(beta0, beta1start, beta2, meanS1=3, varS1=1, 0.25))
# the marginal probability of infection in the vaccine group set to 0.05
print(beta1 <- getBeta1(beta0, beta2, beta3, meanS1=3, varS1=1, 0.05, scenario="A2", threshold=1.5))

q <- qnorm(c(0.1,0.9), mean=3, sd=1)
risk0(q[2], beta0, beta2, 2, 1, 3, 1, 0.7, scenario="A2", threshold=1.5)/risk0(q[1], beta0, beta2, 2, 1, 3, 1, 0.7, scenario="A2", threshold=1.5)
pnorm(beta0+beta1+beta3*q[2])/pnorm(beta0+beta1+beta3*q[1])

# plotting parameters
cexLab <- 1.4
cexAxis <- 1.3

# saves the below .pdf file in 'currentDir'
currentDir <- dirname(getSourceEditorContext()$path)
pdf(file.path(currentDir,"trueTEcurves.pdf"), width=6, height=6)

xLim <- c(1.5, qnorm(0.95, mean=3, sd=1))
s <- seq(xLim[1], xLim[2], length.out=200)
par(mar=c(5,5.2,1.5,0.8), cex.lab=cexLab, cex.axis=cexAxis, las=1)
plot(s, trueTE(s, beta0, beta1, beta2, beta3, 2, 1, 3, 1, 0.7, scenario="A2", threshold=1.5), type="n", lwd=3, xlab=expression(s[1]), ylab="", xlim=xLim, 
     ylim=0:1, yaxt="n")
mtext(expression(paste("True ",TE(s[1]), sep="")), side=2, line=3.7, cex=cexLab, las=0)
axis(side=2, at=seq(0,1,by=0.1), cex.axis=cexAxis)

# the true TE(s1) curve used only for assessing size of the tests of H01 and H02
lines(s, rep(0.5, length(s)), lwd=3, lty="dashed")

# the true TE(s1) curve used throughout the simulation study except for assessing size of the tests of H01 and H02
lines(s, trueTE(s, beta0, beta1, beta2, beta3, 2, 1, 3, 1, 0.7, scenario="A2", threshold=1.5), lwd=3)

# these regression coefficients produce the below dot-dashed TE(s1) curve used only for assessing power of the test of H04
beta2.trial2 <- -0.02
# the marginal probability of infection in the placebo group set to 0.1
print(beta0.trial2 <- getBeta0(beta2.trial2, meanS0=2, varS0=1, 0.1, scenario="A2", threshold=1.5))
print(beta1start.trial2 <- startBeta1(beta0.trial2, 0.75))
print(beta3.trial2 <- getBeta3(beta0.trial2, beta1start.trial2, beta2.trial2, meanS1=3, varS1=1, 0.5))
# the marginal probability of infection in the vaccine group set to 0.05
print(beta1.trial2 <- getBeta1(beta0.trial2, beta2.trial2, beta3.trial2, meanS1=3, varS1=1, 0.05, scenario="A2", threshold=1.5))

q <- qnorm(c(0.1,0.9), mean=3, sd=1)
risk0(q[2], beta0.trial2, beta2.trial2, 2, 1, 3, 1, 0.7, scenario="A2", threshold=1.5)/risk0(q[1], beta0.trial2, beta2.trial2, 2, 1, 3, 1, 0.7, scenario="A2", threshold=1.5)
pnorm(beta0.trial2+beta1.trial2+beta3.trial2*q[2])/pnorm(beta0.trial2+beta1.trial2+beta3.trial2*q[1])

# the true TE(s1) curve used only for assessing power of the test of H04
lines(s, trueTE(s, beta0.trial2, beta1.trial2, beta2.trial2, beta3.trial2, 2, 1, 3, 1, 0.7, scenario="A2", threshold=1.5), lwd=3, lty="dotdash")

dev.off()
