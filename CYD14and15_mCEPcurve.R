# this code conducts the analysis of the dengue vaccine trial data reported in the application section
# it's a pooled analysis of 9-16 year-olds in CYD14+15 who are at risk and free of VCD at month 13
# point estimates of the risk curve in each treatment arm are computed to allow estimation of two mCEP curves: log RR and additive difference
# the proposed NP estimator is employed

rm(list=ls(all=TRUE))

library(RCurl)
library(rstudioapi)

currentDir <- dirname(getSourceEditorContext()$path)

# source 'CYD14and15_myFunctions.R'
script <- getURL("https://raw.githubusercontent.com/mjuraska/mCEPcurve-three-phase/master/CYD14and15_myFunctions.R", ssl.verifypeer = FALSE)
eval(parse(text=script))

# target set: ITT set at-risk at month 13 with no prior infection
# the below code assumes that 'data' is the data frame for analysis

riskCurve(data, markerName="AUC", saveFile="riskCurves_CYD14and15_9to16_AUCMB_hingePoint.RData", saveDir=currentDir)

# compute bootstrap estimates of the risk curve in each treatment arm
# caveat: a parallelized version of 'bRiskCurve' is required for nMC=1000 and nBoot=500
bRiskCurve(data, markerName="AUC", iter=1000, saveFile="bRiskCurves_CYD14and15_9to16_AUCMB_hingePoint.RData", saveDir=currentDir)

# load the output from the above call of 'riskCurve' (a list named 'oList')
load(file.path(currentDir,"riskCurves_CYD14and15_9to16_AUCMB_hingePoint.RData"))

# load the output from the parallelized version of the above call of 'bRiskCurve' (a list named 'result')
load(file.path(currentDir,"bRiskCurves_CYD14and15_9to16_AUCMB_hingePoint.RData"))

plaRiskCurveBootEst <- sapply(result, "[[", "plaRiskCurveBootEst")
txRiskCurveBootEst <- sapply(result, "[[", "txRiskCurveBootEst")

# FIGURE 5 (left column)
load(file.path(currentDir, "CYD14and15_PSNestimator.Rdata"))

pdf(file.path(currentDir, "logRRcurve_CYD14and15_9to16_AUCMB_hingePoint_v4.pdf"), width=6, height=5)
plotLogRRcurve(oList$markerVals, oList$plaRiskCurve, oList$txRiskCurve, plaRiskCurveBootEst, txRiskCurveBootEst, 
               title="Proposed NP Estimator", hingePoint=round(10^min(oList$cpointP, oList$cpointV, na.rm=TRUE),0), 
               yLim=range(c(low.ci.LRR, high.ci.LRR, low.cb.LRR, high.cb.LRR), na.rm=TRUE))
dev.off()

pdf(file.path(currentDir, "riskDiffCurve_CYD14and15_9to16_AUCMB_hingePoint_v4.pdf"), width=6, height=5)
plotRiskDiffCurve(oList$markerVals, oList$plaRiskCurve, oList$txRiskCurve, plaRiskCurveBootEst, txRiskCurveBootEst, 
                  title="Proposed NP Estimator", plotLegend=FALSE, yLim=range(c(low.ci.RD, high.ci.RD, low.cb.RD, high.cb.RD), na.rm=TRUE))
dev.off()
# END OF FIGURE 5 (left column)

# TESTS OF H01 and H02
# compute the overall risk of VCD in each treatment arm
fit <- glm(ofstatus_m13 ~ VACC, data=data, family=binomial)
prob <- predict(fit, newdata=data.frame(VACC=0:1), type="response")

testConstancy2(oList$plaRiskCurve, oList$txRiskCurve, plaRiskCurveBootEst, txRiskCurveBootEst, plaRiskOverall=prob[1], txRiskOverall=prob[2], 
               MCEPcontrast="multiplicativeTE", null="H01")

testConstancy2(oList$plaRiskCurve, oList$txRiskCurve, plaRiskCurveBootEst, txRiskCurveBootEst, plaRiskOverall=prob[1], txRiskOverall=prob[2], 
               MCEPcontrast="additiveTE", null="H01")

testConstancy2(oList$plaRiskCurve, oList$txRiskCurve, plaRiskCurveBootEst, txRiskCurveBootEst, tMCEPconstantNull=0, 
               MCEPcontrast="multiplicativeTE", null="H02", S1=oList$markerVals, limS1=c(0.6,log10(57)))

testConstancy2(oList$plaRiskCurve, oList$txRiskCurve, plaRiskCurveBootEst, txRiskCurveBootEst, tMCEPconstantNull=0, 
               MCEPcontrast="additiveTE", null="H02", S1=oList$markerVals, limS1=c(0.6,log10(57)))
# END OF TESTS OF H01 and H02
