# re-analysis of Aharoni et al. (2013) PNAS paper
# using cross-validation to examine out-of-sample generalization

library(survival)
library(rms)
library(ROCR)
library(pec)

# load data
alldata=read.table('aharoni_data_fixed.txt',header=TRUE, na.strings=".")

use_IQ=FALSE
# optionally, throw out subjects without IQ measure
if (use_IQ == TRUE) {
	alldata=alldata[!is.na(alldata$Wais_IQ),]
	alldata$Wais_IQ_centered=alldata$Wais_IQ - mean(alldata$Wais_IQ)
}

# get overall survival curve
s=survival::Surv(time=alldata$MinMonthsDV_noPVs,event=alldata$AnyChargeSinceScanExclPVs,type='right')


# plot the Kaplan-Meier curve if you want to
plot_km=FALSE
if ( plot_km == TRUE) {
	km=survfit(s~1)
	plot(km)
}

# set up models to test
if (use_IQ == TRUE) {
models=list("Age"=coxph(Surv(MinMonthsDV_noPVs,AnyChargeSinceScanExclPVs)~releaseAge_centered, 
                       data=alldata, y=TRUE), 
            "dACC"=coxph(Surv(MinMonthsDV_noPVs,AnyChargeSinceScanExclPVs)~dACC_centered, 
                       data=alldata, y=TRUE), 
         	"Age.dACC"=coxph(Surv(MinMonthsDV_noPVs,AnyChargeSinceScanExclPVs)~releaseAge_centered+ dACC_centered, data=alldata, y=TRUE),                   
            "IQ"=coxph(Surv(MinMonthsDV_noPVs,AnyChargeSinceScanExclPVs)~Wais_IQ_centered, 
                       data=alldata, y=TRUE))
  } else {
models=list("Age"=coxph(Surv(MinMonthsDV_noPVs,AnyChargeSinceScanExclPVs)~releaseAge_centered, 
                       data=alldata, y=TRUE), 
            "dACC"=coxph(Surv(MinMonthsDV_noPVs,AnyChargeSinceScanExclPVs)~dACC_centered, 
                       data=alldata, y=TRUE), 
         	"Age.dACC"=coxph(Surv(MinMonthsDV_noPVs,AnyChargeSinceScanExclPVs)~releaseAge_centered+ dACC_centered, data=alldata, y=TRUE))
}
	
nrep=100   # number of crossvalidation runs to average

pred_err=pec(models,formula=Surv(MinMonthsDV_noPVs,AnyChargeSinceScanExclPVs)~1, 
                       data=alldata,
                        splitMethod="cv10", B=nrep,
                       	verbose=FALSE,
                       	keep.index=TRUE,
                       	keep.matrix=TRUE)

print (pred_err)

# plot the data for the dACC model

plot(pred_err, models=c(1,3), col=c("black", "green"),lwd=2)
for (i in 1:nrep){
	lines(pred_err$time, pred_err$CrossValErrMat[[i]]$Reference, hex(RGB(0.2,0.2,0.2)),type='s',lwd=0.1)
	lines(pred_err$time, pred_err$CrossValErrMat[[i]]$dACC, col=hex(RGB(0.0,0.7,0,0)),  type='s',lwd=0.1)
}


