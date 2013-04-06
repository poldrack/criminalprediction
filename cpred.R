# re-analysis of Aharoni et al. (2013) PNAS paper
# using cross-validation to examine out-of-sample generalization

library(survival)
library(rms)
library(ROCR)
library(pec)

# load data
alldata=read.table('aharoni_data_fixed.txt',header=TRUE, na.strings=".")
alldata=alldata[!is.na(alldata$Wais_IQ),]
alldata$Wais_IQ_centered=alldata$Wais_IQ - mean(alldata$Wais_IQ)

# get survival curve
s=survival::Surv(time=alldata$MinMonthsDV_noPVs,event=alldata$AnyChargeSinceScanExclPVs,type='right')


# plot the Kaplan-Meier curve if you want to
km=survfit(s~1)
#plot(km)

# compute Cox proportional hazards regression model
# using dACC activation as the predictor

c_exp=psm(s~ dACC_peak_centered,data=alldata,dist='exponential',x=TRUE,y=TRUE)
c_exp_null=psm(s~ 1,data=alldata,dist='exponential',x=TRUE,y=TRUE)

c_weib= psm(s~ dACC_peak_centered,data=alldata,dist='weibull')
c_gaussian= psm(s~ dACC_peak_centered,data=alldata,dist='gaussian')
c_logistic= psm(s~ dACC_peak_centered,data=alldata,dist='logistic')
c_lognormal= psm(s~ dACC_peak_centered,data=alldata,dist='lognormal')
c_cox=coxph(s~ dACC_peak_centered,data=alldata)

AIC=c()
AIC['exp']=extractAIC(c_exp)[2]
AIC['exp_null']=extractAIC(c_exp_null)[2]
AIC['weib']=extractAIC(c_weib)[2]
AIC['gaussian']=extractAIC(c_gaussian)[2]
AIC['logistic']=extractAIC(c_logistic)[2]
AIC['lognormal']=extractAIC(c_lognormal)[2]
AIC['cox']=extractAIC(c_cox)[2]

print(AIC)

v=validate(c_exp,method='.632')
print(v)
modeltype='cox'
if (modeltype=='exp') {
models=list("Age"=psm(Surv(MinMonthsDV_noPVs,AnyChargeSinceScanExclPVs)~releaseAge_centered, 
                       data=alldata, y=TRUE,dist='exponential'), 
                    "dACC"=psm(Surv(MinMonthsDV_noPVs,AnyChargeSinceScanExclPVs)~dACC_centered, 
                       data=alldata, y=TRUE,dist='exponential'), 
                    "Age.dACC"=psm(Surv(MinMonthsDV_noPVs,AnyChargeSinceScanExclPVs)~releaseAge_centered+ dACC_centered, 
                       data=alldata, y=TRUE,dist='exponential'))
} else
{
models=list("Age"=coxph(Surv(MinMonthsDV_noPVs,AnyChargeSinceScanExclPVs)~releaseAge_centered, 
                       data=alldata, y=TRUE), 
                    "dACC"=coxph(Surv(MinMonthsDV_noPVs,AnyChargeSinceScanExclPVs)~dACC_centered, 
                       data=alldata, y=TRUE), 
                    "Age.dACC"=coxph(Surv(MinMonthsDV_noPVs,AnyChargeSinceScanExclPVs)~releaseAge_centered+ dACC_centered, 
                       data=alldata, y=TRUE),
                       
"IQ"=coxph(Surv(MinMonthsDV_noPVs,AnyChargeSinceScanExclPVs)~Wais_IQ_centered, 
                       data=alldata, y=TRUE))
	
	}
	
nrep=100

if (1==1) {
  pred_err=pec(models,formula=Surv(MinMonthsDV_noPVs,AnyChargeSinceScanExclPVs)~1, 
                       data=alldata,
                        splitMethod="cv10", B=nrep,
                       	verbose=TRUE,
                       	keep.index=TRUE,
                       	keep.matrix=TRUE)
}

#par(mfrow=c(3,1))
plot(pred_err)

#looks to me that the age+dACC model does have improved accuracy over the model with no covariates, but only over the time window around 20-28 months
library(colorspace)
#plot the average lines with the 100 separate boot strap estimates for each
plot(pred_err, models=c(1,4), col=c("black", "blue"),lwd=2)
for (i in 1:nrep){
	lines(pred_err$time, pred_err$CrossValErrMat[[i]]$Reference, col=hex(RGB(0.2,0.2,0.2)), type='s',lwd=0.1)
	lines(pred_err$time, pred_err$CrossValErrMat[[i]]$Age.dACC, col=hex(RGB(0.0,0.0,0.7)), type='s',lwd=0.1)
}

#of course I think the model of interest was dACC alone and there's a lot of overlap in those
plot(pred_err, models=c(1,3), col=c("black", "green"),lwd=2)
for (i in 1:nrep){
	lines(pred_err$time, pred_err$CrossValErrMat[[i]]$Reference, hex(RGB(0.2,0.2,0.2)),type='s',lwd=0.1)
	lines(pred_err$time, pred_err$CrossValErrMat[[i]]$dACC, col=hex(RGB(0.0,0.7,0,0)),  type='s',lwd=0.1)
}

# results from exponential psm
#Integrated Brier score (crps):
#
#     IBS[0;time=48]
#[1,]          0.216
#[2,]          0.214
#[3,]          0.209
#[4,]          0.205

