# re-analysis of Aharoni et al. (2013) PNAS paper
# using cross-validation to examine out-of-sample generalization

library(survival)

# load data
alldata=read.table('aharoni_data_fixed.txt',header=TRUE)

# get survival curve
s=Surv(time=alldata$MinMonthsDV_noPVs,event=alldata$AnyChargeSinceScanExclPVs,type='right')

# plot the Kaplan-Meier curve if you want to
#km=survfit(s~alldata$dACC_14mm_split)
# plot(km)

# compute Cox proportional hazards regression model
# using dACC activation as the predictor

cph=coxph(s~dACC_centered,data=alldata)
pred_cph=predict(cph,type='risk')

# print the in-sample accuracy at a range of risk thresholds
show_thresholds=0
if (show_thresholds==1) {
  for (i in c(-15:15)/100 + 1) {
	print(c(i,sum((pred_cph>i)==alldata$AnyChargeSinceScanExclPVs)/96))
	}

}

# set to 1 to shuffle data for empirical null analysis
shuffle_data=0

# use more runs for shuffled data

if (shuffle_data==1){
	nruns=100
	} else {
	nruns=100
	}
predacc=rep(0,nruns)
pred=rep(0,96)

folds=rep(c(1,2,3,4,5,6,7,8),12)

for (r in 1:nruns) {
  folds=folds[sample(96)]
  for (fold in 1:8) {
	if (shuffle_data==1){
		alldata$AnyChargeSinceScanExclPVs=alldata$AnyChargeSinceScanExclPVs[sample(96)]
		}
	train=which(folds!=fold)
	test=which(folds==fold)
	traindata=alldata[train,]
	testdata=data=alldata[test,]
	s_train=Surv(time=traindata$MinMonthsDV_noPVs,event=traindata$AnyChargeSinceScanExclPVs,type='right')
	cph_train=coxph(s_train~dACC_centered,data=traindata)
	pred[test]=predict(cph_train,newdata=testdata,type='risk')
	}
	# compute predictive accuracy relative to threshold
	predacc[r]=sum((pred>1.02)==alldata$AnyChargeSinceScanExclPVs)/96
}

print(c('mean_predacc:',mean(predacc)))

if (shuffle_data==1){
	print(c('95th pctile:',sort(predacc)[nruns*0.95]))
	}
