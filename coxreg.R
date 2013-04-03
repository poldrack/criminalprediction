# re-analysis of Aharoni et al. (2013) PNAS paper
# using cross-validation to examine out-of-sample generalization

library(survival)
library(rms)
library(ROCR)

# load data
alldata=read.table('aharoni_data_fixed.txt',header=TRUE,stringsAsFactors=FALSE)
alldata$mPFC_14mm_beta=alldata$mPFC_14mm_beta-mean(alldata$mPFC_14mm_beta)
alldata$CE_mean = alldata$CE_mean  - mean(alldata$CE_mean)
# get survival curve
s=survival::Surv(time=alldata$MinMonthsDV_noPVs,event=alldata$AnyChargeSinceScanExclPVs,type='right')


# plot the Kaplan-Meier curve if you want to
#km=survfit(s~alldata$dACC_14mm_split)
# plot(km)

# compute Cox proportional hazards regression model
# using dACC activation as the predictor

c=cph(s~dACC_peak_centered,data=alldata,surv=TRUE)
#v=val.surv(c,S=s)
psp=predictSurvProb(c,times=alldata$MinMonthsDV_noPVs,newdata=alldata)
psp_pred=rep(0,96)

for (x in 1:96) {
	psp_pred[x]=psp[x,x]
	}
	

p=prediction(psp_pred,alldata$AnyChargeSinceScanExclPVs)
perf <- performance(p,"tpr","fpr")
auc.tmp <- performance(p,"auc"); auc <- as.numeric(auc.tmp@y.values)
print(c('In-sample AUC:',auc))

# print the in-sample accuracy at a range of risk thresholds
prbe=performance(p,'prbe')
print(c('PPV at PRBE point',prbe@y.values[[1]]))
npv=performance(p,'npv')
cutoff=which(npv@x.values[[1]]==prbe@x.values[[1]])[[1]]

print(c('NPV at PRBE point',npv@y.values[[1]][cutoff]))

# set to 1 to shuffle data for empirical null analysis
#shuffle_data=1

# use more runs for shuffled data

for (shuffle_data in c(0,1)) {
	
if (shuffle_data==1){
	print('shuffling labels for empirical null analysis')
	nruns=100
	} else {
		
	print('Generating crossvalidated predictive accuracy')
	nruns=10
	}
	
predacc=rep(0,nruns)
auc_gen=rep(nruns)
prbe_cutoff=rep(nruns)
prbe_ppv=rep(nruns)

for (r in 1:nruns) {
		
  if (shuffle_data==1){
		alldata$AnyChargeSinceScanExclPVs=alldata$AnyChargeSinceScanExclPVs[sample(96)]
		}
  # create balanced folds
  folds=rep(c(1:8),12)
  for (foldset in 1:12) {
		folds[((foldset-1)*8 +1):(foldset*8)]=c(1:8)[sample(8)]
		}
  folds_rand=rep(0,96)
  folds_rand[alldata$AnyChargeSinceScanExclPVs==1]=folds[1:length(which(alldata$AnyChargeSinceScanExclPVs==1))]
  folds_rand[alldata$AnyChargeSinceScanExclPVs==0]=folds[(length(which(alldata$AnyChargeSinceScanExclPVs==1))+1):96]
  
  
  for (fold in 1:8) {
	train=which(folds_rand!=fold)
	test=which(folds_rand==fold)
	traindata=alldata[train,]
	testdata=alldata[test,]
	s_train=survival::Surv(time=traindata$MinMonthsDV_noPVs,event=traindata$AnyChargeSinceScanExclPVs,type='right')
	s_test=survival::Surv(time=testdata$MinMonthsDV_noPVs,event=testdata$AnyChargeSinceScanExclPVs,type='right')
	cph_train=cph(s_train~dACC_peak_centered ,data=traindata,surv=TRUE)
	psp=predictSurvProb(cph_train,times=alldata$MinMonthsDV_noPVs,newdata=testdata)
	psp_pred=rep(0,length(testdata[,1]))

	for (x in 1:length(testdata[,1])) {
		psp_pred[x]=psp[x,x]
		}
	

	p=prediction(psp_pred,alldata$AnyChargeSinceScanExclPVs[test])
	
	perf <- performance(p,"tpr","fpr")
	auc.tmp <- performance(p,"auc")
	auc_gen[r]= as.numeric(auc.tmp@y.values)
	prbe=performance(p,'prbe')
	# if there are multiple breakeven points, just take the last one
	# this can be an issue for the shuffled datasets
	if (length(prbe@x.values[[1]])>1) {
	
		prbe_cutoff[r]=prbe@x.values[[1]][length(prbe@x.values[[1]])]
		prbe_ppv[r]=prbe@y.values[[1]][length(prbe@y.values[[1]])]
		}
	else {
		prbe_cutoff[r]=prbe@x.values[[1]]
		prbe_ppv[r]=prbe@y.values[[1]]
		
				}
	}
	# compute best predictive balanced accuracy 
	
	#predacc[r]=max((p@tp[[1]]/p@n.pos[[1]] + p@tn[[1]]/p@n.neg[[1]])/2,na.rm=TRUE)
}

print(c('mean_auc:',mean(auc_gen)))
if (shuffle_data==1){
	print(c('95th pctile of AUC null dist:',sort(auc_gen)[nruns*0.95]))
	}

print(c('mean_ppv at PRBE:',mean(prbe_ppv)))

if (shuffle_data==1){
	print(c('95th pctile of PPV null dist:',sort(prbe_ppv)[nruns*0.95]))
	}

	
# make ROC
}
