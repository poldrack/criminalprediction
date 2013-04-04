predacc <- function(time,status,predictor){

alldata=data.frame(time,status,predictor)
nobs=length(time)
results=c()
s=survival::Surv(time=alldata$time,event=alldata$status,type='right')
c=cph(s~predictor)
#print(c)
results$insample_cph=c

for (shuffle_data in c(0)) {
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
			alldata$status=alldata$status[sample(nobs)]
			}
	  # create balanced folds
	  nreps=ceil(nobs/8)
	  folds=rep(c(1:8),nreps)
	  for (foldset in 1:nreps) {
			folds[((foldset-1)*8 +1):(foldset*8)]=c(1:8)[sample(8)]
			}
	  folds_rand=rep(nobs)
	  folds_rand[alldata$status==1]=folds[1:length(which(alldata$status==1))]
	  folds_rand[alldata$status==0]=folds[(length(which(alldata$status==1))+1):nobs]
	  
	  
	  for (fold in 1:8) {
		train=which(folds_rand!=fold)
		test=which(folds_rand==fold)
		traindata=alldata[train,]
		testdata=alldata[test,]
		s_train=survival::Surv(time=traindata$time,event=traindata$status,type='right')
		s_test=survival::Surv(time=testdata$time,event=testdata$status,type='right')
		cph_train=psm(s_train~predictor,data=traindata,dist='exponential')
		psp=predictSurvProb(cph_train,newdata=testdata,times=testdata$time,train.data=traindata)
		psp_pred=rep(0,length(testdata[,1]))
		
		for (x in 1:length(testdata[,1])) {
			#if (psp[x]<testdata$time[x]) {psp_pred[x]=1}
			psp_pred[x]=psp[x,x]
			}
		
		#print(psp_pred)
		#print(testdata$status)
		p=prediction(psp_pred,testdata$status)
		
		perf <- performance(p,"tpr","fpr")
		auc.tmp <- performance(p,"auc")
		auc_gen[r]= as.numeric(auc.tmp@y.values)
		prbe=performance(p,'prbe')
		# if there are multiple breakeven points, just take the last one
		# this can be an issue for the shuffled datasets
		#print(p)
		if (length(prbe@x.values[[1]])>1) {
		
			prbe_cutoff[r]=prbe@x.values[[1]][length(prbe@x.values[[1]])]
			prbe_ppv[r]=prbe@y.values[[1]][length(prbe@y.values[[1]])]
			}
		else {
			prbe_cutoff[r]=prbe@x.values[[1]]
			prbe_ppv[r]=prbe@y.values[[1]]
			
					}
		}
	}


	print(c('mean_auc:',mean(auc_gen)))
	if (shuffle_data==1){
		print(c('95th pctile of AUC null dist:',sort(auc_gen)[nruns*0.95]))
		results$auc_gen_rand=auc_gen
		} else {
			results$auc_gen=auc_gen
		
			}
	
	print(c('mean_ppv at PRBE:',mean(prbe_ppv)))
	
	if (shuffle_data==1){
		print(c('95th pctile of PPV null dist:',sort(prbe_ppv)[nruns*0.95]))
	results$prbe_ppv_rand=prbe_ppv
		} else {
		results$prbe_ppv=prbe_ppv
			}
		

}
return(results)
}