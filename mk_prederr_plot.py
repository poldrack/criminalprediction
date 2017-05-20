"""
load data created by Criminal_prediction.Rmd and make plots for paper
"""

import numpy,pandas
import seaborn
import matplotlib.pyplot as plt
plt.clf()

# load data and drop rows with nans
cverr=pandas.read_csv('cverr_data.txt',index_col=0,sep=' ')
plt.plot(cverr.time.values,cverr.baseline.values)
plt.plot(cverr.time.values,cverr.dacc.values)
plt.xlabel('Time (months)',fontsize=18)
plt.ylabel('Prediction error',fontsize=18)
plt.legend(['Baseline (age)','Age + dACC'])
plt.savefig('prediction_error.png',dpi=300)
plt.axis([0,44,0,1])
#plt.show()

# make survival plot split by dacc
plt.clf()
surv=pandas.read_csv('km_accsplit.txt',sep=' ')
plt.plot(surv['km_acc.time'][:48],surv['km_acc.surv'][:48])
plt.plot(surv['km_acc.time'][48:],surv['km_acc.surv'][48:],'red')
plt.axis([0,44,0,1])
plt.legend(['above dACC median','below dACC median'],fontsize=18)
plt.ylabel('Survival (proportion not rearrested)',fontsize=18)
plt.xlabel('Time (months)',fontsize=18)

plt.savefig('survival_by_dACC.png',dpi=300)
