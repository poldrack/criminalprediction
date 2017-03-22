"""
load data created by Criminal_prediction.Rmd and make plots for paper
"""

import numpy,pandas
import seaborn
import matplotlib.pyplot as plt
plt.clf()

# load data and drop rows with nans
baseline=pandas.DataFrame.from_csv('baseline_cverr.txt',sep=' ').dropna()
dacc=pandas.DataFrame.from_csv('dacc_cverr.txt',sep=' ').dropna()
dacc['Index']=[int(i) for i in dacc.index]
baseline['Index']=[int(i) for i in baseline.index]

dacc_long=pandas.melt(dacc,id_vars='Index')
base_long=pandas.melt(baseline,id_vars='Index')
base_long['Model']='Age'
dacc_long['Model']='Age+dACC'
longdata=base_long.append(dacc_long)
longdata['variable']=[int(i.replace('V','')) for i in longdata.variable]
seaborn.tsplot(data=longdata,time='variable',condition='Model',
                value='value',unit='Index',err_style='ci_band')
plt.xlabel('Time (months)',fontsize=18)
plt.ylabel('Prediction error',fontsize=18)
plt.savefig('prediction_error.png',dpi=300)
#plt.show()

# make survival plot split by dacc
plt.clf()
surv=pandas.read_csv('km_accsplit_nonvio.txt',sep=' ')
plt.plot(surv['km_acc.time'][:48],surv['km_acc.surv'][:48])
plt.plot(surv['km_acc.time'][48:],surv['km_acc.surv'][48:],'red')
plt.axis([0,48,0,1])
plt.legend(['above dACC median','below dACC median'],fontsize=18)
plt.ylabel('Survival (proportion not rearrested)',fontsize=18)
plt.xlabel('Time (months)',fontsize=18)

plt.savefig('survival_by_dACC.png',dpi=300)
