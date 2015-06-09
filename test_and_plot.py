# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import h5py
import pymc3
import numpy as np
import seaborn as sns
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib as mpl

data = h5py.File('/triton/becs/scratch/braindata/shared/2pnhyperEMO/BioPacTesting/data_processed.mat')
val = np.array(data.get('val'))
aro = np.array(data.get('aro'))
br = np.array(data.get('br'))
hr = np.array(data.get('hr'))

plt.grid(True, color='w', linestyle='-', linewidth=2)
ax1.yaxis.grid(False)





biopacdf = []
for subj in range(val.shape[1]):
    biopacdf.append(pd.DataFrame(dict(valence=val[:,subj],
                                      arousal=aro[:,subj],
                                      heartrate=stats.mstats.zscore(hr[:,subj]),
                                      breath=stats.mstats.zscore(br[:,subj]),
                                      subject=['subject{0}'.format(subj)] * val.shape[0],
                                      timevec=[i for i in range(val.shape[0])]), dtype=np.float))
biopacdf = pd.concat(biopacdf)

ax1 = sns.tsplot(biopacdf, time="timevec", unit="subject", value="valence", ci=95)
ax1.set_xticks(np.arange(0,1015,29))
ax1.set_xlabel("time (TR)");
ax2 = sns.tsplot(biopacdf, time="timevec", unit="subject", value="heartrate", ci=95)

df=pd.DataFrame({'valence':val, 'heartbeat freq':stats.mstats.zscore(hr,axis=0)})
df.plot()
df=pd.DataFrame({'valence':val[:,1], 'heartbeat freq':stats.mstats.zscore(hr[:,1])})
sns.set(style="ticks")
sns.jointplot("valence", "heartbeat freq", df, kind="hex", color="#4CB391")


sns.set(style="darkgrid")
color = sns.color_palette()[2]
ax = sns.jointplot("valence", "heartrate", data=biopacdf, kind="reg",
                  color=color, size=7)
                  
# PYMC3 GLM
                  