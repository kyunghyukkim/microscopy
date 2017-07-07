# -*- coding: utf-8 -*-

# Creates contour plot with scatter plot
# Kiri Choi, 2017

from __future__ import print_function, division

import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage
import seaborn
import h5py

clistpathlist = ['Path to uclist.mat']

def idx(string, char):
  for key, x in enumerate(string):
    if char in x:
      return key

#%% Reading data
datalist = []
data3dlist = []
deflist = []
def3Dlist = []
exclist = []

for i in range(len(clistpathlist)):
    row_datalist = []
    row_data3Dlist = []
    f = h5py.File(clistpathlist[i],'r') 
    data = f.get(u'data') 
    data = np.array(data) # For converting to numpy array
    datalist.append(data)
    data3d = f.get(u'data3D')
    data3d = np.array(data3d) # For converting to numpy array
    data3dlist.append(data3d)
    try:
        exc = f.get(u'idExclude')
        exc = np.array(exc)
        for i in range(len(exc)):
            exc[i] = exc[i] - 1
        exclist.append(exc)
    except:
        exclist.append([])
    for column in f['def']:
        row_data = []
        for row_number in range(len(column)):            
            row_data = ''.join(map(unichr, f[column[row_number]][:]))
        row_datalist.append(row_data.split(': ')[1])
    deflist.append(row_datalist)
    for column in f['def3d']:
        row_data3D = []
        for row_number in range(len(column)):            
            row_data3D = ''.join(map(unichr, f[column[row_number]][:]))
        row_data3Dlist.append(row_data3D.split(': ')[1])
    def3Dlist.append(row_data3Dlist)
    f.close()
    
#%% Exclude
if len(exclist) > 0:
    for i in range(len(exclist)):
        datalist[i] = np.delete(datalist[i], exclist[i], 1)
        for j in range(len(exclist[i])):
            data3dlist[i][j][:,np.int(exclist[i][j])] = np.nan
            
#%% Populate list
f1 = []
f2 = []

for i in range(len(datalist)):
    f1.append(datalist[i][idx(deflist[i], 'fluor1 mean'),:])
    f2.append(datalist[i][idx(deflist[i], 'fluor2 mean'),:])

mcherry = [item for sublist in f1 for item in sublist]
venus = [item for sublist in f2 for item in sublist]

#%% Plotting
seaborn.set_style("white")
sb = seaborn.hls_palette(8, l=.3, s=.8)

scale=1

fig = plt.gcf()
fig.set_size_inches(8*scale,5*scale)
ax = fig.add_subplot(111)

counts,xbins,ybins = np.histogram2d(venus, mcherry, bins=100, 
                                range = ([90, 160],[0, 20000]), normed=True)
counts = ndimage.gaussian_filter(counts, sigma=3.0, order=0)

plt.contourf(np.transpose(counts),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],
             levels = np.linspace(np.max(counts)/50, np.max(counts), 6), cmap=plt.cm.Reds, alpha=1.)
plt.scatter(venus, mcherry, marker = '.', color = sb[5], alpha = 0.3)

plt.xticks(fontsize = 15*scale)
plt.yticks(fontsize = 15*scale)
plt.xlabel("Venus (AU)", fontsize=20*scale)
plt.ylabel("mCherry (AU)", fontsize=20*scale)
plt.axis([95, 160, 0, 16000])

plt.show()
    

