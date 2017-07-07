# -*- coding: utf-8 -*-

# Creates multiple scatter plots with different phase-shift values
# Kiri Choi, 2017

from __future__ import print_function, division

import os
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
def3dlist = []
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
    def3dlist.append(row_data3Dlist)
    f.close()
    
#%% Exclude
if len(exclist) > 0:
    for i in range(len(exclist)):
        #datalist[i] = np.delete(datalist[i], exclist[i], 1)
        for j in range(len(exclist[i])):
            datalist[i][:,np.int(exclist[i][j])] = np.nan
            data3dlist[i][j][:,np.int(exclist[i][j])] = np.nan
            
#%% Populate list
mcherry = []
venus = []

for i in range(len(clistpathlist)):
    mcherry.append(data3dlist[i][:,idx(def3dlist[i], 'fluor1 mean')].T)
    venus.append(data3dlist[i][:,idx(def3dlist[i], 'fluor2 mean')].T),

#%% Check original cells
cell_o_list = []

for i in range(len(clistpathlist)):
    cell_o_list.append((np.isnan(data3dlist[i][:,idx(def3dlist[i], 
                                    'fluor1 mean')][0,:]) == False).nonzero())
    
print(cell_o_list)

#%% Track lineage
exp_no = 0
cell_no = 0
line_no = 0
gen_no = 8
timediff = 75 # Delay amount
timescale = 5 # Measurement interval
numframe = timediff/timescale

gen_list_all = []
lin_list_all = []

for k in range(len(cell_o_list)):
    gen_list_temp1 = []
    lin_list_temp1 = []
    for l in range(len(cell_o_list[k][0])):
        d_idx1 = idx(deflist[k], 'daughter1 ID')
        d_idx2 = idx(deflist[k], 'daughter2 ID')
        m_idx = idx(deflist[k], 'mother ID')

        gen_list = []
        gen_list.append([l])
        d1_id = datalist[k][d_idx1][l]
        d2_id = datalist[k][d_idx2][l]
        if np.isnan(d1_id):
            continue
        if np.isnan(d2_id):
            continue
        gen_list.append([int(d1_id) - 1, int(d2_id) - 1])

        if gen_no != 0:
            while (len(gen_list) < gen_no):
                templist = []
                for i in range(len(gen_list[-1])):
                    d1_id = datalist[k][d_idx1][gen_list[-1][i]]
                    d2_id = datalist[k][d_idx2][gen_list[-1][i]]
                    if np.isnan(d1_id):
                        pass
                    else:
                        templist.append(int(d1_id) - 1)
                    if np.isnan(d1_id):
                        pass
                    else:
                        templist.append(int(d2_id) - 1)
                if len(templist) == 0:
                    break
                gen_list.append(templist)

        lin_list = []
        
        for i in range(len(gen_list[-1])):
            lin_list_temp = []
            lin_list_temp.append(gen_list[-1][i])
            m_id = int(datalist[k][m_idx][gen_list[-1][i]] - 1)
            lin_list_temp.append(m_id)
            
            while m_id != -1:
                m_id = int(datalist[k][m_idx][m_id] - 1)
                if m_id != -1:
                    lin_list_temp.append(m_id)
                else:
                    pass
            else:
                pass
            lin_list.append(lin_list_temp)
        gen_list_temp1.append(gen_list)
        lin_list_temp1.append(lin_list)
    gen_list_all.append(gen_list_temp1)
    lin_list_all.append(lin_list_temp1)
        

#%% Time Delay Correlation
mcherry_link = []
venus_link = []

for i in range(len(lin_list_all)):
    for j in range(len(lin_list_all[i])):
        for k in range(len(lin_list_all[i][j])):
            lin_list_all[i][j][k].reverse()

for i in range(len(lin_list_all)):
    mcherry_link_temp1 = []
    venus_link_temp1 = []
    for j in range(len(lin_list_all[i])):
        mcherry_link_temp = []
        venus_link_temp = []
        for k in range(len(lin_list_all[i][j])):
            mcherry_temp = []
            venus_temp = []
            for m in range(len(lin_list_all[i][j][k])):
                mcherry_temp.append(mcherry[i][lin_list_all[i][j][k][m]][~np.isnan(mcherry[i][lin_list_all[i][j][k][m]])].tolist())
                venus_temp.append(venus[i][lin_list_all[i][j][k][m]][~np.isnan(venus[i][lin_list_all[i][j][k][m]])].tolist())
            mcherry_temp = [item for sublist in mcherry_temp for item in sublist]
            venus_temp = [item for sublist in venus_temp for item in sublist]
            mcherry_link_temp.append(mcherry_temp)
            venus_link_temp.append(venus_temp)
        mcherry_link_temp1.append(mcherry_link_temp)
        venus_link_temp1.append(venus_link_temp)
    mcherry_link.append(mcherry_link_temp1)
    venus_link.append(venus_link_temp1)

mcherry_shift_all = []
venus_shift_all = []

for i in range(len(lin_list_all)):
    mcherry_shift_temp1 = []
    venus_shift_temp1 = []
    for j in range(len(lin_list_all[i])):
        mcherry_shift = []
        venus_shift = []
        for k in range(len(lin_list_all[i][j])):
            if numframe > 0:
                mcherry_shift.append(mcherry_link[i][j][k][numframe:])
                venus_shift.append(venus_link[i][j][k][:-numframe])
            elif numframe < 0:
                mcherry_shift.append(mcherry_link[i][j][k][:numframe])
                venus_shift.append(venus_link[i][j][k][-numframe:])
            elif numframe == 0:
                mcherry_shift.append(mcherry_link[i][j][k])
                venus_shift.append(venus_link[i][j][k])
        mcherry_shift_temp1.append(mcherry_shift)
        venus_shift_temp1.append(venus_shift)
    mcherry_shift_all.append(mcherry_shift_temp1)
    venus_shift_all.append(venus_shift_temp1)

#%% Plotting
scale = 2

imagepath = 'Path to folder where images will be stored'

if os.path.exists(imagepath) == False:
    os.mkdir(imagepath)

seaborn.set_style("white")
sb = seaborn.hls_palette(8, l=.3, s=.8)

fig = plt.gcf()
fig.set_size_inches(8*scale,5*scale)
ax = fig.add_subplot(111)

counts,xbins,ybins = np.histogram2d(venus_shift_all[exp_no][cell_no][line_no], 
                                mcherry_shift_all[exp_no][cell_no][line_no], 
                                bins=100, range = ([90, 160],[0, 20000]), normed=True)
counts = ndimage.gaussian_filter(counts, sigma=3.0, order=0)

plt.contourf(np.transpose(counts),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],
             levels = np.linspace(np.max(counts)/50, np.max(counts), 6), cmap=plt.cm.Reds, alpha=1.)
plt.scatter(venus_shift_all[exp_no][cell_no][line_no], 
            mcherry_shift_all[exp_no][cell_no][line_no], marker = 'o', 
            color = sb[5], alpha = 1, s = 50)

plt.xticks(fontsize = 15*scale)
plt.yticks(fontsize = 15*scale)
plt.xlabel("Venus (AU)", fontsize=20*scale)
plt.ylabel("mCherry (AU)", fontsize=20*scale)
fig.set_dpi(1200)
#plt.savefig(os.path.join(imagepath, '051116' + '_' + str(exp_no) + '_' + str(cell_no) + '_' + str(line_no) + '_' + str(timediff) + '.png')
#                             , bbox_inches='tight')
plt.close()

