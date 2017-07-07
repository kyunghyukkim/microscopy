# -*- coding: utf-8 -*-

# Creates line plots for each branches for both mCherry and Venus
# Kiri Choi, 2017

from __future__ import print_function, division

import os
import matplotlib.pyplot as plt
import numpy as np
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
    row_data3dlist = []
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
        row_data3d = []
        for row_number in range(len(column)):            
            row_data3d = ''.join(map(unichr, f[column[row_number]][:]))
        row_data3dlist.append(row_data3d.split(': ')[1])
    def3dlist.append(row_data3dlist)
    f.close()

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
cell_o = 3
gen_no = 8

d_idx1 = idx(deflist[exp_no], 'daughter1 ID')
d_idx2 = idx(deflist[exp_no], 'daughter2 ID')
m_idx = idx(deflist[exp_no], 'mother ID')

seaborn.set_style("white")
sb = seaborn.hls_palette(8, l=.3, s=.8)

gen_list = []
gen_list.append([cell_o])
d1_id = datalist[exp_no][d_idx1][cell_o]
d2_id = datalist[exp_no][d_idx2][cell_o]
if np.isnan(d1_id):
    raise Exception()
if np.isnan(d1_id):
    raise Exception()
gen_list.append([int(d1_id) - 1, int(d2_id) - 1])

if gen_no != 0:
    while (len(gen_list) < gen_no):
        templist = []
        for i in range(len(gen_list[-1])):
            d1_id = datalist[exp_no][d_idx1][gen_list[-1][i]]
            d2_id = datalist[exp_no][d_idx2][gen_list[-1][i]]
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
    m_id = int(datalist[exp_no][m_idx][gen_list[-1][i]] - 1)
    lin_list_temp.append(m_id)
    
    while m_id != -1:
        m_id = int(datalist[exp_no][m_idx][m_id] - 1)
        if m_id != -1:
            lin_list_temp.append(m_id)
        else:
            pass
    else:
        pass
    lin_list.append(lin_list_temp)
    
gen_list_flat = [item for sublist in gen_list for item in sublist]

#%% Plotting

scale=1

plt.ioff()

imagepath = 'Path to folder where images will be stored'

if os.path.exists(imagepath) == False:
    os.mkdir(imagepath)

for i in range(len(lin_list)):
    testlist = lin_list[i]

    fig, ax1 = plt.subplots(figsize = (16, 10))
    ax2 = ax1.twinx()
    
    for j in range(len(testlist)):
        ax1.plot(mcherry[exp_no][testlist[j]], color = sb[4], label="mCherry")
        ax2.plot(venus[exp_no][testlist[j]], color = sb[3], label="Venus")
    
    plt.title("mCherry and Venus", fontsize=20*scale)
    #fig.savefig(os.path.join(imagepath, str(exp_no) + '_' + str(cell_o) + '_' + str(i) + '.png')
    #                         , bbox_inches='tight')
    plt.close(fig)
    
