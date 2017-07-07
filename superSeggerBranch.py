# -*- coding: utf-8 -*-
"""
Creates line plots for each branches for both mCherry and Venus
"""

import sys, os
import random
import scipy
import scipy.stats
import numpy as np
import math
import scipy.optimize as opt
from scipy import ndimage
import seaborn
import h5py
import tables
import re

#clistpathlist = [r'J:\Kyung\060916\crop2\xy1\uclist.mat']

clistpathlist = [r'J:\Kyung\051116\crop1\xy1\uclist.mat']#, 
                 #r'J:\Kyung\051116\crop2\xy1\uclist.mat',
                 #r'J:\Kyung\051116\crop3\xy1\uclist.mat',
                 #r'J:\Kyung\051116\crop4\xy1\uclist.mat',
                 #r'J:\Kyung\051116\crop5\xy1\uclist.mat',
                 #r'J:\Kyung\051116\crop6\xy1\uclist.mat',
                 #r'J:\Kyung\051116\crop7\xy1\uclist.mat']
                 
#clistpathlist = [r'J:\Kyung\051616\crop1\xy1\uclist.mat',
#                 r'J:\Kyung\051616\crop2\xy1\uclist.mat',
#                 r'J:\Kyung\051616\crop3\xy1\uclist.mat',
#                 r'J:\Kyung\051616\crop4\xy1\uclist.mat']

#clistpathlist = [r'J:\Kyung\051716\crop1\xy1\uclist.mat',
#                 r'J:\Kyung\051716\crop2\xy1\uclist.mat',
#                 r'J:\Kyung\051716\crop3\xy1\uclist.mat',
#                 r'J:\Kyung\051716\crop4\xy1\uclist.mat',
#                 r'J:\Kyung\051716\crop5\xy1\uclist.mat']

#clistpathlist = [r'J:\Kyung\052516\p1\crop1\xy1\uclist.mat',
#                 r'J:\Kyung\052516\p1\crop2\xy1\uclist.mat',
#                 r'J:\Kyung\052516\p2\crop1\xy2\uclist.mat']

#clistpathlist = [r'J:\Kyung\060916\crop1\xy1\uclist.mat']

#r'J:\Kyung\050316\Venusmcherry\xy2\crop1\xy1\uclist.mat'

#data = scipy.io.loadmat(clistpath)

def idx(string, char):
  for key, x in enumerate(string):
    if char in x:
      return key

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

#%% exclude
#if len(exclist) > 0:
#    for i in range(len(exclist)):
#        datalist[i] = np.delete(datalist[i], exclist[i], 1)
#        for j in range(len(exclist[i])):
#            data3dlist[i][j][:,np.int(exclist[i][j])] = np.nan

#%% Populate list
mcherry = []
venus = []

for i in range(len(clistpathlist)):
    mcherry.append(data3dlist[i][:,idx(def3dlist[i], 'fluor1 mean')].T)
    venus.append(data3dlist[i][:,idx(def3dlist[i], 'fluor2 mean')].T),


#%% Check original cells
cell_o_list = []

for i in range(len(clistpathlist)):
    cell_o_list.append((np.isnan(data3dlist[i][:,idx(def3dlist[i], 'fluor1 mean')][0,:]) == False).nonzero())
    
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

#%% Figures
# Export all into a folder

scale=1

import matplotlib as mlt
import matplotlib.pyplot as plt

plt.ioff()

imagepath = r'J:\Images\051116_branch_pdf'

if os.path.exists(imagepath) == False:
    os.mkdir(imagepath)

for i in range(len(lin_list)):
    testlist = lin_list[i]

    fig, ax1 = plt.subplots(figsize = (16, 10))
    ax2 = ax1.twinx()
    
    for j in range(len(testlist)):
        ax1.plot(mcherry[exp_no][testlist[j]], color = sb[4], label="mCherry")
        ax2.plot(venus[exp_no][testlist[j]], color = sb[3], label="Venus")
        #plt.plot(mcherry[exp_no][gen_list_flat[i]])
    
    plt.title("mCherry and Venus", fontsize=20*scale)
    fig.savefig(os.path.join(imagepath, str(exp_no) + '_' + str(cell_o) + '_' + str(i) + '.png')
                             , bbox_inches='tight')
    plt.close(fig)
    

#%%
##testlist = [3, 7, 26, 39]
##testlist = [2, 5, 11, 15, 25]
##testlist = [1, 3, 7, 10, 19]
#testlist = lin_list[0]
##testlist = [1, 2, 9, 12, 26, 36, 86]
##testlist = [0, 2, 3, 6, 11, 37, 42, 75, 125, 217, 438]
#
##fig = plt.gcf()
##fig.set_size_inches(8*scale,5*scale)
##ax = fig.add_subplot(111)
#fig, ax1 = plt.subplots()
#ax2 = ax1.twinx()
#
#for i in range(len(testlist)):
#    ax1.plot(mcherry[exp_no][testlist[i]], color = sb[4], label="mCherry")#, alpha = 0.5)
#    ax2.plot(venus[exp_no][testlist[i]], color = sb[3], label="Venus")#, alpha = 0.5)
#    #plt.plot(mcherry[exp_no][gen_list_flat[i]])
#
#plt.title("mCherry and Venus", fontsize=20*scale)
#plt.show()
    
    
    
#for i in range(len(testlist)):
#    plt.plot(mcherry[exp_no][testlist[i]])
#
#plt.title("mCherry sum", fontsize=20*scale)
#plt.xticks(fontsize = 15*scale)
#plt.yticks(fontsize = 15*scale)
##plt.title("KK61", fontsize=80)
#plt.xlabel("Time (5 mins)", fontsize=20*scale)
#plt.ylabel("Fluoresence Intensity", fontsize=20*scale)
##plt.axis([0, 85, 100, 115])
##fig.set_dpi(1200)
#
#plt.show()
##ax.xaxis.set_tick_params(length=5*scale)
##ax.yaxis.set_tick_params(length=5*scale)
##ax.tick_params(axis='x', pad=10)
##ax.tick_params(axis='y', pad=20)
##ax.tick_params('both', length=15, width=3, which='major')
##ax.tick_params('both', length=10, width=3, which='minor')
##ax.spines['top'].set_linewidth(3)
##ax.spines['right'].set_linewidth(3)
##ax.spines['left'].set_linewidth(3)
##ax.spines['bottom'].set_linewidth(3)
#
#for i in range(len(testlist)):
#    plt.plot(venus[exp_no][testlist[i]])
#
#plt.title("Venus mean", fontsize=20*scale)
#plt.xticks(fontsize = 15*scale)
#plt.yticks(fontsize = 15*scale)
##plt.title("KK61", fontsize=80)
#plt.xlabel("Time (5 mins)", fontsize=20*scale)
#plt.ylabel("Fluoresence Intensity", fontsize=20*scale)
##plt.axis([0, 85, 100, 115])
##fig.set_dpi(1200)
#
#plt.show()

print(gen_list)
