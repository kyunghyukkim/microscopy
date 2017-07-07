"""
Scatter plots with interpolation
"""

import sys, os
import random
import matplotlib as mlt
import matplotlib.pyplot as plt
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
import itertools

#clistpathlist = [r'J:\Kyung\051116\crop6\xy1\uclist.mat']

clistpathlist = [r'J:\Kyung\051116\crop1\xy1\uclist.mat',
                 r'J:\Kyung\051116\crop2\xy1\uclist.mat',
                 r'J:\Kyung\051116\crop3\xy1\uclist.mat',
                 r'J:\Kyung\051116\crop4\xy1\uclist.mat',
                 r'J:\Kyung\051116\crop5\xy1\uclist.mat',
                 r'J:\Kyung\051116\crop6\xy1\uclist.mat',
                 r'J:\Kyung\051116\crop7\xy1\uclist.mat']



#                 r'J:\Kyung\051716\crop1\xy1\uclist.mat',
#                 r'J:\Kyung\051716\crop2\xy1\uclist.mat',
#                 r'J:\Kyung\051716\crop3\xy1\uclist.mat',
#                 r'J:\Kyung\051716\crop4\xy1\uclist.mat',
#                 r'J:\Kyung\051716\crop5\xy1\uclist.mat']

#clistpathlist = [r'J:\Kyung\051716\crop1\xy1\uclist.mat',
#                 r'J:\Kyung\051716\crop2\xy1\uclist.mat',
#                 r'J:\Kyung\051716\crop3\xy1\uclist.mat',
#                 r'J:\Kyung\051716\crop4\xy1\uclist.mat',
#                 r'J:\Kyung\051716\crop5\xy1\uclist.mat']

#                 r'J:\Kyung\052516\p1\crop1\xy1\uclist.mat',
#                 r'J:\Kyung\052516\p1\crop2\xy1\uclist.mat',
                 #r'J:\Kyung\052516\p2\crop1\xy2\uclist.mat']
#clistpathlist = [r'J:\Kyung\051616\crop1\xy1\uclist.mat',
#                 r'J:\Kyung\051616\crop2\xy1\uclist.mat',
#                 r'J:\Kyung\051616\crop3\xy1\uclist.mat',
#                 r'J:\Kyung\051616\crop4\xy1\uclist.mat']

#r'J:\Kyung\050316\Venusmcherry\xy2\crop1\xy1\uclist.mat'

def idx(string, char):
  for key, x in enumerate(string):
    if char in x:
      return key

def quin(x, a, b, c, d, e):
    return a*np.power(x,4) + b*np.power(x,3) + c*np.square(x) + d*x + e

def cubic(x, a, b, c, d):
    return a*np.power(x,3) + b*np.square(x) + c*x + d
    
def quad(x, a, b, c):
    return a*np.square(x) + b*x + c

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
    
#%% exclude
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
    cell_o_list.append((np.isnan(data3dlist[i][:,idx(def3dlist[i], 'fluor1 mean')][0,:]) == False).nonzero())
    
print(cell_o_list)

#%% Track lineage
exp_no = 0
cell_no = 0
line_no = 0
gen_no = 8
timediff1 = 0
timediff2 = 60
timediff3 = 75
timescale = 5
numframe1 = timediff1/timescale
numframe2 = timediff2/timescale
numframe3 = timediff3/timescale


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
        

#%%
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

#%% Interpolation
#params = []
#cov = []

diff_N = 100.0

mcherry_link_int = []
mcherry_link_int_diff = []

for i in range(len(lin_list_all)):
    mcherry_link_int_temp = []
    mcherry_link_int_diff_temp = []
    #params_temp = []
    #cov_temp = []
    for j in range(len(lin_list_all[i])):
        mcherry_link_int_temp1 = []
        mcherry_link_int_diff_temp1 = []
        #params_temp1 = []
        #cov_temp1 = []
        for k in range(len(lin_list_all[i][j])):
            #popt, pcov = scipy.optimize.curve_fit(quin, np.arange(len(mcherry_link[i][j][k])), mcherry_link[i][j][k])
            x = np.linspace(0, len(mcherry_link[i][j][k]) - 1, len(mcherry_link[i][j][k]))
            #y = quin(np.arange(len(mcherry_link[i][j][k])), popt[0], popt[1], popt[2], popt[3], popt[4])
            f_int = scipy.interpolate.interp1d(x, mcherry_link[i][j][k], kind='cubic')
            int_mcherry = f_int(np.linspace(0, len(mcherry_link[i][j][k]) - 1, 100*len(mcherry_link[i][j][k])))
            #plt.plot(np.linspace(0, len(mcherry_link[i][j][k]) - 1, 100*len(mcherry_link[i][j][k])), int_mcherry)
            #plt.plot(np.linspace(0, len(mcherry_link[i][j][k]) - 1, len(mcherry_link[i][j][k])), mcherry_link[i][j][k])
            #plt.show()
            int_mcherry_diff = np.gradient(int_mcherry)
            #plt.plot()
            #params_temp1.append(popt.tolist())
            #cov_temp1.append(pcov.tolist())

            #plt.plot(mcherry_link[i][j][k])
            #plt.plot(quin(np.arange(len(mcherry_link[i][j][k])), popt[0], popt[1], popt[2], popt[3], popt[4]))
            #plt.show()
            mcherry_link_int_temp1.append(int_mcherry)
            mcherry_link_int_diff_temp1.append(int_mcherry_diff)
        mcherry_link_int_temp.append(mcherry_link_int_temp1)
        mcherry_link_int_diff_temp.append(mcherry_link_int_diff_temp1)
        #params_temp.append(params_temp1)
        #cov_temp.append(cov_temp1)
    mcherry_link_int.append(mcherry_link_int_temp)
    mcherry_link_int_diff.append(mcherry_link_int_diff_temp)
    #params.append(params_temp)
    #cov.append(cov_temp)
    
    
#%%
mcherry_link_int_diff_reduced = []

for i in range(len(mcherry_link_int_diff)):
    mcherry_link_int_diff_temp = []
    for j in range(len(mcherry_link_int_diff[i])):
        mcherry_link_int_diff_temp1 = []
        for k in range(len(mcherry_link_int_diff[i][j])):
            mcherry_link_int_diff_temp2 = []
            for l in range(len(mcherry_link_int_diff[i][j][k])):
                if l/diff_N == int(l/diff_N):
                    mcherry_link_int_diff_temp2.append(mcherry_link_int_diff[i][j][k][l])
                else:
                    pass
            mcherry_link_int_diff_temp1.append(mcherry_link_int_diff_temp2)
        mcherry_link_int_diff_temp.append(mcherry_link_int_diff_temp1)
    mcherry_link_int_diff_reduced.append(mcherry_link_int_diff_temp)
    
mcherry_shift_all1 = []
venus_shift_all1 = []
b = 0.02

for i in range(len(lin_list_all)):
    mcherry_shift_temp1 = []
    venus_shift_temp1 = []
    for j in range(len(lin_list_all[i])):
        mcherry_shift = []
        venus_shift = []
        for k in range(len(lin_list_all[i][j])):
            if numframe1 > 0:
                mcherry_shift.append(np.add(np.multiply(b, mcherry_link[i][j][k][numframe1:]), mcherry_link_int_diff_reduced[i][j][k][numframe1:]))
                        #mcherry_link_int_diff_reduced[i][j][k][numframe1:])
                venus_shift.append(venus_link[i][j][k][:-numframe1])
            elif numframe1 < 0:
                mcherry_shift.append(np.add(np.multiply(b, mcherry_link[i][j][k][:numframe1]), mcherry_link_int_diff_reduced[i][j][k][:numframe1]))
                        #mcherry_link_int_diff_reduced[i][j][k][:numframe1])
                venus_shift.append(venus_link[i][j][k][-numframe1:])
            elif numframe1 == 0:
                mcherry_shift.append(np.add(np.multiply(b, mcherry_link[i][j][k]), mcherry_link_int_diff_reduced[i][j][k]))
                        #mcherry_link_int_diff_reduced[i][j][k])
                venus_shift.append(venus_link[i][j][k])
        mcherry_shift_temp1.append(mcherry_shift)
        venus_shift_temp1.append(venus_shift)
    mcherry_shift_all1.append(mcherry_shift_temp1)
    venus_shift_all1.append(venus_shift_temp1)
    
mcherry_shift_all2 = []
venus_shift_all2 = []

for i in range(len(lin_list_all)):
    mcherry_shift_temp1 = []
    venus_shift_temp1 = []
    for j in range(len(lin_list_all[i])):
        mcherry_shift = []
        venus_shift = []
        for k in range(len(lin_list_all[i][j])):
            if numframe2 > 0:
                mcherry_shift.append(np.add(np.multiply(b, mcherry_link[i][j][k][numframe2:]), mcherry_link_int_diff_reduced[i][j][k][numframe2:]))
                        #np.multiply(b, mcherry_link[i][j][k][numframe2:]))
                venus_shift.append(venus_link[i][j][k][:-numframe2])
            elif numframe2 < 0:
                mcherry_shift.append(np.add(np.multiply(b, mcherry_link[i][j][k][:numframe2]), mcherry_link_int_diff_reduced[i][j][k][:numframe2]))
                        #np.multiply(b, mcherry_link[i][j][k][:numframe2]))
                venus_shift.append(venus_link[i][j][k][-numframe2:])
            elif numframe2 == 0:
                mcherry_shift.append(np.add(np.multiply(b, mcherry_link[i][j][k]), mcherry_link_int_diff_reduced[i][j][k]))
                        #np.multiply(b, mcherry_link[i][j][k]))
                venus_shift.append(venus_link[i][j][k])
        mcherry_shift_temp1.append(mcherry_shift)
        venus_shift_temp1.append(venus_shift)
    mcherry_shift_all2.append(mcherry_shift_temp1)
    venus_shift_all2.append(venus_shift_temp1)
    
mcherry_shift_all3 = []
venus_shift_all3 = []

for i in range(len(lin_list_all)):
    mcherry_shift_temp1 = []
    venus_shift_temp1 = []
    for j in range(len(lin_list_all[i])):
        mcherry_shift = []
        venus_shift = []
        for k in range(len(lin_list_all[i][j])):
            if numframe3 > 0:
                mcherry_shift.append(np.add(np.multiply(b, mcherry_link[i][j][k][numframe3:]), mcherry_link_int_diff_reduced[i][j][k][numframe3:]))
                venus_shift.append(venus_link[i][j][k][:-numframe3])
            elif numframe3 < 0:
                mcherry_shift.append(np.add(np.multiply(b, mcherry_link[i][j][k][:numframe3]), mcherry_link_int_diff_reduced[i][j][k][:numframe3]))
                venus_shift.append(venus_link[i][j][k][-numframe3:])
            elif numframe3 == 0:
                mcherry_shift.append(np.add(np.multiply(b, mcherry_link[i][j][k]), mcherry_link_int_diff_reduced[i][j][k]))
                venus_shift.append(venus_link[i][j][k])
        mcherry_shift_temp1.append(mcherry_shift)
        venus_shift_temp1.append(venus_shift)
    mcherry_shift_all3.append(mcherry_shift_temp1)
    venus_shift_all3.append(venus_shift_temp1)

#mcherry_shift_all_flat = [item for sublist in mcherry_shift_all for item in sublist]
#venus_shift_all_flat = [item for sublist in venus_shift_all for item in sublist]
#
#for i in range(2):
#    mcherry_shift_all_flat = [item for sublist in mcherry_shift_all_flat for item in sublist]
#    venus_shift_all_flat = [item for sublist in venus_shift_all_flat for item in sublist]
#
#test, unique_ind = np.unique(mcherry_shift_all_flat, return_index = True)
#re_ind = np.arange(len(venus_shift_all_flat))
#re_ind = np.delete(re_ind, np.sort(unique_ind))
#mcherry_shift_all_flat = np.delete(mcherry_shift_all_flat, re_ind)
#venus_shift_all_flat = np.delete(venus_shift_all_flat, re_ind)
#%%
scale = 2

imagepath = r'J:\Images\shifted_dots_051116_interp_iden_pdf'

if os.path.exists(imagepath) == False:
    os.mkdir(imagepath)

seaborn.set_style("white")
sb = seaborn.hls_palette(8, l=.3, s=.8)

#plt.ioff()

for i in range(len(venus_shift_all1)):
    for j in range(len(venus_shift_all1[i])):
        for k in range(len(venus_shift_all1[i][j])):
            fig = plt.gcf()
            fig.set_size_inches(8*scale,5*scale)
            ax = fig.add_subplot(111)
            #ax.xaxis.set_tick_params(length=5*scale)
            #ax.yaxis.set_tick_params(length=5*scale)
            #ax.tick_params(axis='x', pad=10)
            #ax.tick_params(axis='y', pad=20)
            #ax.tick_params('both', length=15, width=3, which='major')
            #ax.tick_params('both', length=10, width=3, which='minor')
            #ax.spines['top'].set_linewidth(3)
            #ax.spines['right'].set_linewidth(3)
            #ax.spines['left'].set_linewidth(3)
            #ax.spines['bottom'].set_linewidth(3)
            
            #counts,xbins,ybins = np.histogram2d(venus_shift_all[i][j][k], mcherry_shift_all[i][j][k], bins=100, range = ([90, 160],[0, 20000]), normed=True)
            #counts = ndimage.gaussian_filter(counts, sigma=3.0, order=0)
            
            #plt.contourf(np.transpose(counts),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],
            #             levels = np.linspace(np.max(counts)/50, np.max(counts), 6), cmap=plt.cm.Reds, alpha=1.)
            
            #plt.scatter(venus_shift_all1[i][j][k], mcherry_shift_all1[i][j][k], marker = 'o', color = sb[5], alpha = 1, s = 50)
            #plt.scatter(venus_shift_all2[i][j][k], mcherry_shift_all2[i][j][k], marker = 'o', color = sb[3], alpha = 1, s = 50)
            plt.scatter(venus_shift_all3[i][j][k], mcherry_shift_all3[i][j][k], marker = 'o', color = sb[0], alpha = 1, s = 50)
            
            #plt.plot(venus_shift_all[exp_no][cell_no][line_no], mcherry_shift_all[exp_no][cell_no][line_no], ls = '', marker = 'o', color = sb[5], alpha = 1)
            
            plt.xticks(fontsize = 15*scale)
            plt.yticks(fontsize = 15*scale)
            #plt.title("KK61", fontsize=80)
            plt.xlabel("Venus (AU)", fontsize=20*scale)
            plt.ylabel("mCherry (AU)", fontsize=20*scale)
            #fig.set_dpi(1200)
            plt.axis([100, 145, 0, 300])
            #fig.savefig(os.path.join(imagepath, '051116' + '_' + str(i) + '_' + str(j) + '_' + str(k) + '_75.pdf')
            #                             , bbox_inches='tight')
            plt.close(fig)
            #plt.show()

