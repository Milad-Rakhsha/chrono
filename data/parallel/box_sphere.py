import csv,os,sys
import subprocess,re


import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
from decimal import Decimal
from collections import OrderedDict
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 18})
plt.rc('xtick',labelsize=24)
plt.rc('ytick',labelsize=24)
MARKERSIZE=5
path = str(sys.argv[1])
title = str(sys.argv[2])

frame=10

def prepare(path):
    cmd="ls " + path
    print(cmd)
    process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    text = output.decode('utf-8')
    out = re.sub(r',', r'\n', text)
    cases = []
    for l in out.split("\n"):
        if (l not in cases and l != ""):
            cases.append(l)

    files = OrderedDict([])
    for i in cases:
        param = pd.read_csv(path+i+"/params.txt")
        files[i]=dict([])
        files[i]["regularize"]=param["reg"][0]
        files[i]["solver"]=param["solver"][0]
        files[i]["alpha_0"]=param["reg_alpha0"][0]
        files[i]["NC"]=param["NC"][0]**2
        files[i]["randomized"]=param["rand"][0]
        files[i]["regularize_t"]=param["reg_tan"][0]
        files[i]["compliance"]=param["comp"][0]
        # print (files[i])
    return cases,files

cases,files=prepare(path)
N=12
mean_n=np.zeros((N,2))
mean_t=np.zeros((N,2))
std_n=np.zeros((N,2))
std_t=np.zeros((N,2))
x_ax=np.zeros((N,2))

idx_reg=0
idx=0
W=10
for i in files:
    reg=files[i]["regularize"]
    nc=int(np.sqrt(files[i]["NC"]))
    data = pd.read_csv(path+i+"/_forces"+str(frame)+".txt")
    FN_ex=W / nc**2
    Ft_ex=W /2 * 0.5 / nc**2
    Fn_num=np.mean(data["Fn"])
    Ft_num=np.mean(np.sqrt(np.power(data["Ft1"],2)+np.power(data["Ft2"],2)))

    Fn_std=np.std(data["Fn"]-FN_ex)
    Ft_std=np.std(np.sqrt(np.power(data["Ft1"],2)+np.power(data["Ft2"],2)))

    x_ax[nc,reg]=files[i]["NC"]
    mean_n[nc,reg]=np.abs(Fn_num-FN_ex)/FN_ex*100
    mean_t[nc,reg]=np.abs(Ft_num-Ft_ex)/Ft_ex*100
    std_n[nc,reg]=Fn_std
    std_t[nc,reg]=Ft_std



x_ax=x_ax[2:,:]
std_n=std_n[2:,:]
std_t=std_t[2:,:]
mean_n=mean_n[2:,:]
mean_t=mean_t[2:,:]


BLUE='tab:blue'
RED='tab:red'


fig = plt.figure(num=None,figsize=(8, 6),  dpi=100, facecolor='w', edgecolor='k')
ax1 = fig.add_subplot(111)
ax1.plot(x_ax[:,0],mean_n[:,0],
        "ro-",
        linewidth=1, markersize=MARKERSIZE,label='Normal'
        )
ax1.plot(x_ax[:,0],mean_t[:,0],
        "r+-",
        linewidth=1, markersize=MARKERSIZE,label='Tangential'
        )
ax2 = ax1.twinx()
ax2.plot(x_ax[:,1],mean_n[:,1],
        "bo-",
        linewidth=1, markersize=MARKERSIZE,label='Normal (Tikhonov)'
        )
ax2.plot(x_ax[:,1],mean_t[:,1],
        "b+-",
        linewidth=1, markersize=MARKERSIZE,label='Tangential (Tikhonov)'
        )
ax1.set_xlabel('$n$')
ax1.set_ylabel(r'err%$=\frac{\|x-\bar{x}\|}{\bar{x}}\times 100$', color=RED,fontsize=22)
ax2.set_ylabel(r'err%$=\frac{\|x-\bar{x}\|}{\bar{x}}\times 100$, Tikhonov',color=BLUE,fontsize=22)
ax1.tick_params(axis='y', labelcolor=RED)
ax2.tick_params(axis='y', labelcolor=BLUE)
ax2.legend(loc=1)
ax1.legend(loc=0)
# plt.legend( fancybox=True, shadow=True, ncol=1)
plt.title(title )
plt.savefig(title+'_mean.png', bbox_inches='tight')
plt.show()

plt.cla()   # Clear axis
plt.clf()   # Clear figure
plt.close() # Close a figure window

fig = plt.figure(num=None,figsize=(8, 6),  dpi=100, facecolor='w', edgecolor='k')
ax1 = fig.add_subplot(111)
ax1.plot(x_ax[:,0],std_n[:,0],
        "ro-",
        linewidth=1, markersize=MARKERSIZE,label='Normal'
        )
ax1.plot(x_ax[:,0],std_t[:,0],
        "r+-",
        linewidth=1, markersize=MARKERSIZE,label='Tangential'
        )
ax2 = ax1.twinx()
ax2.plot(x_ax[:,1],std_n[:,1],
        "bo-",
        linewidth=1, markersize=MARKERSIZE,label='Normal (Tikhonov)'
        )
ax2.plot(x_ax[:,1],std_t[:,1],
        "b+-",
        linewidth=1, markersize=MARKERSIZE,label='Tangential (Tikhonov)'
        )

ax1.set_xlabel('$n$')
ax1.set_ylabel(r'$\sigma^2=\frac{\Sigma^n_i(x_i-\bar{x})^2}{n}$', color=RED,fontsize=22)
ax2.set_ylabel(r'$\sigma^2=\frac{\Sigma^n_i(x_i-\bar{x})^2}{n}$, Tikhonov',color=BLUE,fontsize=22)
ax1.tick_params(axis='y', labelcolor=RED)
ax2.tick_params(axis='y', labelcolor=BLUE)
ax2.legend(loc=1)
ax1.legend(loc=2)
# plt.legend( fancybox=True, shadow=True, ncol=1)
plt.title(title )

plt.savefig(title+'_std.png', bbox_inches='tight')
plt.show()
