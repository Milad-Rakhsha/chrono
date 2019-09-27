import csv,os,sys
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

path = str(sys.argv[1])

DEM=path+"CannonballSMC/F_SCM_49.txt"
DVI=path+"CannonballNSC/F_NSC_49.txt"


files = OrderedDict([#
    (DVI,   {"axis": "x", "shift": +0.65, "dT": 0.05, "lineStyle": 'ko-', "markerEvery": 1, 'markersize':10, 'label': 'DVI'}),
    (DEM,   {"axis": "x", "shift": +0.65, "dT": 0.05, "lineStyle": 'b-',  "markerEvery": 1, 'markersize':5, 'label': 'DEM'}),
])
fig = plt.figure(num=None,figsize=(18, 12),  dpi=300, facecolor='w', edgecolor='k')
ax1 = fig.add_subplot(211)
ax1.autoscale(enable=True, axis='x', tight=True)

dem = pd.read_csv(DEM)
dvi = pd.read_csv(DVI)
dem["Fn"]

ax1.plot(dem["Fn"], dvi["Fn"],
          "ko-",
          linewidth=3, markersize=5,
          )
# ax1.set_xticks(np.linspace(0, 0.4, 9))
# ax1.set_yticks(np.linspace(0, 50, 11))
ax1.grid(which='both', linestyle='--', linewidth=0.5)
# ax1.set_xlim(0, 0.35)
# ax1.set_ylim(0, 20)
# ax1.grid(which='minor', alpha=0.2)
ax1.grid(which='major', alpha=0.3)
ax1.grid(color='k', linestyle='-', linewidth=0.2)
ax1.set_ylabel('$DVI$(mm)', fontsize=36)
ax1.set_ylabel('$DEM$(mm)', fontsize=36)

plt.savefig('CannonBall.png', bbox_inches='tight', fontsize=36)
plt.show()
