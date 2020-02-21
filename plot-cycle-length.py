
import methods
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import ezodf
import pandas as pd
import re
import glob
import sys

FileName = "./optical/BeatingRate.ods"
try:
    xlsx = pd.ExcelFile(FileName)
    allDrugs = xlsx.sheet_names
except:
    ods = ezodf.opendoc(FileName)
    allDrugs = [a.name for a in ods.sheets]

CONC_ID = [(1,5),(6,10),(11,15),(16,20)]
plotchange = not True

for (i,drug) in enumerate(allDrugs):
    try:
        sheet = xlsx.parse(drug)
        idx = np.array(sheet)[:,0]
        control = np.array(sheet)[:,1]
        applied = np.array(sheet)[:,2]
    except:
        sheet = ods.sheets[drug]
        idx = np.array([int(i[0].value) for i in sheet.rows()])
        control = np.array([i[1].value for i in sheet.rows()])
        applied = np.array([i[2].value for i in sheet.rows()])
    plt.figure()
    for (j,conc_id) in enumerate(CONC_ID):
        conc_idx = np.arange(np.max(idx))[((idx>=conc_id[0]) & (idx<=conc_id[1]))]
        plt.xlabel(r"concentration")
        plt.title(drug)
        if plotchange:
            plt.plot([j+1]*len(conc_idx), applied[conc_idx] - control[conc_idx], 'ro')
            plt.ylabel(r"$\Delta$ cycle length [ms]")
            plt.xlim([0.1,4.9])
            methods.save('cyclelength_Diff_%s'%drug)
        else:
            plt.plot([j+1]*len(conc_idx), applied[conc_idx], 'kx')
            plt.plot([0]*len(conc_idx), control[conc_idx], 'kx')
            plt.ylabel(r"cycle length [ms]")
            plt.xlim([-0.9,4.9])
            methods.save('cyclelength_%s'%drug)
#plt.show()
