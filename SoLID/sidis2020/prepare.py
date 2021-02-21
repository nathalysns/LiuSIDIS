#!/usr/bin/env python

import sys
import numpy as np
import scipy as sp
import pandas as pd

import tmdlib.tmd as tmd

pseudodata = {}
pseudodata['sbs01'] = pd.read_csv('datacollins/sbs01.dat', delim_whitespace=True)
pseudodata['sbs02'] = pd.read_csv('datacollins/sbs02.dat', delim_whitespace=True)
pseudodata['clas01'] = pd.read_csv('datacollins/clas01.dat', delim_whitespace=True)
pseudodata['clas02'] = pd.read_csv('datacollins/clas02.dat', delim_whitespace=True)
pseudodata['basePpip'] = pd.read_csv('datacollins/basePpip.csv', delim_whitespace=False)
pseudodata['basePpim'] = pd.read_csv('datacollins/basePpim.csv', delim_whitespace=False)
pseudodata['baseNpip'] = pd.read_csv('datacollins/baseNpip.csv', delim_whitespace=False)
pseudodata['baseNpim'] = pd.read_csv('datacollins/baseNpim.csv', delim_whitespace=False)
pseudodata['enhancedPpip'] = pd.read_csv('datacollins/enhancedPpip.csv', delim_whitespace=False)
pseudodata['enhancedPpim'] = pd.read_csv('datacollins/enhancedPpim.csv', delim_whitespace=False)
pseudodata['enhancedNpip'] = pd.read_csv('datacollins/enhancedNpip.csv', delim_whitespace=False)
pseudodata['enhancedNpim'] = pd.read_csv('datacollins/enhancedNpim.csv', delim_whitespace=False)

sbs = pd.concat([pseudodata['sbs01'],pseudodata['sbs02']], axis=0, ignore_index=True)
clas = pd.concat([pseudodata['clas01'],pseudodata['clas02']], axis=0, ignore_index=True)
base = pd.concat([pseudodata['basePpip'],pseudodata['basePpim'],pseudodata['baseNpip'],pseudodata['baseNpim']], axis=0, ignore_index=True)
enhanced = pd.concat([pseudodata['enhancedPpip'],pseudodata['enhancedPpim'],pseudodata['enhancedNpip'],pseudodata['enhancedNpim']], axis=0, ignore_index=True)
base3he = pd.concat([pseudodata['baseNpip'],pseudodata['baseNpim']], axis=0, ignore_index=True)
enhanced3he = pd.concat([pseudodata['enhancedNpip'],pseudodata['enhancedNpim']], axis=0, ignore_index=True)

def simulatecollins(data, var):
    val = []
    for i in range(0,len(data)):
        x, y, z, Q2, pT, target, hadron = data.loc[i]['x'], data.loc[i]['y'], data.loc[i]['z'], data.loc[i]['Q2'], data.loc[i]['pT'], data.loc[i]['target'], data.loc[i]['hadron']
        val.append(tmd.AUTCollins(x, y, Q2, z, pT, target, hadron, var))
    data['value'] = val
    return

def simulatesivers(data, var):
    val = []
    for i in range(0,len(data)):
        x, y, z, Q2, pT, target, hadron = data.loc[i]['x'], data.loc[i]['y'], data.loc[i]['z'], data.loc[i]['Q2'], data.loc[i]['pT'], data.loc[i]['target'], data.loc[i]['hadron']
        val.append(tmd.AUTSivers(x, Q2, z, pT, target, hadron, var))
    data['value'] = val
    return

def preparecollins():
    global sbs, clas, base, enhanced
    par0 = {'Nu':0.4, 'Nd':-0.45, 'a':1.0, 'b':3.0, 'c':0., 'kt2':0.25}
    simulatecollins(sbs,par0)
    simulatecollins(clas,par0)
    simulatecollins(base,par0)
    simulatecollins(enhanced,par0)
    simulatecollins(base3he,par0)
    simulatecollins(enhanced3he,par0)
    sbsrep = sbs.copy()
    clasrep = clas.copy()
    baserep = base.copy()
    baserep['error'] = base['stat']
    basesystrep = base.copy()
    basesystrep['error'] = (base['stat']**2+base['systabs']**2+base['value']**2*base['systrel']**2)**0.5
    enhancedrep = enhanced.copy()
    enhancedrep['error'] = enhanced['stat']
    enhancedsystrep = enhanced.copy()
    enhancedsystrep['error'] = (enhanced['stat']**2+enhanced['systabs']**2+enhanced['value']**2*enhanced['systrel']**2)**0.5
    base3herep = base3he.copy()
    base3herep['error'] = base3he['stat']
    base3hesystrep = base3he.copy()
    base3hesystrep['error'] = (base3he['stat']**2+base3he['systabs']**2+base3he['value']**2*base3he['systrel']**2)**0.5
    enhanced3herep = enhanced3he.copy()
    enhanced3herep['error'] = enhanced3he['stat']
    enhanced3hesystrep = enhanced3he.copy()
    enhanced3hesystrep['error'] = (enhanced3he['stat']**2+enhanced3he['systabs']**2+enhanced3he['value']**2*enhanced3he['systrel']**2)**0.5
    sbsrep.to_csv('datacollins/simsbs.dat', sep='\t', index=False)
    clasrep.to_csv('datacollins/simclas.dat', sep='\t', index=False)
    baserep.to_csv('datacollins/simbase.dat', sep='\t', index=False)
    basesystrep.to_csv('datacollins/simbasesyst.dat', sep='\t', index=False)
    enhancedrep.to_csv('datacollins/simenhanced.dat', sep='\t', index=False)
    enhancedsystrep.to_csv('datacollins/simenhancedsyst.dat', sep='\t', index=False)
    base3herep.to_csv('datacollins/simbase3he.dat', sep='\t', index=False)
    base3hesystrep.to_csv('datacollins/simbase3hesyst.dat', sep='\t', index=False)
    enhanced3herep.to_csv('datacollins/simenhanced3he.dat', sep='\t', index=False)
    enhanced3hesystrep.to_csv('datacollins/simenhanced3hesyst.dat', sep='\t', index=False)
    return

def preparesivers():
    global sbs, clas, base, enhanced
    par1 = {'Nu': -0.03851696777000435,'au': 0.6624828400702693,'bu': 4.103081356761233,'cu': 0.0,'Nd': 0.05141101415846694,'ad': 0.3984321190429335,'bd': 3.1988043553441288,'cd': 0.0,'Nub': 0.0,'Ndb': 0.0,'kt2': 0.16000000088572233}
    simulatesivers(sbs,par1)
    simulatesivers(clas,par1)
    simulatesivers(base,par1)
    simulatesivers(enhanced,par1)
    simulatesivers(base3he,par1)
    simulatesivers(enhanced3he,par1)
    sbsrep = sbs.copy()
    clasrep = clas.copy()
    baserep = base.copy()
    baserep['error'] = base['stat']
    basesystrep = base.copy()
    basesystrep['error'] = (base['stat']**2+base['systabs']**2+base['value']**2*base['systrel']**2)**0.5
    enhancedrep = enhanced.copy()
    enhancedrep['error'] = enhanced['stat']
    enhancedsystrep = enhanced.copy()
    enhancedsystrep['error'] = (enhanced['stat']**2+enhanced['systabs']**2+enhanced['value']**2*enhanced['systrel']**2)**0.5
    base3herep = base3he.copy()
    base3herep['error'] = base3he['stat']
    base3hesystrep = base3he.copy()
    base3hesystrep['error'] = (base3he['stat']**2+base3he['systabs']**2+base3he['value']**2*base3he['systrel']**2)**0.5
    enhanced3herep = enhanced3he.copy()
    enhanced3herep['error'] = enhanced3he['stat']
    enhanced3hesystrep = enhanced3he.copy()
    enhanced3hesystrep['error'] = (enhanced3he['stat']**2+enhanced3he['systabs']**2+enhanced3he['value']**2*enhanced3he['systrel']**2)**0.5
    sbsrep.to_csv('data/simsbs.dat', sep='\t', index=False)
    clasrep.to_csv('data/simclas.dat', sep='\t', index=False)
    baserep.to_csv('data/simbase.dat', sep='\t', index=False)
    basesystrep.to_csv('data/simbasesyst.dat', sep='\t', index=False)
    enhancedrep.to_csv('data/simenhanced.dat', sep='\t', index=False)
    enhancedsystrep.to_csv('data/simenhancedsyst.dat', sep='\t', index=False)
    base3herep.to_csv('data/simbase3he.dat', sep='\t', index=False)
    base3hesystrep.to_csv('data/simbase3hesyst.dat', sep='\t', index=False)
    enhanced3herep.to_csv('data/simenhanced3he.dat', sep='\t', index=False)
    enhanced3hesystrep.to_csv('data/simenhanced3hesyst.dat', sep='\t', index=False)
    return




if __name__ == "__main__":
    opt = sys.argv[1]
    if opt == 'collins':
        print("Preparing pseudodata sets for Collins asymmetry ...", end='\n')
        preparecollins()
    elif opt == 'sivers':
        print("Preparing pseudodata sets for Sivers asymmetry ...", end='\n')
        preparesivers()

    exit
