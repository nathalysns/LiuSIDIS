#!/usr/bin/env python

import sys
import numpy as np
import scipy as sp
import pandas as pd
from iminuit import Minuit
from scipy import random

import tmdlib.tmd as tmd

world = pd.read_csv('datacollins/colworld.dat', delim_whitespace=True)
sbs = pd.read_csv('datacollins/simsbs.dat', delim_whitespace=True)
clas = pd.read_csv('datacollins/simclas.dat', delim_whitespace=True)
base = pd.read_csv('datacollins/simbase.dat', delim_whitespace=True)
basesyst = pd.read_csv('datacollins/simbasesyst.dat', delim_whitespace=True)
enhanced = pd.read_csv('datacollins/simenhanced.dat', delim_whitespace=True)
enhancedsyst = pd.read_csv('datacollins/simenhancedsyst.dat', delim_whitespace=True)

worldrep = 0
simdata = 0
simdatarep = 0

def fitfunc0(var):
    par = {'Nu':var[0], 'Nd':var[1], 'a':var[2], 'b':var[3], 'c':var[4], 'kt2':var[5]}
    res = []
    for i in range(len(worldrep)):
        tmp = worldrep.loc[i]
        x, y, Q2, z, pT, target, hadron, value, error = tmp['x'], tmp['y'], tmp['Q2'], tmp['z'], tmp['pT'], tmp['target'], tmp['hadron'], tmp['value'], tmp['error']
        res.append( (tmd.AUTCollins(x,y,Q2,z,pT,target,hadron,par) - value)**2 / error**2)
    return np.sum(res)

def fitfunc(var):
    par = {'Nu':var[0], 'Nd':var[1], 'a':var[2], 'b':var[3], 'c':var[4], 'kt2':var[5]}
    res = []
    for i in range(len(worldrep)):
        tmp = worldrep.loc[i]
        x, y, Q2, z, pT, target, hadron, value, error = tmp['x'], tmp['y'], tmp['Q2'], tmp['z'], tmp['pT'], tmp['target'], tmp['hadron'], tmp['value'], tmp['error']
        res.append( (tmd.AUTCollins(x,y,Q2,z,pT,target,hadron,par) - value)**2 / error**2)
    for i in range(len(simdatarep)):
        tmp = simdatarep.loc[i]
        x, y, Q2, z, pT, target, hadron, value, error = tmp['x'], tmp['y'], tmp['Q2'], tmp['z'], tmp['pT'], tmp['target'], tmp['hadron'], tmp['value'], tmp['error']
        res.append( (tmd.AUTCollins(x,y,Q2,z,pT,target,hadron,par) - value)**2 / error**2)
    return np.sum(res)

def fitworld(Nrep, filename):
    global worldrep
    worldrep = world.copy()
    out = []
    var0 = [0.4, -0.45, 1.0, 3.0, 0.0, 0.25]
    for i in range(Nrep):
        print(i, end='\r')
        worldrep['value'] = random.normal(world['value'], world['error'])
        Min = Minuit.from_array_func(fitfunc0, start=var0,\
                name=['Nu','Nd','a','b','c','kt2'],\
                error=[1e-4,1e-4,1e-4,1e-4,1e-4,1e-4],\
                fix=[False,False,False,False,True,False],\
                limit=[None,None,None,None,None,None],\
                errordef=1)
        Min.print_level=0
        Min.strategy=1
        Min.migrad()
        out.append([Min.values['Nu'],Min.values['Nd'],Min.values['a'],Min.values['b'],Min.values['c'],Min.values['kt2'],Min.fval])
        del Min
    fs = pd.DataFrame(out, columns=['Nu','Nd','a','b','c','kt2','chi2'])
    fs.to_csv(filename, sep='\t', index=False)
    return

def simulate(data):
    global simdata, worldrep
    worldrep = world.copy()
    var0 = [0.4, -0.45, 1.0, 3.0, 0.0, 0.25]
    Min = Minuit.from_array_func(fitfunc0, start=var0,\
                name=['Nu','Nd','a','b','c','kt2'],\
                error=[1e-4,1e-4,1e-4,1e-4,1e-4,1e-4],\
                fix=[False,False,False,False,True,False],\
                limit=[None,None,None,None,None,None],\
                errordef=1)
    Min.print_level=0
    Min.strategy=1
    Min.migrad()
    par0 = Min.values
    var0 = Min.np_values()
    val = []
    for i in range(len(simdata)):
        tmp = simdata.loc[i]
        x, y, Q2, z, pT, target, hadron = tmp['x'], tmp['y'], tmp['Q2'], tmp['z'], tmp['pT'], tmp['target'], tmp['hadron']
        val.append(tmd.AUTCollins(x,y,Q2,z,pT,target,hadron,par0))
    simdata['value'] = val
    return var0

def fitsim(Nrep, filename):
    global world, worldrep, simdata, simdatarep
    var0 = simulate(simdata)
    simdatarep = simdata.copy()
    out = []
    for i in range(Nrep):
        print(i, end='\r')
        simdatarep['value'] = random.normal(simdata['value'], simdata['error'])
        Min = Minuit.from_array_func(fitfunc, start=var0,\
                name=['Nu','Nd','a','b','c','kt2'],\
                error=[1e-4,1e-4,1e-4,1e-4,1e-4,1e-4],\
                fix=[False,False,False,False,False,False],\
                limit=[None,None,None,None,None,None],\
                errordef=1)
        Min.print_level=0
        Min.strategy=1
        Min.migrad()
        out.append([Min.values['Nu'],Min.values['Nd'],Min.values['a'],Min.values['b'],Min.values['c'],Min.values['kt2'],Min.fval])
        del Min
    fs = pd.DataFrame(out, columns=['Nu','Nd','a','b','c','kt2','chi2'])
    fs.to_csv(filename, sep='\t', index=False)

if __name__ == "__main__":
    print('Running the fit to transversity function ...', end='\n')
    opt = sys.argv[1]
    if opt == 'world':
        print('fitting world data ...', end='\n')
        fitworld(50, 'outputcollins/out-world.dat')
    elif opt == 'sbs':
        print('fitting SBS ...', end='\n')
        simdata = sbs.copy()
        fitsim(50, 'outputcollins/out-sbs.dat')
    elif opt == 'clas':
        print('fitting CLAS12 ...', end='\n')
        simdata = clas.copy()
        fitsim(50, 'outputcollins/out-clas.dat')
    elif opt == 'sbs+clas':
        print('fitting SBS+CLAS12 ...', end='\n')
        simdata = pd.concat([sbs,clas], axis=0, ignore_index=True)
        fitsim(50, 'outputcollins/out-sbsclas.dat')
    elif opt == 'base':
        print('fitting SoLID baseline ...', end='\n')
        simdata = base.copy()
        fitsim(50, 'outputcollins/out-base.dat')
    elif opt == 'basesyst':
        print('fitting SoLID baseline (including syst) ...', end='\n')
        simdata = basesyst.copy()
        fitsim(50, 'outputcollins/out-basesyst.dat')
    elif opt == 'enhanced':
        print('fitting SoLID enhanced ...', end='\n')
        simdata = enhanced.copy()
        fitsim(50, 'outputcollins/out-enhanced.dat')
    elif opt == 'enhancedsyst':
        print('fitting SoLID enhanced (including syst) ...', end='\n')
        simdata = enhancedsyst.copy()
        fitsim(50, 'outputcollins/out-enhancedsyst.dat')
    elif opt == 'sbs+clas+base':
        print('fitting SBS+CLAS12+SoLID baseline ...', end='\n')
        simdata = pd.concat([sbs,clas,base], axis=0, ignore_index=True)
        fitsim(50, 'outputcollins/out-sbsclasbase.dat')
    elif opt == 'sbs+clas+basesyst':
        print('fitting SBS+CLAS12+SoLID baseline (including syst) ...', end='\n')
        simdata = pd.concat([sbs,clas,basesyst], axis=0, ignore_index=True)
        fitsim(50, 'outputcollins/out-sbsclasbasesyst.dat')
    elif opt == 'sbs+clas+enhanced':
        print('fitting SBS+CLAS12+SoLID enhanced ...', end='\n')
        simdata = pd.concat([sbs,clas,enhanced], axis=0, ignore_index=True)
        fitsim(50, 'outputcollins/out-sbsclasenhanced.dat')
    elif opt == 'sbs+clas+enhancedsyst':
        print('fitting SBS+CLAS12+SoLID enhancedsyst (including syst) ...', end='\n')
        simdata = pd.concat([sbs,clas,enhancedsyst], axis=0, ignore_index=True)
        fitsim(50, 'outputcollins/out-sbsclasenhancedsyst.dat')

    exit
    
    
    
    
        
        
