#!/usr/bin/env python

import sys
import numpy as np
import scipy as sp
import pandas as pd
from iminuit import Minuit
from scipy import random

import tmdlib.tmd as tmd

world = pd.read_csv('data/colworld.dat', delim_whitespace=True)
sbs = pd.read_csv('data/simsbs.dat', delim_whitespace=True)
clas = pd.read_csv('data/simclas.dat', delim_whitespace=True)
base = pd.read_csv('data/simbase.dat', delim_whitespace=True)
basesyst = pd.read_csv('data/simbasesyst.dat', delim_whitespace=True)
enhanced = pd.read_csv('data/simenhanced.dat', delim_whitespace=True)
enhancedsyst = pd.read_csv('data/simenhancedsyst.dat', delim_whitespace=True)
base3he = pd.read_csv('data/simbase3he.dat', delim_whitespace=True)
base3hesyst = pd.read_csv('data/simbase3hesyst.dat', delim_whitespace=True)
enhanced3he = pd.read_csv('data/simenhanced3he.dat', delim_whitespace=True)
enhanced3hesyst = pd.read_csv('data/simenhanced3hesyst.dat', delim_whitespace=True)


worldrep = 0
simdata = 0
simdatarep = 0

def fitfunc0(var):
    par = {'Nu':var[0], 'au':var[1], 'bu':var[2], 'cu':var[3],\
           'Nd':var[4], 'ad':var[5], 'bd':var[6], 'cd':var[7],\
           'Nub':var[8], 'Ndb':var[9], 'kt2':var[10]}
    res = []
    for i in range(len(worldrep)):
        x, Q2, z, pT, target, hadron, value, error = worldrep.loc[i]['x'], worldrep.loc[i]['Q2'], worldrep.loc[i]['z'], worldrep.loc[i]['pT'], worldrep.loc[i]['target'], worldrep.loc[i]['hadron'], worldrep.loc[i]['value'], worldrep.loc[i]['error']
        res.append( (tmd.AUTSivers(x,Q2,z,pT,target,hadron,par) - value)**2 / error**2 )
    return np.sum(res)

def fitfunc(var):
    par = {'Nu':var[0], 'au':var[1], 'bu':var[2], 'cu':var[3],\
           'Nd':var[4], 'ad':var[5], 'bd':var[6], 'cd':var[7],\
           'Nub':var[8], 'Ndb':var[9], 'kt2':var[10]}
    res = []
    for i in range(len(worldrep)):
        x, Q2, z, pT, target, hadron, value, error = worldrep.loc[i]['x'], worldrep.loc[i]['Q2'], worldrep.loc[i]['z'], worldrep.loc[i]['pT'], worldrep.loc[i]['target'], worldrep.loc[i]['hadron'], worldrep.loc[i]['value'], worldrep.loc[i]['error']
        res.append( (tmd.AUTSivers(x,Q2,z,pT,target,hadron,par) - value)**2 / error**2 )
    for i in range(len(simdatarep)):
        x, Q2, z, pT, target, hadron, value, error = simdatarep.loc[i]['x'], simdatarep.loc[i]['Q2'], simdatarep.loc[i]['z'], simdatarep.loc[i]['pT'], simdatarep.loc[i]['target'], simdatarep.loc[i]['hadron'], simdatarep.loc[i]['value'], simdatarep.loc[i]['error']
        res.append( (tmd.AUTSivers(x,Q2,z,pT,target,hadron,par) - value)**2 / error**2 )
    return np.sum(res)
    

def fitworld(Nrep, filename):
    global worldrep
    worldrep = world.copy()
    out = []
    var0 = [-0.03851697, 0.66248284, 4.10308136,  0.0, 0.05141101, 0.39843212, 3.19880436, 0.0, 0.0, 0.0, 0.16]
    for i in range(Nrep):
        print(i, end='\r')
        worldrep['value'] = random.normal(world['value'], world['error'])
        Min = Minuit.from_array_func(fitfunc0, start=var0,\
                name=['Nu','au','bu','cu','Nd','ad','bd','cd','Nub','Ndb','kt2'],\
                error=[1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4],\
                fix=[False,False,False,True,False,False,False,True,True,True,False],\
                limit=[(-1,1),(0,3),(0,10),None,(-1,1),(0,3),(0,10),None,(-1,1),(-1,1),(0.1,1)],\
                errordef=1)
        Min.print_level=0
        Min.strategy=1
        Min.migrad()
        out.append([Min.values['Nu'],Min.values['au'],Min.values['bu'],Min.values['cu'],\
                Min.values['Nd'],Min.values['ad'],Min.values['bd'],Min.values['cd'],\
                Min.values['Nub'],Min.values['Ndb'],Min.values['kt2'],Min.fval])
        del Min
    fs = pd.DataFrame(out, columns=['Nu','au','bu','cu','Nd','ad','bd','cd','Nub','Ndb','kt2','chi2'])
    fs.to_csv(filename, sep='\t', index=False)
    return

def simulate(data):
    global simdata, worldrep
    worldrep = world.copy()
    var0 = [-0.03851697, 0.66248284, 4.10308136,  0.0, 0.05141101, 0.39843212, 3.19880436, 0.0, 0.0, 0.0, 0.16]
    Min = Minuit.from_array_func(fitfunc0, start=var0,\
                name=['Nu','au','bu','cu','Nd','ad','bd','cd','Nub','Ndb','kt2'],\
                error=[1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4],\
                fix=[False,False,False,True,False,False,False,True,True,True,False],\
                limit=[(-1,1),(0,3),(0,10),None,(-1,1),(0,3),(0,10),None,(-1,1),(-1,1),(0.1,1)],\
                errordef=1)
    Min.print_level=0
    Min.strategy=1
    Min.migrad()
    par0 = Min.values
    var0 = Min.np_values()
    val = []
    for i in range(len(simdata)):
        x, y, z, Q2, pT, target, hadron = simdata.loc[i]['x'], simdata.loc[i]['y'], simdata.loc[i]['z'], simdata.loc[i]['Q2'], simdata.loc[i]['pT'], simdata.loc[i]['target'], simdata.loc[i]['hadron']
        val.append(tmd.AUTSivers(x, Q2, z, pT, target, hadron, par0))
    simdata['value'] = val
    return var0
        
def fitsim(Nrep, filename):
    global world, worldrep, simdata, simdatarep
    var0 = simulate(simdata)
    simdatarep = simdata.copy()
    out=[]
    for i in range(Nrep):
        print(i, end='\r')
        simdatarep['value'] = random.normal(simdata['value'], simdata['error'])
        Min = Minuit.from_array_func(fitfunc, start=var0,\
                name=['Nu','au','bu','cu','Nd','ad','bd','cd','Nub','Ndb','kt2'],\
                error=[1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4],\
                fix=[False,False,False,False,False,False,False,False,True,True,False],\
                limit=[(-1,1),(0,3),(0,10),None,(-1,1),(0,3),(0,10),None,(-1,1),(-1,1),(0.1,1)],\
                errordef=1)
        Min.print_level=0
        Min.strategy=1
        Min.migrad()
        out.append([Min.values['Nu'],Min.values['au'],Min.values['bu'],Min.values['cu'],\
                Min.values['Nd'],Min.values['ad'],Min.values['bd'],Min.values['cd'],\
                Min.values['Nub'],Min.values['Ndb'],Min.values['kt2'],Min.fval])
        del Min
    fs = pd.DataFrame(out, columns=['Nu','au','bu','cu','Nd','ad','bd','cd','Nub','Ndb','kt2','chi2'])
    fs.to_csv(filename, sep='\t', index=False)
    return


if __name__ == "__main__":
    print('Running the fit to Sivers function ...', end = '\n')
    opt = sys.argv[1]
    if opt == 'world':
        print('fitting world data ...', end='\n')
        fitworld(50, 'output/out-world.dat')
    elif opt == 'sbs':
        print('fitting SBS ...', end='\n')
        simdata = sbs.copy()
        fitsim(30, 'output/out-sbs.dat')
    elif opt == 'clas':
        print('fitting CLAS12 ...', end='\n')
        simdata = clas.copy()
        fitsim(30, 'output/out-clas.dat')
    elif opt == 'sbs+clas':
        print('fitting SBS+CLAS12 ...', end='\n')
        simdata = pd.concat([sbs,clas], axis=0, ignore_index=True)
        fitsim(30, 'output/out-sbsclas.dat')
    elif opt == 'base':
        print('fitting SoLID baseline ...', end='\n')
        simdata = base.copy()
        fitsim(30, 'output/out-base.dat')
    elif opt == 'basesyst':
        print('fitting SoLID baseline (including syst) ...', end='\n')
        simdata = basesyst.copy()
        fitsim(30, 'output/out-basesyst.dat')
    elif opt == 'enhanced':
        print('fitting SoLID enhanced ...', end='\n')
        simdata = enhanced.copy()
        fitsim(30, 'output/out-enhanced.dat')
    elif opt == 'enhancedsyst':
        print('fitting SoLID enhanced (including syst) ...', end='\n')
        simdata = enhancedsyst.copy()
        fitsim(30, 'output/out-enhancedsyst.dat')
    elif opt == 'sbs+clas+base':
        print('fitting SBD+CLAS12+SoLID baseline ...', end='\n')
        simdata = pd.concat([sbs,clas,base], axis=0, ignore_index=True)
        fitsim(30, 'output/out-sbsclasbase.dat')
    elif opt == 'sbs+clas+basesyst':
        print('fitting SBD+CLAS12+SoLID baseline (including syst) ...', end='\n')
        simdata = pd.concat([sbs,clas,basesyst], axis=0, ignore_index=True)
        fitsim(30, 'output/out-sbsclasbasesyst.dat')
    elif opt == 'sbs+clas+enhanced':
        print('fitting SBD+CLAS12+SoLID enhanced ...', end='\n')
        simdata = pd.concat([sbs,clas,enhanced], axis=0, ignore_index=True)
        fitsim(30, 'output/out-sbsclasenhanced.dat')
    elif opt == 'sbs+clas+enhancedsyst':
        print('fitting SBD+CLAS12+SoLID enhanced (including syst) ...', end='\n')
        simdata = pd.concat([sbs,clas,enhancedsyst], axis=0, ignore_index=True)
        fitsim(30, 'output/out-sbsclasenhancedsyst.dat')
    elif opt == 'base3he':
        print('fitting SoLID baseline 3he ...', end='\n')
        simdata = base3he.copy()
        fitsim(30, 'output/out-base3he.dat')
    elif opt == 'base3hesyst':
        print('fitting SoLID baseline 3he (including syst) ...', end='\n')
        simdata = base3hesyst.copy()
        fitsim(30, 'output/out-base3hesyst.dat')
    elif opt == 'enhanced3he':
        print('fitting SoLID enhanced 3he ...', end='\n')
        simdata = enhanced3he.copy()
        fitsim(30, 'output/out-enhanced3he.dat')
    elif opt == 'enhanced3hesyst':
        print('fitting SoLID enhanced 3he (including syst) ...', end='\n')
        simdata = enhanced3hesyst.copy()
        fitsim(30, 'output/out-enhanced3hesyst.dat')
    exit
