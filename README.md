# LiuSIDIS
A sidis generator

# Authors
Tianbo Liu, Zhiwen Zhao, Weizhi Xiong, Haiyan Gao

# How to use it
A quick start to use this generator
## Initialization
* __#include "Lsidis.h"__ //including the head file
* __Lsidis mysidis;__     //create an object
* __TLorentzVector l;__   //initial state electron
* __TLorentzVector P;__   //initial state nucleon
* __mysidis.SetNucleus(Np, Nn);__  //set proton and neutron numbers 
* __mysidis.SetHadron("pi+");__    //set the final state hadron
* __mysidis.SetInitialState(l,P);__ //set initial state kinematics

