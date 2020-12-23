#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# to generate root data files to test my module
print("******************************* start root Module ***************************************")

from basf2 import *
from simulation import add_simulation
from L1trigger import add_tsim
from reconstruction import add_reconstruction
import sys
from mdst import add_mdst_output

eventnum = sys.argv[1]
pdg=sys.argv[3]

set_random_seed(12345)

# create path
main = create_path()
# specify number of events to be generated
main.add_module('EventInfoSetter', evtNumList=[int(eventnum)])

particleGunModule = register_module('ParticleGun')
particleGunModule.param({
    'pdgCodes': [int(pdg)],#[211],
    'nTracks': 1,
    'varyNTracks': False,
    'momentumGeneration': 'uniformPt',
    'momentumParams': [0.05, float(sys.argv[2])],#0.05
    'thetaGeneration': 'uniform',
    'thetaParams': [17., 150.]
})

main.add_module(particleGunModule)

# detector simulation
add_simulation(main)
# trigger simulation
add_tsim(main) 
# reconstruction
add_reconstruction(main)


# full output
p_up = ""
for i in str(sys.argv[2]):
    if i !=".":
        p_up+=i
    else:
        p_up+="_"
tag = "gun_" + str(eventnum) +"_"+p_up+"_w_all"+f"_pdg{pdg}"#+"_pdown4"+"_electron"_vert _looper +"pdown2
outputname = f"{tag}.root"
add_mdst_output(main,filename=str(sys.argv[4])+"/D/Data/modul_test/simulation/"+outputname)

process(main)
#print(statistics)
print("******************************* end root Module ***************************************")
