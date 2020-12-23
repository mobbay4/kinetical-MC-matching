# to generate root data files to test my module
print("******************************* start root Module ***************************************")

from basf2 import *
from simulation import add_simulation
from L1trigger import add_tsim
from reconstruction import add_reconstruction
import sys
from mdst import add_mdst_output

eventnum = sys.argv[1]

set_random_seed(12345)


# full output
if len(sys.argv)==3:
    basepath = sys.argv[2]
else:
    basepath = "~"


# create path
main = create_path()
# specify number of events to be generated
main.add_module('EventInfoSetter', evtNumList=[int(eventnum)])
# generate BBbar events
main.add_module('EvtGenInput'
                ,userDECFile=basepath+"/D/functions/Python_Modul/Dstar.dec"#/afs/desy.de/user/f/floschw/
                )
# detector simulation
add_simulation(main)
# trigger simulation
add_tsim(main)
# reconstruction
add_reconstruction(main)


tag = "sim_" + (eventnum)
outputname = f"{tag}.root"
add_mdst_output(main,filename=basepath+"/D/Data/modul_test/simulation/"+outputname)

process(main)
#print(statistics)
print("******************************* end root Module ***************************************")
