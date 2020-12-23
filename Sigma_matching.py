from basf2 import *
from ROOT import Belle2
import pandas as pd
import numpy as np
#import time
import sys

"""
To Do:
------
- implement chargedstabelset 
"""
print("******************************* start Module ***************************************")

class mc_matching_new(Module):
    def __init__(self):
        super().__init__()
    def initialize(self):
        
        # get requiered classes and check if they exist
        self.mcparticles = Belle2.PyStoreArray("MCParticles")
        self.mcparticles.isRequired()
        
        self.tracks = Belle2.PyStoreArray("Tracks")
        self.tracks.isRequired()
        self.tracks.registerRelationTo(self.mcparticles)

    def event(self): 
        #print("#entry_mc = ",self.mcparticles.getEntries())
        #print("#entry_track = ",self.tracks.getEntries())
        #print(dir(Belle2.Const.chargedStableSet))
        #print(Belle2.Const.chargedStableSet.find(11))        
        #print(Belle2.Const.chargedStableSet.find(12))
        
        #chargedstable list
        chargedstable_ids = [11,13,211,321,2212,1000010020]

        
        
        for track in self.tracks:
            controlsignal = 0
            for mcparticle in self.mcparticles:
                abspdg = abs(mcparticle.getPDG())
                # check if current mc_particle dose traces
                if abspdg in chargedstable_ids:
                    
                    # find most likely Trackfitresult
                    trackfit_wcm = track.getTrackFitResultWithClosestMass(Belle2.Const.ChargedStable(abspdg))
                    if type(trackfit_wcm)!= type(Belle2.TrackFitResult()):
                        print("..................................................")
                        print("==========")
                        print("trackfit_wcm = ",trackfit_wcm)
                        print("type(trackfit_wcm) = ",type(trackfit_wcm))
                        print("==========")
                        break
                    # get trackfit_wcm variables
                    uh = trackfit_wcm.getUncertainHelix()
                    convariance = uh.getCovariance()

                    track_val = {}
                    track_val["TanLambdaErr"]   = np.sqrt(convariance(4,4))
                    track_val["Z0Err"]          = np.sqrt(convariance(3,3))
                    track_val["OmegaErr"]       = np.sqrt(convariance(2,2))
                    track_val["Phi0Err"]        = np.sqrt(convariance(1,1))
                    track_val["D0Err"]          = np.sqrt(convariance(0,0))
                    track_val["TanLambda"]      = uh.getTanLambda()
                    track_val["Z0"]             = uh.getZ0()
                    track_val["Omega"]          = uh.getOmega()
                    track_val["Phi0"]           = uh.getPhi0()
                    track_val["D0"]             = uh.getD0() 

                    # get mc mcparticle variables
                    charge_sign = (-1 if mcparticle.getCharge() < 0 else 1)
                    b_field = Belle2.BFieldManager.getField(mcparticle.getVertex()).Z() / Belle2.Unit.T
                    h = Belle2.Helix(mcparticle.getVertex(), mcparticle.getMomentum(), charge_sign, b_field)
                    
                    par_val = {}
                    par_val["TanLambda"] = h.getTanLambda()
                    par_val["Z0"]        = h.getZ0()
                    par_val["Omega"]     = h.getOmega()
                    par_val["Phi0"]      = h.getPhi0()
                    par_val["D0"]        = h.getD0()
                    

                    #if track_val["Phi0"]<0:
                    #    track_val["Phi0"]=track_val["Phi0"]+2*np.pi
                    #if par_val["Phi0"]<0:
                    #    par_val["Phi0"]=par_val["Phi0"]+2*np.pi

                    # compare the variables and if all fullfil the 3 Sigma interval we call it a signal
                    overall = 1
                    for j in ["Z0","D0","Omega","TanLambda","Phi0"]:
                        if not abs(track_val[j]-par_val[j]) < int(num_sigma)*abs(track_val[(j+"Err")]):
                            overall=0                            
                            break
                        
                    if overall==1:
                        if controlsignal<1:
                            track.addRelationTo(mcparticle)
                        controlsignal+=1
            

#============================================
#***************Testcall*********************
#============================================
##
num_sigma = sys.argv[4]
eventnum  = str(sys.argv[3]) 
typ       = sys.argv[2]
server    = sys.argv[1]

##
 
set_random_seed(12345)
main = create_path()
if typ == "sim":
    loadfile = f"sim_{eventnum}.root"
    tag = f"sigma_{num_sigma}_match_{eventnum}_{typ}_newreco_D0cut_fillfromparticle"
if typ == "gun":
    pdg=sys.argv[5]
    loadfile = f"gun_20000_5_w_all_pdg{pdg}.root"
    eventnum =str(20000)
    tag = f"sigma_{num_sigma}_match_{typ}_{eventnum}_pdg{pdg}"
if server=="kuhrios":
    basepath = "/srv/data/floschw/D/Data/modul_test/"
elif server=="naf":
    basepath = "~/D/Data/kuhrios/Data/modul_test/"

main.add_module("RootInput",inputFileName =f"{basepath}simulation/{loadfile}",excludeBranchNames=["TracksToMCParticles"])
main.add_module(mc_matching_new())

##
filename=f"{basepath}analysis/tuple_{tag}.root"
##

if typ=="sim":
    # example analysis D*+:pi reconstruction
    ##########################################################
    from stdCharged import stdK, stdPi, stdE, stdMu
    import modularAnalysis as ma
    import variables.collections as vc
    import variables.utils as vu

    stdPi("all",path=main)
    stdK("good",path=main)
    decaystring = "anti-B0 -> [D*+ -> [[D0 -> K-:good pi+:good] pi+:good]e-:good anti-nu_e"
    ma.reconstructDecay("D0 -> K-:good pi+:all", "1.85466<M<1.875", path=main)
    ma.matchMCTruth("D0", path=main)

    ma.reconstructDecay("D*+ -> D0 pi+:all", "", path=main)
    ma.matchMCTruth("D*+", path=main)

    particle= "D*+"
    listofvariables  = ["E","M", "PDG","isSignalAcceptMissing"] 
    listofvariables += ["d0","d0Err","omega","omegaErr","phi0","phi0Err","phi","phiErr","tanlambda","tanlambdaErr","z0","z0Err","theta","thetaErr",
                        "p","pErr","pt","ptErr","px","pxErr","py","pyErr", "pz","pzErr",
                        "d0Pull","omegaPull","phi0Pull","tanLambdaPull","z0Pull"]
    listofvariables += ["mcP","mcPT","mcPX","mcPY","mcPZ","mcDecayVertexX","mcDecayVertexY","mcDecayVertexZ","isSignal","mcErrors",
                        "mcPDG","charge"]  + vc.mc_vertex
    variables = vu.create_aliases(listofvariables,"{variable}","Dstar")
    variables += vu.create_aliases(listofvariables,"daughter(0,{variable})","D0")
    variables += vu.create_aliases(listofvariables,"daughter(1,{variable})","Dstarpi")

    ma.variablesToNtuple(particle, variables = variables, filename = filename , treename = "tree", path=main)

    ma.fillParticleListFromMC("pi+:fun",cut="",skipNonPrimaryDaughters=True,path= main)
    particle2 = "pi+:fun"
    variables2 = vu.create_aliases(listofvariables,"{variable}","fun")
    ma.variablesToNtuple(particle2, variables = variables2, filename = filename , treename = "fun", path=main)
    #############################################################

if typ=="gun":
        # example analysis for particle gun
    ##########################################################
    from stdCharged import stdK, stdPi, stdE, stdMu
    import stdCharged
    import modularAnalysis as ma
    import variables.collections as vc
    import variables.utils as vu
    listofvariables  = ["E","M", "PDG","isSignalAcceptMissing"] 
    listofvariables += ["d0","d0Err","omega","omegaErr","phi0","phi0Err","phi","phiErr","tanlambda","tanlambdaErr","z0","z0Err","theta","thetaErr",
                        "p","pErr","pt","ptErr","px","pxErr","py","pyErr", "pz","pzErr",
                        "dx","dy","dz",
                        "d0Pull","omegaPull","phi0Pull","tanLambdaPull","z0Pull"]
    listofvariables += ["mcP","mcPT","mcPX","mcPY","mcPZ","isSignal","mcErrors", "mcPDG","charge",
                        "mcDecayVertexFromIPX","mcDecayVertexFromIPY","mcDecayVertexFromIPZ","nTracks"] + vc.mc_vertex
        
    if int(pdg)==211:#part=="":
        stdPi("all",path=main)

        ma.fillParticleListFromMC("pi+:fun",cut="",skipNonPrimaryDaughters=True,path= main)

        particle= "pi+:all" #sys.argv[2] 
        variables = vu.create_aliases(listofvariables,"{variable}","all")
        particle2 = "pi+:fun"
        variables2 = vu.create_aliases(listofvariables,"{variable}","fun")
      
        ma.variablesToNtuple(particle, variables = variables, filename = filename , treename = "all", path=main)
        ma.variablesToNtuple(particle2, variables = variables2, filename = filename , treename = "fun", path=main)
    elif int(pdg)==11:#part=="_electron":
        stdCharged.stdE("all",path=main)

        particle= "e-:all" 
        
        variables = vu.create_aliases(listofvariables,"{variable}","all")
        ma.variablesToNtuple(particle, variables = variables, filename = filename , treename = "all", path=main)
    elif int(pdg)==13:#part=="_myon":
        stdCharged.stdMu("all",path=main)

        particle= "mu-:all"  
        
        variables = vu.create_aliases(listofvariables,"{variable}","all")

        ma.variablesToNtuple(particle, variables = variables, filename = filename , treename = "all", path=main)
    elif int(pdg)==2212:#part=="_proton":
        stdCharged.stdPr("all",path=main)

        particle= "p+:all" 
        variables = vu.create_aliases(listofvariables,"{variable}","all")

        ma.variablesToNtuple(particle, variables = variables, filename = filename , treename = "all", path=main)
    #############################################################

process(main)

print("******************************* end Module ***************************************")
