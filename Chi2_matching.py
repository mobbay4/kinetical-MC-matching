from basf2 import *
from ROOT import Belle2
import pandas as pd
import numpy.matlib
import numpy as np

#import time
import sys
import ROOT
import copy



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
        #print(dir(Belle2.Const.chargedStableSet))
        #print(Belle2.Const.chargedStableSet.find(11))        
        #print(Belle2.Const.chargedStableSet.find(12))

        #chargedstable list
        chargedstable_ids = [11,13,211,321,2212,1000010020]

        for track in self.tracks:

            #test ob track schon relation zu mc hat dann track überspringen
            # getRelated(mcPart) und auf null pointer testen 

            mccount = 0
            tracknum +=1
            chi2_l=[]
            chi2_n=[]
            for mcparticle in self.mcparticles:
                mcnum +=1
                abspdg = abs(mcparticle.getPDG())
                # check if current mc_particle dose traces
                if abspdg in chargedstable_ids:
                    # find most likely Trackfitresult
                    trackfit_wcm = track.getTrackFitResultWithClosestMass(Belle2.Const.ChargedStable(abspdg))
                    d =trackfit_wcm.getParticleType()
                    
                    if type(trackfit_wcm)!= type(Belle2.TrackFitResult()):

                        print("..................................................")
                        print("==========")
                        print("trackfit_wcm = ",trackfit_wcm)
                        print("type(trackfit_wcm) = ",type(trackfit_wcm))
                        print("==========")
                        break
                    # get trackfit_wcm variables
                    covariance = trackfit_wcm.getCovariance5()#uh.getCovariance()
                    cov_C = trackfit_wcm.getCovariance5()
                    covariance_inv = cov_C.Invert()
                    track_val = {} 
                    track_val["TanLambdaErr"]   = np.sqrt(covariance(4,4))
                    track_val["Z0Err"]          = np.sqrt(covariance(3,3))
                    track_val["OmegaErr"]       = np.sqrt(covariance(2,2))
                    track_val["Phi0Err"]        = np.sqrt(covariance(1,1))
                    track_val["D0Err"]          = np.sqrt(covariance(0,0))
                    track_val["TanLambda"]      = trackfit_wcm.getTanLambda()
                    track_val["Z0"]             = trackfit_wcm.getZ0()
                    track_val["Omega"]          = trackfit_wcm.getOmega()
                    track_val["Phi0"]           = trackfit_wcm.getPhi0()
                    track_val["D0"]             = trackfit_wcm.getD0()

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

                    

                    cov = np.matlib.zeros((5,5))
                    for i in range(0,5):
                        for j in range(0,5):
                            cov[i,j] =covariance(i,j)
                    try:
                        cov_inv=np.linalg.inv(cov)
                    except:
                        continue
                    #chi2
                    Delta= np.matlib.zeros((1,5))
                    Delta = Delta.T
                    Delta[4,0]=   track_val["TanLambda"]-par_val["TanLambda"]
                    Delta[3,0]=   track_val["Z0"]-par_val["Z0"]
                    Delta[2,0] =   track_val["Omega"]-par_val["Omega"]
                    Delta[1,0] =   track_val["Phi0"]-par_val["Phi0"]
                    Delta[0,0] =   track_val["D0"]-par_val["D0"]
                    chi2= np.matmul(Delta.T,np.matmul(cov_inv,Delta))
                    chi2_l.append(chi2[0,0])
                    chi2_n.append(mccount)
                    self.count +=1      
                    
                mccount+=1 # to get right order in self.mcparticles
            # set relation        
            if len(chi2_l)!=0:
                min_chi=min(chi2_l)
                if min_chi< float(CutOff):#Cut Off
                    track.addRelationTo(self.mcparticles[chi2_n[chi2_l.index(min_chi)]])
        print("eventnum = ",self.eventnumber)
  
                

#============================================
#***************Testcall*********************
#============================================

set_random_seed(12345)

# create path
main = create_path()

eventnum=sys.argv[3]
typ = sys.argv[2]
server=sys.argv[1]
if typ == "gun":
    CutOff=sys.argv[5]
    pdg=sys.argv[4]
    loadfile=f"gun_20000_5_w_all_pdg{pdg}.root"
    tag=f"chi2_match_pdg{pdg}_20000_gun"
if typ == "sim":
    CutOff=sys.argv[4]
    loadfile = f"sim_{eventnum}.root"
    tag = f"chi2_match_CutOff{CutOff}_{eventnum}_{typ}_newreco_D0cut_fillfromparticle"
if server=="kuhrios":
    basepath = "/srv/data/floschw/D/Data/modul_test/"
elif server=="naf":
    basepath = "~/D/Data/kuhrios/Data/modul_test/"


main.add_module("RootInput",inputFileName = f"{basepath}simulation/{loadfile}"
                           ,excludeBranchNames=["TracksToMCParticles"])


# add our module
main.add_module(mc_matching_new())

filename=f"{basepath}analysis/tuple_{tag}.root"

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
    
    
    #filename = f"~/D/Data/kuhrios/Data/modul_test/analysis/tuple_{tag}.root"
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

#process the thing
process(main)
#print(statistics)


print("******************************* end Module ***************************************")
