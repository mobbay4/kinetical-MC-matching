from basf2 import *
from ROOT import Belle2
import pandas as pd
import numpy.matlib
import numpy as np#import time
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

        # control variables
        self.controlsignal_all = 0
        self.notsofull_all = 0
        self.morethanonerelation = 0

        self.eventnumber = 0   
        self.count = 0
        self.df = pd.DataFrame(columns=["TanLambda_track","Z0_track","Omega_track","Phi0_track","D0_track",
                                        "TanLambda_trackErr","Z0_trackErr","Omega_trackErr","Phi0_trackErr","D0_trackErr",
                                        "TanLambda_mc","Z0_mc","Omega_mc","Phi0_mc","D0_mc",
                                        "TanLambda_rej","Z0_rej","Omega_rej","Phi0_rej","D0_rej",
                                        "eventnum","tracknum","mcnum","isRelation","setRelation","cova"
                                        ])
        self.norelationscount=0
    def event(self): 
        self.eventnumber +=1
        #print("#entry_mc = ",self.mcparticles.getEntries())
        #print("#entry_track = ",self.tracks.getEntries())
        #print(dir(Belle2.Const.chargedStableSet))
        #print(Belle2.Const.chargedStableSet.find(11))        
        #print(Belle2.Const.chargedStableSet.find(12))

        #chargedstable list
        chargedstable_ids = [11,13,211,321,2212,1000010020]
        tracknum=0
        mcnum=0

        
        for track in self.tracks:
            controlsignal = 0
            notsignal = 0
            tracknum +=1
            #print("relationsto = ", track.getRelationsTo("MCParticles"))
            #print("relationstotype = ",type(track.getRelationsTo("MCParticles")))
            #print("relationstodir = ",dir(track.getRelationsTo("MCParticles")))
            #print("relationstorelations = ",(track.getRelationsTo("MCParticles")).relations())
            #print("relationstosize = ",(track.getRelationsTo("MCParticles")).size())
            #print("relationstolänge = ",len(track.getRelationsTo("MCParticles")))

            if track.getRelationsTo("MCParticles").size()==0:
                self.norelationscount+=1
                chi2_l=[]
                chi2_n=[]
                mccount = 0
                breakcheck=True
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
                        try:
                            # get trackfit_wcm variables
                            covariance = trackfit_wcm.getCovariance5()#uh.getCovariance()
                            cov_C = trackfit_wcm.getCovariance5()
                            covariance_inv = cov_C.Invert()
                            #print(covariance)
                            #print(covariance_inv)
                            track_val = {} 
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

                            #if track_val["Phi0"]<0:
                            #    track_val["Phi0"]=track_val["Phi0"]+2*np.pi
                            #if par_val["Phi0"]<0:
                            #    par_val["Phi0"]=par_val["Phi0"]+2*np.pi
                            


                            cov = np.matlib.zeros((5,5))
                            for i in range(0,5):
                                for j in range(0,5):
                                    cov[i,j] =covariance(i,j)
                            #print(cov)
                            cov_inv=np.linalg.inv(cov)

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
                        except:
                            print("break appeared")
                            breakcheck=False
                            break
                            

                    mccount+=1
                # set relation
                if breakcheck:
                    min_chi=min(chi2_l)
                    if min_chi< 20:
                        
                        track.addRelationTo(self.mcparticles[chi2_n[chi2_l.index(min_chi)]])
                    else:
                        print("=================")
                        print("cutoff")
                        print("=================")
                
        print("eventnum = ",self.eventnumber)
        print("norelationscount =",self.norelationscount)
        """
                for mcparticle in self.mcparticles:
                    mcnum +=1
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

                        #print("======================")
                        #convariance.Print()
                        #print("======================")

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


                        self.df.loc[self.count,"eventnum"] = self.eventnumber
                        self.df.loc[self.count,"tracknum"] = tracknum
                        self.df.loc[self.count,"mcnum"]    = mcnum
                        df_cov = pd.DataFrame(columns=["0","1","2","3","4"])
                        cov={}
                        for i in range(0,5):
                            for j in range(0,5):
                                cov[str(i)+str(j)] =convariance(i,j) 
                                
                        self.df.at[self.count,"cova"]=cov
                        # compare the variables and if all fullfil the 3 Sigma interval we call it a signal
                        overall = 1
                        for j in ["Z0","D0","Omega","TanLambda","Phi0"]:
                            self.df.loc[self.count,j+"_track"] = track_val[j]
                            self.df.loc[self.count,j+"_mc"]    = par_val[j]
                            self.df.loc[self.count,j+"_trackErr"] = track_val[j+"Err"]
                            # check if the mcparticle variable is in a 3 sigma interval of the trackfit variable
                            if not abs(track_val[j]-par_val[j]) < int(sys.argv[2])*abs(track_val[(j+"Err")]):
                                overall=0
                                self.df.loc[self.count,j+"_rej"]=1
                            else:
                                self.df.loc[self.count,j+"_rej"]=0
                                
                                #break
                            
                        if overall==1:
                            self.df.loc[self.count,"isRelation"] = 1
                            self.controlsignal_all +=1
                            if controlsignal<1:
                                self.df.loc[self.count,"setRelation"] = 1
                                track.addRelationTo(mcparticle)
                            else:
                                self.df.loc[self.count,"setRelation"] = 0
                            controlsignal+=1
                        else:
                            self.df.loc[self.count,"isRelation"] = 0
                            self.notsofull_all +=1
                            notsignal +=1
                        self.df.loc[self.count,"controlsignal"] = controlsignal
                        self.count +=1
                
                if controlsignal>1:
                    self.morethanonerelation+=1
        print("eventnum = ",self.eventnumber)
        """
        

    def terminate(self):
        if len(sys.argv)==5:
            basepath = sys.argv[4]
        else:
            basepath = "~"
        tag = str(sys.argv[1]) + "_" + str(sys.argv[2]) + "_" + str(sys.argv[3])+"_hittune"+"_pdown4"+f"{part}"#+"_phi0correct"#+"_helixreverse"
        filename = f"{basepath}/D/Data/modul_test/analysis/df_{tag}.csv"
        #self.df.to_csv(filename)
        print("..................................................")
        print("==========")
        print("controlsignal = ", self.controlsignal_all) # how much signals
        print("==========")                    
        print("notsofull = ",self.notsofull_all) # how much posebilitys
        print("==========")
        #print("Signal/Background = ", self.controlsignal_all / self.notsofull_all)
        print("==========")
        print("==========")
        print("#morethanonerelation = ",self.morethanonerelation)
        print("==========")

#============================================
#***************Testcall*********************
#============================================

set_random_seed(12345)

# create path
main = create_path()
# load test data

# load test data
if len(sys.argv)==5:
    basepath = sys.argv[4]
else:
    basepath = "~"

part=""#sys.argv[6]
print(part)
gentype = str(sys.argv[3])
eventnum = str(sys.argv[1])

tag = gentype + "_" + eventnum #sys.argv[1] = eventnumber, 2 Sigma, 3 sim/gun
if gentype=="gun":
    p_up = ""
    for i in str(sys.argv[5]):
        if i !=".":
            p_up+=i
        else:
            p_up+="_"
    tag+="_"+p_up+"_w_all"+"_pdown4"+f"{part}"#+"_elektron"#+"_pdown3"
outputname = f"{tag}.root"
main.add_module("RootInput",inputFileName = "~/D/Data/vergleichsdaten/mdst_All_Dstar_04.root"#basepath+"/D/Data/modul_test/simulation/"+outputname#
                           )#,excludeBranchNames=["TracksToMCParticles"])
num_sigma = str(sys.argv[2])

# add our module
main.add_module(mc_matching_new())

if sys.argv[3] == "sim":
    # example analysis D*+:pi reconstruction
    ##########################################################
    from stdCharged import stdK, stdPi, stdE, stdMu
    import modularAnalysis as ma
    import variables.collections as vc
    import variables.utils as vu

    stdPi("good",path=main)
    stdPi("all",path=main)
    stdK("good",path=main)
    stdE("good",path=main)
    decaystring = "anti-B0 -> [D*+ -> [[D0 -> K-:good pi+:good] pi+:good]e-:good anti-nu_e"
    ma.reconstructDecay("D0 -> K-:good pi+:all", "", path=main)
    ma.matchMCTruth("D0", path=main)

    ma.reconstructDecay("D*+ -> D0 pi+:all", "", path=main)
    ma.matchMCTruth("D*+", path=main)

    ma.reconstructDecay("anti-B0 -> D*+ e-:good ", "", path=main) #anti-nu_e
    ma.matchMCTruth("anti-B0", path=main)

    particle= "anti-B0"
    listofvariables  = ["E","M", "PDG","isSignalAcceptMissing"] 
    listofvariables += ["d0","d0Err","omega","omegaErr","phi0","phi0Err","phi","phiErr","tanlambda","tanlambdaErr","z0","z0Err","theta","thetaErr",
                        "p","pErr","pt","ptErr","px","pxErr","py","pyErr", "pz","pzErr",
                        "d0Pull","omegaPull","phi0Pull","tanLambdaPull","z0Pull"]
    listofvariables += ["mcP","mcPT","mcPX","mcPY","mcPZ","mcDecayVertexX","mcDecayVertexY","mcDecayVertexZ","isSignal","mcErrors",
                        "mcPDG","charge"]  + vc.mc_vertex

    variables = vu.create_aliases(listofvariables,"{variable}","antiB0")
    variables += vu.create_aliases(listofvariables,"daughter(0,{variable})","Dstar")
    variables += vu.create_aliases(listofvariables,"daughter(0,daughter(0,{variable}))","D0")
    variables += vu.create_aliases(listofvariables,"daughter(0,daughter(1,{variable}))","Dstarpi")

    tag = "hitmatch_mdst_All_Dstar_04"+"_" + num_sigma+"_chi2"#"original_" + eventnum + "_" + num_sigma + "_" + gentype+"_hittune"#+"_nowinkel" ##
    filename = f"{basepath}/D/Data/modul_test/analysis/tuple_{tag}.root"
    ma.variablesToNtuple(particle, variables = variables, filename = filename , treename = "tree", path=main)
    #############################################################

if sys.argv[3]=="gun":
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
    
    tagend="_w_all"+"_hittune"+"_pdown4"#"_pdown4"#+"_pdown3"#+"_nowinkel"
    
    if part=="":
        stdPi("all",path=main)

        ma.fillParticleListFromMC("pi+:fun",cut="",skipNonPrimaryDaughters=True,path= main)

        particle= "pi+:all" #sys.argv[2] 
        variables = vu.create_aliases(listofvariables,"{variable}","all")
        particle2 = "pi+:fun"
        variables2 = vu.create_aliases(listofvariables,"{variable}","fun")
        p_up = ""
        for i in str(sys.argv[5]):
            if i !=".":
                p_up+=i
            else:
                p_up+="_"

        maintag = "original_2pi_"

        tag = maintag + eventnum +"_"+ num_sigma+ "_" + gentype+"_"+p_up+tagend#"_w_all"+"_chi2"+"_pdown4"#+"_pdown3"#+"_nowinkel"
        if maintag == "original_justpi_":
            filename = f"{basepath}/D/Data/modul_test/analysis/tuple_{tag}.root"

            ma.variablesToNtuple(particle, variables = variables, filename = filename , treename = "tree", path=main)

        if maintag=="original_2pi_":

            filename = f"{basepath}/D/Data/modul_test/analysis/tuple_{tag}.root"
            ma.variablesToNtuple(particle, variables = variables, filename = filename , treename = "all", path=main)
            ma.variablesToNtuple(particle2, variables = variables2, filename = filename , treename = "fun", path=main)
    elif part=="_electron":
        stdCharged.stdE("all",path=main)

        particle= "e-:all" #sys.argv[2] 
        
        #variables=listofvariables
        variables = vu.create_aliases(listofvariables,"{variable}","all")
        p_up = ""
        for i in str(sys.argv[5]):
            if i !=".":
                p_up+=i
            else:
                p_up+="_"

        maintag = "original_electron_"
    
        tag = maintag + eventnum +"_"+ num_sigma+ "_" + gentype+"_"+p_up+tagend
        filename = f"{basepath}/D/Data/modul_test/analysis/tuple_{tag}.root"
        ma.variablesToNtuple(particle, variables = variables, filename = filename , treename = "all", path=main)
    elif part=="_myon":
        stdCharged.stdMu("all",path=main)

        particle= "mu-:all" #sys.argv[2] 
        
        #variables=listofvariables
        variables = vu.create_aliases(listofvariables,"{variable}","all")
        p_up = ""
        for i in str(sys.argv[5]):
            if i !=".":
                p_up+=i
            else:
                p_up+="_"

        maintag = "original_myon_"
    
        tag = maintag + eventnum +"_"+ num_sigma+ "_" + gentype+"_"+p_up+tagend#"_w_all"+"_chi2"+"_pdown4"#+"_pdown3"#+"_nowinkel"
        filename = f"{basepath}/D/Data/modul_test/analysis/tuple_{tag}.root"
        ma.variablesToNtuple(particle, variables = variables, filename = filename , treename = "all", path=main)
    elif part=="_proton":
        stdCharged.stdPr("all",path=main)

        particle= "p+:all" #sys.argv[2] 
        
        #variables=listofvariables
        variables = vu.create_aliases(listofvariables,"{variable}","all")
        p_up = ""
        for i in str(sys.argv[5]):
            if i !=".":
                p_up+=i
            else:
                p_up+="_"

        maintag = "original_proton_"
    
        tag = maintag + eventnum +"_"+ num_sigma+ "_" + gentype+"_"+p_up+tagend#"_w_all"+"_chi2"+"_pdown4"#+"_pdown3"#+"_nowinkel"
        filename = f"{basepath}/D/Data/modul_test/analysis/tuple_{tag}.root"
        ma.variablesToNtuple(particle, variables = variables, filename = filename , treename = "all", path=main)
    #############################################################


#process the thing
process(main)
#print(statistics)


print("******************************* end Module ***************************************")
