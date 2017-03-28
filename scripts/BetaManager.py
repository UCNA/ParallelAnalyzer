#!/usr/bin/python

import os
import sys
import shutil
from optparse import OptionParser
from math import *
import MButils
#import CalibrationManager

##### Set up list of runs which are to be omitted from the beta data sets
##### ANYTIME RUNS ARE ADDED/SUBTRACTED FROM HERE, YOU MUST RERUN CREATEOCTETLISTS in this script
omittedBetaRuns = [19232]
omittedBetaRanges = [] 

betaRunTypes = ["A2","A5","A7","A10","B2","B5","B7","B10"]
bgRunTypes = ["A1","A4","A9","A12","B1","B4","B9","B12"]

for Range in omittedBetaRanges:
    for run in range(Range[0],Range[1]+1,1):
        omittedBetaRuns.append(run)


class BetaReplayManager:
    
    def __init__(self,AnalysisType="MB"):
        self.AnalyzerPath = "../"
        self.runListPath = self.AnalyzerPath + "run_lists/"
        self.AnalysisDataPath = os.getenv("PARALLEL_DATA_PATH")
        self.srcPositionsPath = os.getenv("SOURCE_POSITIONS")
        self.srcPeakPath = os.getenv("SOURCE_PEAKS")
        self.replayPass1 = os.getenv("REPLAY_PASS1")
        self.replayPass2 = os.getenv("REPLAY_PASS2")
        self.replayPass3 = os.getenv("REPLAY_PASS3")
        self.replayPass4 = os.getenv("REPLAY_PASS4")
        self.srcListPath = os.getenv("SOURCE_LIST")
        self.gainBismuthPath = os.getenv("GAIN_BISMUTH")
        self.gainCathodesPath = os.getenv("GAIN_CATHODES")
        self.mwpcCalPath = os.getenv("MWPC_CALIBRATION")
        self.nPEweightsPath = os.getenv("NPE_WEIGHTS")
        self.octetListPath = os.getenv("OCTET_LIST")
        self.triggerFuncPath = os.getenv("TRIGGER_FUNC")
        self.revCalSimPath = os.getenv("REVCALSIM")
        self.UKspecReplayPath = os.getenv("UK_SPEC_REPLAY")
        self.AnalysisResultsPath = os.getenv("ANALYSIS_RESULTS")
        self.SimAnalysisResultsPath = os.getenv("SIM_ANALYSIS_RESULTS")
        self.MPMAnalysisResultsPath = os.getenv("MPM_ANALYSIS_RESULTS")

    def makeAllDirectories(self):
        os.system("mkdir -p %s"%os.getenv("PEDESTALS"))
        os.system("mkdir -p %s"%os.getenv("BASIC_HISTOGRAMS"))
        os.system("mkdir -p %s"%self.srcPositionsPath)
        os.system("mkdir -p %s"%self.srcPeakPath)
        os.system("mkdir -p %s"%self.srcListPath)
        os.system("mkdir -p %s"%self.replayPass1)
        os.system("mkdir -p %s"%self.replayPass2)
        os.system("mkdir -p %s"%self.replayPass3)
        os.system("mkdir -p %s"%self.replayPass4)
        os.system("mkdir -p %s"%self.gainBismuthPath)
        os.system("mkdir -p %s"%self.nPEweightsPath)
        os.system("mkdir -p %s"%self.octetListPath)
        os.system("mkdir -p %s"%self.triggerFuncPath)
        os.system("mkdir -p %s"%self.revCalSimPath)
        os.system("mkdir -p %s"%self.UKspecReplayPath)
        os.system("mkdir -p %s"%self.AnalysisResultsPath)
        os.system("mkdir -p %s"%self.SimAnalysisResultsPath)
        os.system("mkdir -p %s"%self.MPMAnalysisResultsPath)
        os.system("mkdir -p %s"%self.mwpcCalPath)
        os.system("mkdir -p %s"%self.gainCathodesPath)
        os.system("mkdir -p %s/position_maps"%self.mwpcCalPath)
        os.system("mkdir -p %s/octets"%self.mwpcCalPath)
        os.system("mkdir -p %s/runs"%self.mwpcCalPath)
        os.system("mkdir -p %s/cathode_model"%self.mwpcCalPath)        
        os.system("mkdir -p %s/octets"%self.gainCathodesPath)
        os.system("mkdir -p %s/runs"%self.gainCathodesPath)

    def createOctetLists(self): #This creates lists of runs and type of run for each octet. There are 122 octets in the combined datasets
        octet=0
	for year in ["20112012","20122013"]:
            with open(self.octetListPath+"OctetList_%s.txt"%year) as fp:
                for line in fp:
                    if (line[:6]=="Octet:"):
			with open(self.octetListPath+"%s-%s/octet_list_%i.dat"%(year[:4],year[4:],octet), "wt") as out:
                            words=line.split()
                            for idx, w in enumerate(words):
                                if w[:1]=='A' or w[:1]=='B' and int(words[idx+2]) not in omittedBetaRuns:
                                    out.write(words[idx] + ' ' + words[idx+2] + '\n')
                        shutil.copy(self.octetListPath+"%s-%s/octet_list_%i.dat"%(year[:4],year[4:],octet),self.octetListPath+"All_Octets/octet_list_%i.dat"%(octet))
                        octet+=1
     
    def makeBasicHistograms(self, runORoctet):
        runs = []
        if runORoctet > 16000:
            print "Making basic histograms for run %i"%runORoctet
            os.system("cd ../basic_histograms/; ./basic_histograms.exe %i"%runORoctet)
            os.system("root -l -b -q '../basic_histograms/plot_basic_histograms.C(\"%i\")'"%runORoctet)
        else: 
            filename = "All_Octets/octet_list_%i.dat"%(runORoctet)
            infile = open(self.octetListPath+filename,'r')
        
            for line in infile:  
                words=line.split()
                runs.append(int(words[1]))
        
            for run in runs:
                print "Making basic histograms for run %i"%run
                os.system("cd ../basic_histograms/; ./basic_histograms.exe %i"%run)
                os.system("root -l -b -q '../basic_histograms/plot_basic_histograms.C(\"%i\")'"%run)
        print "DONE"
    

    def findTriggerFunctions(self, runORoctet):
        runs = []

        if runORoctet > 16000:
            print "Running trigger functions for run %i"%runORoctet
            if runORoctet>20000:
                os.system("cd ../trigger_functions/; ./findADCthreshold_singleRun.exe %i"%runORoctet)

        else: 
            filename = "All_Octets/octet_list_%i.dat"%(runORoctet)
            infile = open(self.octetListPath+filename,'r')
        
            for line in infile:  
                words=line.split()
                if words[0] in betaRunTypes: # Avoids background and depol runs
                    runs.append(int(words[1]))
        
            for run in runs:
                print "Running Trigger Functions for run %i"%run
                os.system("cd ../trigger_functions/; ./findADCthreshold_singleRun.exe %i"%run)

        print "DONE"


    def findPedestals(self, runORoctet):
        runs = []
        if runORoctet > 16000:
            print "Running pedestals for run %i"%runORoctet
            os.system("cd ../pedestals/; ./pedestals.exe %i"%runORoctet)
        else: 
            filename = "All_Octets/octet_list_%i.dat"%(runORoctet)
            infile = open(self.octetListPath+filename,'r')
        
            for line in infile:  
                words=line.split()
                if words[0] in betaRunTypes or words[0] in bgRunTypes: # Avoids depol runs
                    runs.append(int(words[1]))
        
            for run in runs:
                print "Running pedestals for run %i"%run
                os.system("cd ../pedestals/; ./pedestals.exe %i"%run)
                #os.system("cd ../pedestals/; ./pedestal_widths.exe %i"%run)
        print "DONE"

    def findBeamDrops(self, runORoctet):
        runs = []
        if runORoctet > 16000:
            print "Running beamDropCuts for run %i"%runORoctet
            os.system("cd ../cuts/; ./beamDropCuts.exe %i"%runORoctet)
        else: 
            filename = "All_Octets/octet_list_%i.dat"%(runORoctet)
            infile = open(self.octetListPath+filename,'r')
        
            for line in infile:  
                words=line.split()
                if words[0] in betaRunTypes or words[0] in bgRunTypes: # Avoids depol runs
                    runs.append(int(words[1]))
        
            for run in runs:
                print "Running pedestals for run %i"%run
                os.system("cd ../cuts/; ./beamDropCuts.exe %i"%run)
                #os.system("cd ../pedestals/; ./pedestal_widths.exe %i"%run)
        print "DONE"


    def runGainBismuth(self,runORoctet):
        runs = []
        if runORoctet > 16000:
            print "Running gain_bismuth for run %i"%runORoctet
            os.system("cd ../gain_bismuth/; ./gain_bismuth.exe %i"%runORoctet)
            os.system("root -l -b -q '../gain_bismuth/plot_gain_bismuth.C(\"%i\")'"%runORoctet)
        else: 
            filename = "All_Octets/octet_list_%i.dat"%(runORoctet)
            infile = open(self.octetListPath+filename,'r')
        
            for line in infile:      
                words=line.split()
                if words[0] in betaRunTypes or words[0] in bgRunTypes: # Avoids depol runs
                    runs.append(int(words[1]))
                    
        
            for run in runs:
                print "Running gain_bismuth for run %i"%run
                os.system("cd ../gain_bismuth/; ./gain_bismuth.exe %i"%run)
                os.system("root -l -b -q '../gain_bismuth/plot_gain_bismuth.C(\"%i\")'"%run)
        print "DONE"


    def runGainCathodes(self,octet):

        print "Running gain_cathodes for octet %i"%octet
        os.system("cd ../mwpc_cal/; ./cathode_gain.exe %i"%octet)

        runs = []
         
        filename = "All_Octets/octet_list_%i.dat"%(octet)
        infile = open(self.octetListPath+filename,'r')
        
        for line in infile:      
            words=line.split()
            if words[0] in betaRunTypes or words[0] in bgRunTypes: # Avoids depol runs
                runs.append(int(words[1]))
                    
        
        for run in runs:
            shutil.copy(self.gainCathodesPath+"/octets/gain_cathodes_octet_%i.dat"%octet,
                        self.gainCathodesPath+"/runs/gain_cathodes_%i.dat"%run)
                           
        print "DONE"


    def runCathodeModel(self,octet):

        print "Running cathode_model for octet %i"%octet

        runs = []
         
        filename = "All_Octets/octet_list_%i.dat"%(octet)
        infile = open(self.octetListPath+filename,'r')
        
        for line in infile:      
            words=line.split()
            if words[0] in betaRunTypes: # Avoids depol runs and bg runs since only used for sims
                runs.append(int(words[1]))
                    
        
        for run in runs:
            os.system("cd ../mwpc_cal/; ./cathode_model.exe %i"%run)
    
        print "DONE"




    def runWirechamberCal(self,runORoctet):

        runs = []
        bgruns = []
        
        if runORoctet > 16000:
            print "Running MWPC cal for run %i"%runORoctet
            os.system("cd ../mwpc_cal/; ./MWPC_calibration.exe %i"%runORoctet)
        else: 
            filename = "All_Octets/octet_list_%i.dat"%(runORoctet)
            infile = open(self.octetListPath+filename,'r')
        
            for line in infile:      
                words=line.split()
                if words[0] in betaRunTypes: # get FG runs
                    runs.append( int(words[1]) )
                if words[0] in bgRunTypes: # get BG runs
                    bgruns.append( int(words[1]) )

            numBG = len(bgruns)

            for run in runs:
                print "Running MWPC cal for run %i"%run
                os.system("cd ../mwpc_cal/; ./MWPC_calibration.exe %i"%run)
                
            for i in range(0,numBG,1): #as long as each FG has a BG, they will line up.
                    shutil.copy(self.mwpcCalPath+"/runs/MWPC_cal_%i.dat"%runs[i],
                                self.mwpcCalPath+"/runs/MWPC_cal_%i.dat"%bgruns[i])
                           
        print "DONE"
    

    def runReplayPass1(self,runORoctet):
        runs = []
        if runORoctet > 16000:
            print "Running replay_pass1 for run %i"%runORoctet
            os.system("cd ../replay_pass1/; ./replay_pass1.exe %i"%runORoctet)
        else: 
            filename = "All_Octets/octet_list_%i.dat"%(runORoctet)
            infile = open(self.octetListPath+filename,'r')
        
            for line in infile:  
                words=line.split()
                if words[0] in betaRunTypes or words[0] in bgRunTypes: # Avoids depol runs
                #if words[0] in bgRunTypes: # Avoids depol runs
                    runs.append(int(words[1]))

        
            for run in runs:
                print "Running replay_pass1 for run %i"%run
                os.system("cd ../replay_pass1/; ./replay_pass1.exe %i"%run)
        print "DONE"

        
    def runReplayPass2(self,runORoctet):
        runs = []
        if runORoctet > 16000:
            print "Running replay_pass2 for run %i"%runORoctet
            os.system("cd ../replay_pass2/; ./replay_pass2.exe %i"%runORoctet)
        else: 
            filename = "All_Octets/octet_list_%i.dat"%(runORoctet)
            infile = open(self.octetListPath+filename,'r')
        
            for line in infile:     
                words=line.split()
                if words[0] in betaRunTypes or words[0] in bgRunTypes: # Avoids depol runs
                #if words[0] in bgRunTypes: # Avoids depol runs
                    runs.append(int(words[1]))

        
            for run in runs:
                print "Running replay_pass2 for run %i"%run
                os.system("cd ../replay_pass2/; ./replay_pass2.exe %i"%run)
        print "DONE"

        
    def runReplayPass3(self,runORoctet):
        runs = []
        if runORoctet > 16000:
            print "Running replay_pass3 for run %i"%runORoctet
            os.system("cd ../replay_pass3/; ./replay_pass3.exe %i"%runORoctet)
        else: 
            filename = "All_Octets/octet_list_%i.dat"%(runORoctet)
            infile = open(self.octetListPath+filename,'r')
        
            for line in infile:      
                words=line.split()
                if words[0] in betaRunTypes or words[0] in bgRunTypes: # Avoids depol runs
                #if words[0] in bgRunTypes: # Avoids depol runs
                    runs.append(int(words[1]))

        
            for run in runs:
                print "Running replay_pass3 for run %i"%run
                os.system("cd ../replay_pass3/; ./replay_pass3.exe %i"%run)
        print "DONE"

    
    def removeDepolRunFiles(self,runORoctet):
        runs = []
        
        
        filename = "All_Octets/octet_list_%i.dat"%(runORoctet)
        infile = open(self.octetListPath+filename,'r')
        
        for line in infile:      
            words=line.split()
            if words[0] not in betaRunTypes and words[0] not in bgRunTypes: # Avoids depol runs
                os.system("rm %s/*%s*"%(os.getenv("PEDESTALS"),words[1]))
                os.system("rm %s/*%s*"%(self.gainBismuthPath,words[1]))
                os.system("rm %s/*%s*"%(self.replayPass1,words[1]))
                os.system("rm %s/*%s*"%(self.replayPass2,words[1]))
                os.system("rm %s/*%s*"%(self.replayPass3,words[1]))

            

    
    def runReverseCalibration(self, runORoctet, veryHighStatistics=False):
        runs = []

        if runORoctet > 16000:
            print "Running reverse calibration for run %i"%runORoctet
            os.system("./../simulation_comparison/revCalSim.exe %i %s"%(runORoctet, "Beta"))
            print "Finished reverse calibration for run %i"%runORoctet

        else: 
            filename = "All_Octets/octet_list_%i.dat"%(runORoctet)
            infile = open(self.octetListPath+filename,'r')
        
            for line in infile:      
                words=line.split()
                if words[0] in betaRunTypes:
                    runs.append(int(words[1]))
            
            
            for run in runs:
                print "Running reverse calibration for run %i"%run
                if veryHighStatistics:
                    os.system("./../simulation_comparison/revCalSim.exe %i %s %i"%(run, "Beta", runORoctet))
                else:
                    os.system("./../simulation_comparison/revCalSim.exe %i %s"%(run, "Beta"))

            print "Finished reverse calibration for octet %s"%runORoctet


    def runRootfileTranslator(self,runORoctet):
        runs = []
        if runORoctet > 16000:
            print "Running rootfile translator for run %i"%runORoctet
            os.system("cd ../rootfile_translator/; ./rootfile_translator %i"%runORoctet)
        else: 
            filename = "All_Octets/octet_list_%i.dat"%(runORoctet)
            infile = open(self.octetListPath+filename,'r')
        
            for line in infile:      
                words=line.split()
                if words[0] in betaRunTypes or words[0] in bgRunTypes:
                    runs.append(int(words[1]))
        
            for run in runs:
                print "Running rootfile translator for run %i"%run
                os.system("cd ../rootfile_translator/; ./rootfile_translator %i"%run)
        print "DONE"
        
    



class BetaAsymmetryManager:
    def __init__(self):
        self.AnalyzerPath = "../"
        self.runListPath = self.AnalyzerPath + "run_lists/"
        self.AnalysisDataPath = os.getenv("PARALLEL_DATA_PATH")
        self.srcPositionsPath = os.getenv("SOURCE_POSITIONS")
        self.srcPeakPath = os.getenv("SOURCE_PEAKS")
        self.replayPass1 = os.getenv("REPLAY_PASS1")
        self.replayPass2 = os.getenv("REPLAY_PASS2")
        self.replayPass3 = os.getenv("REPLAY_PASS3")
        self.replayPass4 = os.getenv("REPLAY_PASS4")
        self.srcListPath = os.getenv("SOURCE_LIST")
        self.gainBismuthPath = os.getenv("GAIN_BISMUTH")
        self.nPEweightsPath = os.getenv("NPE_WEIGHTS")
        self.octetListPath = os.getenv("OCTET_LIST")
        self.triggerFuncPath = os.getenv("TRIGGER_FUNC")
        self.revCalSimPath = os.getenv("REVCALSIM")
        self.UKspecReplayPath = os.getenv("UK_SPEC_REPLAY")
        self.AnalysisResultsPath = os.getenv("ANALYSIS_RESULTS")
        self.SimAnalysisResultsPath = os.getenv("SIM_ANALYSIS_RESULTS")
        self.MPMAnalysisResultsPath = os.getenv("MPM_ANALYSIS_RESULTS")
                
    def makeOctetAnalysisDirectories(self):
        for octet in range(0,122,1):
            os.system("mkdir -p %s/Octet_%i/OctetAsymmetry"%(self.AnalysisResultsPath,octet))
            os.system("mkdir -p %s/Octet_%i/QuartetAsymmetry"%(self.AnalysisResultsPath,octet))
            os.system("mkdir -p %s/Octet_%i/PairAsymmetry"%(self.AnalysisResultsPath,octet))

            os.system("mkdir -p %s/Octet_%i/OctetAsymmetry/Systematics"%(self.AnalysisResultsPath,octet))
            os.system("mkdir -p %s/Octet_%i/QuartetAsymmetry/Systematics"%(self.AnalysisResultsPath,octet))
            os.system("mkdir -p %s/Octet_%i/PairAsymmetry/Systematics"%(self.AnalysisResultsPath,octet))

            os.system("mkdir -p %s/Octet_%i/OctetAsymmetry"%(self.SimAnalysisResultsPath,octet))
            os.system("mkdir -p %s/Octet_%i/QuartetAsymmetry"%(self.SimAnalysisResultsPath,octet))
            os.system("mkdir -p %s/Octet_%i/PairAsymmetry"%(self.SimAnalysisResultsPath,octet))
            
            os.system("mkdir -p %s/Octet_%i/OctetAsymmetry"%(self.MPMAnalysisResultsPath,octet))
            os.system("mkdir -p %s/Octet_%i/QuartetAsymmetry"%(self.MPMAnalysisResultsPath,octet))
            os.system("mkdir -p %s/Octet_%i/PairAsymmetry"%(self.MPMAnalysisResultsPath,octet))
            
        
        os.system("mkdir -p %s/Asymmetries"%(self.AnalysisResultsPath))
        os.system("mkdir -p %s/Asymmetries"%(self.SimAnalysisResultsPath))
        os.system("mkdir -p %s/Asymmetries"%(self.MPMAnalysisResultsPath))
        print "Made all octet analysis directories"



if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option("--createOctetLists",dest="createOctetLists",action="store_true",default=False,
                      help="Creates Octet lists for individual octets stored in {dataPath}/octet_lists_MB")
    parser.add_option("--findPedestals",dest="findPedestals",action="store_true",default=False,
                      help="Finds the pedestals!")  
    parser.add_option("--replayPass1",dest="replayPass1",action="store_true",default=False,
                      help="Run replayPass1 on specified runs or octets (Applies cuts, reconstructs event positions, determines PID,Type,Side,etc.)")
    parser.add_option("--replayPass2",dest="replayPass2",action="store_true",default=False,
                      help="Run replayPass2 on specified runs or octets (Applies Bi pulser gain corrections")
    parser.add_option("--replayPass3",dest="replayPass3",action="store_true",default=False,
                      help="Run replayPass3 on specified runs or octets (Applies Xe position maps to scintillator response)")
    parser.add_option("--replayPass4",dest="replayPass4",action="store_true",default=False,
                      help="Run replayPass4 on specified runs or octets (Applies Energy Calibration)")
    parser.add_option("--runGainBismuth",dest="runGainBismuth",action="store_true",default=False,
                      help="Calculate and plot the Bi pulser spectra and track the gain relative to the reference run")
    parser.add_option("--makeDirectories",dest="makeDirectories",action="store_true",default=False,
                      help="Makes all the analysis directories.")
                      

    options, args = parser.parse_args()


    if options.makeDirectories:
        beta = BetaReplayManager()
        beta.makeAllDirectories()
        asymm = BetaAsymmetryManager()
        asymm.makeOctetAnalysisDirectories()

 
    if options.createOctetLists:
        beta = BetaReplayManager()
        beta.createOctetLists()
    
    if options.findPedestals:
        beta = BetaReplayManager()
        #for octet in range(0,60,1):
        beta.findPedestals(5)

    if 0:
        beta = BetaReplayManager()
        for octet in range(0,122,1):
            beta.makeBasicHistograms(octet)


    if 1: 
        octet_range =[109,121]#[20,28]#[45,50]#[38,40]#[0,59];
        beta = BetaReplayManager()
        for octet in range(octet_range[0],octet_range[1]+1,1):
            #beta.findPedestals(octet)
            #beta.runReplayPass1(octet)
            #beta.findBeamDrops(octet)
            #beta.runGainBismuth(octet)
            #beta.findTriggerFunctions(octet)
            #beta.runReplayPass2(octet)
            beta.runReplayPass3(octet)
            #beta.runRootfileTranslator(octet)
            #beta.removeDepolRunFiles(octet)
            beta.runReverseCalibration(octet)


    # Wirechamber stuff
    if 0:
        octet_range =[109,121]#[20,28]#[45,50]#[38,40]#[0,59];
        beta = BetaReplayManager()
        for octet in range(octet_range[0],octet_range[1]+1,1):
            
            #beta.runGainCathodes(octet)
            beta.runCathodeModel(octet)
            #beta.runWirechamberCal(octet)


    #Running reverse calibrations
    if 0:
        octet_range = [67,79];
        beta = BetaReplayManager()
        for octet in range(octet_range[0],octet_range[1]+1,1):
            #beta.findTriggerFunctions(octet)
            beta.runReverseCalibration(octet)
    
