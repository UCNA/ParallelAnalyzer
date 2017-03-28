#!/usr/bin/python

########################################################################
# This code is meant to control the replaying and energy calibration of 
# the source data for a given data set. It is currently set up for the 
# 2011/2012 data set and needs to be adapted for 2012/2013 unless we 
# enumerate these source runs starting with the next number in line after
# 2011/2012...
########################################################################

import os
import sys
from optparse import OptionParser
from math import *
import MButils

##### Set up list of runs which are to be omitted from the Energy Calibration
omittedRuns = [17588,17950,17953,18749,21298,21605]
#17588, 17950, 17953 are seemingly empty
#18749 - just an outlier on its west side... messes up fits in WPMT3 & 4
#19232 has abnormally high Bi peak energies on West side... no explanation
#21298 has very few events in Bi Pulser
#21605 has very few counts



omittedRanges = [(17923,18055),(20515,21086)] 
#These runs are from Run period 4 and include very long runs and runs with no Sn or Bi
# (20515,21086) had bad West side 2 fold triggers, This is src periods 13,14 and octets 60-67
#     and also the garbage from src period 16 and the bad octet right after that.

for Range in omittedRanges:
    for run in range(Range[0],Range[1]+1,1):
        omittedRuns.append(run)

##### Set up which runs to omit which PMTs.
#### This will be done by having a file with all source runs, where the PMT value
#### is set to 0 or 1 to represent false (don't use) and true (do use)

EPMT1 = [] #These hold individual runs where PMT was flaky or Bi pulser was not working. 
EPMT2 = []
EPMT3 = []
EPMT4 = []
WPMT1 = []
WPMT2 = []
WPMT3 = []
WPMT4 = []

EPMT1_runRanges = [] #These hold chunks of runs where PMT is dead or Bi pulser is not working.
EPMT2_runRanges = []
EPMT3_runRanges = []
EPMT4_runRanges = [(20000,23173)] #(17233,18055) not sure why these used to be removed... 2012-2013 have weird Bi gain and odd lin curves
WPMT1_runRanges = [(17359,18055)]
WPMT2_runRanges = [(16983,17297)] #PMTW2 dead for (16983,17297)
WPMT3_runRanges = []
WPMT4_runRanges = [(18712,19999),(20000,24000)]
#(18712,19999) WPMT4 becomes unreliable and there is no way to avoid affecting all of these runs due to bad Xe maps & Calibrations for this PMT

for Range in EPMT1_runRanges:
    for run in range(Range[0],Range[1]+1,1):
        EPMT1.append(run)
for Range in EPMT2_runRanges:
    for run in range(Range[0],Range[1]+1,1):
        EPMT2.append(run)
for Range in EPMT3_runRanges:
    for run in range(Range[0],Range[1]+1,1):
        EPMT3.append(run)
for Range in EPMT4_runRanges:
    for run in range(Range[0],Range[1]+1,1):
        EPMT4.append(run)
for Range in WPMT1_runRanges:
    for run in range(Range[0],Range[1]+1,1):
        WPMT1.append(run)
for Range in WPMT2_runRanges:
    for run in range(Range[0],Range[1]+1,1):
        WPMT2.append(run)
for Range in WPMT3_runRanges:
    for run in range(Range[0],Range[1]+1,1):
        WPMT3.append(run)
for Range in WPMT4_runRanges:
    for run in range(Range[0],Range[1]+1,1):
        WPMT4.append(run)

for run in omittedRuns:
    EPMT1.append(run)
    EPMT2.append(run)
    EPMT3.append(run)
    EPMT4.append(run)
    WPMT1.append(run)
    WPMT2.append(run)
    WPMT3.append(run)
    WPMT4.append(run)

class CalReplayManager:
    
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
        self.nPEweightsPath = os.getenv("NPE_WEIGHTS")
        self.octetListPath = os.getenv("OCTET_LIST")
        self.triggerFuncPath = os.getenv("TRIGGER_FUNC")
        self.revCalSimPath = os.getenv("REVCALSIM")
        self.UKspecReplayPath = os.getenv("UK_SPEC_REPLAY")
        self.AnalysisResultsPath = os.getenv("ANALYSIS_RESULTS")
        self.runInfoPath = os.getenv("RUN_INFO_FILES")

    def makeAllDirectories(self):
        #os.system("mkdir -p %s"%self.AnalysisDataPath)
        os.system("mkdir -p %s"%self.srcPositionsPath)
        os.system("mkdir -p %s"%os.getenv("BASIC_HISTOGRAMS"))
        os.system("mkdir -p %s"%os.getenv("PEDESTALS"))
        os.system("mkdir -p %s"%os.getenv("CUTS"))
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
        os.system("mkdir -p %s/sources/"%self.revCalSimPath)
        os.system("mkdir -p %s/beta/"%self.revCalSimPath)
        os.system("mkdir -p %s/beta_highStatistics/"%self.revCalSimPath)
        os.system("mkdir -p %s/source_peaks/"%self.revCalSimPath)
        os.system("mkdir -p %s"%self.UKspecReplayPath)
        os.system("mkdir -p %s"%self.AnalysisResultsPath)
        os.system("mkdir -p %s"%self.runInfoPath)

    def runReplayPass1(self,srcRunPeriod=1, sourceORxenon="source"):
        print "Running replay_pass1 for %s run period %i"%(sourceORxenon,srcRunPeriod)
        filename=None
        if sourceORxenon=="source":
            filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        elif sourceORxenon=="xenon":
            filename = "Xenon_Calibration_Run_Period_%i.dat"%srcRunPeriod
        else:
            print "Not a valid source type!! Options: source or xenon"
            exit();
        infile = open(self.runListPath+filename,'r')
        runs = []
        for line in infile:      
            if int(line) not in omittedRuns:
                runs.append(int(line))
        
        for run in runs:
            os.system("cd ../replay_pass1/; ./replay_pass1.exe %i"%run)
        print "DONE"
        
    def runReplayPass2(self,srcRunPeriod=1, sourceORxenon="source"):
        print "Running replay_pass2 for %s run period %i"%(sourceORxenon,srcRunPeriod)
        filename=None
        if sourceORxenon=="source":
            filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        elif sourceORxenon=="xenon":
            filename = "Xenon_Calibration_Run_Period_%i.dat"%srcRunPeriod
        else:
            print "Not a valid source type!! Options: source or xenon"
            exit();
        infile = open(self.runListPath+filename,'r')
        runs = []
        for line in infile:      
            if int(line) not in omittedRuns:
                runs.append(int(line))
        
        for run in runs:
            os.system("cd ../replay_pass2/; ./replay_pass2.exe %i false"%run)
        print "DONE"
        
    def runReplayPass3(self,srcRunPeriod=1, sourceORxenon="source"):
        print "Running replay_pass3 for %s run period %i"%(sourceORxenon,srcRunPeriod)
        filename=None
        if sourceORxenon=="source":
            filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        elif sourceORxenon=="xenon":
            filename = "Xenon_Calibration_Run_Period_%i.dat"%srcRunPeriod
        else:
            print "Not a valid source type!! Options: source or xenon"
            exit();
        infile = open(self.runListPath+filename,'r')
        runs = []
        for line in infile:      
            if int(line) not in omittedRuns:
                runs.append(int(line))
        
        for run in runs:
            os.system("cd ../replay_pass3/; ./replay_pass3.exe %i"%run)
        print "DONE"

    def runReplayPass4(self,srcRunPeriod=1, sourceORxenon="source"):
        print "Running replay_pass4 for %s run period %i"%(sourceORxenon,srcRunPeriod)
        os.system("mkdir -p %s"%(os.getenv("REPLAY_PASS4")))
        filename=None
        if sourceORxenon=="source":
            filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        elif sourceORxenon=="xenon":
            filename = "Xenon_Calibration_Run_Period_%i.dat"%srcRunPeriod
        else:
            print "Not a valid source type!! Options: source or xenon"
            exit();
        infile = open(self.runListPath+filename,'r')
        runs = []
        for line in infile:      
            if int(line) not in omittedRuns:
                runs.append(int(line))
        
        for run in runs:
            os.system("cd ../replay_pass4/; ./replay_pass4.exe %i"%run)
        print "DONE"
        
    def runGainBismuth(self,srcRunPeriod=1, sourceORxenon="source"):
        print "Running gain_bismuth for %s run period %i"%(sourceORxenon,srcRunPeriod)
        filename=None
        if sourceORxenon=="source":
            filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        elif sourceORxenon=="xenon":
            filename = "Xenon_Calibration_Run_Period_%i.dat"%srcRunPeriod
        else:
            print "Not a valid source type!! Options: source or xenon"
            exit();
        infile = open(self.runListPath+filename,'r')
        runs = []
        for line in infile:      
            if int(line) not in omittedRuns:
                runs.append(int(line))
        
        for run in runs:
            os.system("cd ../gain_bismuth/; ./gain_bismuth.exe %i"%run)
            os.system("root -l -b -q '../gain_bismuth/plot_gain_bismuth.C(\"%i\")'"%run)
        print "DONE"


    def makeBasicHistograms(self, srcRunPeriod=1, sourceORxenon="source"):
        print "Making Basic Histograms for %s run period %i"%(sourceORxenon,srcRunPeriod)
        filename=None
        if sourceORxenon=="source":
            filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        elif sourceORxenon=="xenon":
            filename = "Xenon_Calibration_Run_Period_%i.dat"%srcRunPeriod
        else:
            print "Not a valid source type!! Options: source or xenon"
            exit();
        infile = open(self.runListPath+filename,'r')
        runs = []
        for line in infile:      
            if int(line) not in omittedRuns:
                runs.append(int(line))
        
        for run in runs:
            print "Making basic histograms for run %i"%run
            os.system("cd ../basic_histograms/; ./basic_histograms.exe %i"%run)
            os.system("root -l -b -q '../basic_histograms/plot_basic_histograms.C(\"%i\")'"%run)
            
        print "DONE"


   
    def findTriggerFunctions(self, srcRunPeriod=1, sourceORxenon="source"):
        print "Running trigger functions for %s run period %i"%(sourceORxenon,srcRunPeriod)
        filename=None
        if sourceORxenon=="source":
            filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        elif sourceORxenon=="xenon":
            filename = "Xenon_Calibration_Run_Period_%i.dat"%srcRunPeriod
        else:
            print "Not a valid source type!! Options: source or xenon"
            exit();
        infile = open(self.runListPath+filename,'r')
        runs = []
        for line in infile:      
            if int(line) not in omittedRuns:
                runs.append(int(line))
        
        for run in runs:
            os.system("cd ../trigger_functions/; ./findADCthreshold_singleRun.exe %i"%run)
            
        print "DONE"


    def findPedestals(self, srcRunPeriod=1, sourceORxenon="source"):
        print "Running pedestals for %s run period %i"%(sourceORxenon,srcRunPeriod)
        filename=None
        if sourceORxenon=="source":
            filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        elif sourceORxenon=="xenon":
            filename = "Xenon_Calibration_Run_Period_%i.dat"%srcRunPeriod
        else:
            print "Not a valid source type!! Options: source or xenon"
            exit();
        infile = open(self.runListPath+filename,'r')
        runs = []
        for line in infile:      
            if int(line) not in omittedRuns:
                runs.append(int(line))
        
        for run in runs:
            os.system("cd ../pedestals/; ./pedestals.exe %i"%run)
            #os.system("cd ../pedestals/; ./pedestal_widths.exe %i"%run)
            
        print "DONE"


    def runReverseCalibration(self, srcRunPeriod=1, sourceORxenon="source"):
        print "Running reverse calibration for %s run period %i"%(sourceORxenon,srcRunPeriod)
        filename=None
        if sourceORxenon=="source":
            filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        elif sourceORxenon=="xenon":
            filename = "Xenon_Calibration_Run_Period_%i.dat"%srcRunPeriod
        else:
            print "Not a valid source type!! Options: source or xenon"
            exit();
        infile = open(self.runListPath+filename,'r')
        runs = []

        for line in infile:      
          if int(line) not in omittedRuns:
              runs.append(int(line))

        for run in runs:
            filename = self.srcListPath+"source_list_%i.dat"%run
            if not MButils.fileExistsAndNotEmpty(filename):
                continue
            srcFile = open(filename)
            srcFileInput = []
            for line in srcFile:
                srcFileInput.append(line)
            for src in range(1,int(srcFileInput[0])+1):
                source = srcFileInput[src][0:2]
                if source=="Ce":
                    source = source+"139"
                elif source=="Sn":
                    source = source + "113"
                elif source=="Bi":
                    source = source + "207"
                elif source == "In":
                    source = source + "114"
                elif source == "Cd":
                    source = source + "109"
                elif source == "Cs":
                    source = source+"137"

                if source in ["Ce139","Sn113","Bi207"]: #,"In114"
                    os.system("cd ../simulation_comparison/;./revCalSim.exe %i %s"%(run, source))
                    #print "./../simulation_comparison/revCalSim.exe %i %s"%(run, source)

        print "Finished reverse calibration for " + sourceORxenon + "run period %i"%srcRunPeriod


        
class CalibrationManager:
    
    def __init__(self,AnalysisType="MB"):
        self.AnalyzerPath = "../"
        self.runListPath = self.AnalyzerPath + "run_lists/"
        self.AnalysisDataPath = os.getenv("PARALLEL_DATA_PATH")
        self.srcPositionsPath = os.getenv("SOURCE_POSITIONS")
        self.srcPeakPath = os.getenv("SOURCE_PEAKS")
        self.replayPass3 = os.getenv("REPLAY_PASS3")
        self.srcListPath = os.getenv("SOURCE_LIST")
        self.nPEcountPath = os.getenv("NPE_WEIGHTS")
        self.revCalSimPath = os.getenv("REVCALSIM")

    def findCalibrationPeriod(self, runNumber):
        if runNumber <= 17297:
            return 1
        elif runNumber <= 17439:
            return 2
        elif runNumber <= 17734:
            return 3
        elif runNumber <= 17955:
            return 4
        elif runNumber <= 18386:
            return 5
        elif runNumber <= 18683:
            return 6
        elif runNumber <= 18994:
            return 7
        elif runNumber <= 19239:
            return 8
        elif runNumber <= 19544:
            return 9
        elif runNumber < 20000:
            return 11
        elif runNumber <= 20741:
            return 13
        elif runNumber <= 20837:
            return 14
        elif runNumber <= 21237:
            return 16
        elif runNumber <= 21605:
            return 17
        elif runNumber <= 21863:
            return 18
        elif runNumber <= 22118:
            return 19
        elif runNumber <= 22238:
            return 20
        elif runNumber <= 22630:
            return 22
        elif runNumber <= 23173:
            return 23
        else:
            print "Bad Run Number!!! No Calibration Period..."
            exit


    def runSourceCalReplayPeakFitter(self,srcRunPeriod=1, doOnlyCeSnBi=True):
        print "Running SrcCalReplay for run period %i"%(srcRunPeriod)
        filename=None
        
        filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        infile = open(self.runListPath+filename,'r')
        runs = []
        for line in infile:  
            if int(line) not in omittedRuns:
                runs.append(int(line))
            
        for run in runs:
            filename = self.srcListPath+"source_list_%i.dat"%run
            if not MButils.fileExistsAndNotEmpty(filename):
                continue
            srcFile = open(filename)
            sources = []
            for line in srcFile:
                sources.append(line)

            isCeSnBi = False
            for source in sources:
                if source[0:2]=="Ce":
                    isCeSnBi = True
                    break
                if source[0:2]=="Sn":
                    isCeSnBi = True
                    break
                if source[0:2]=="Bi":
                    isCeSnBi = True
                    break
                
            if doOnlyCeSnBi:
                if isCeSnBi:
                    os.system("cd ../source_peaks/; ./srcCalReplay.exe %i"%run)

            else: 
                os.system("cd ../source_peaks/; ./srcCalReplay.exe %i"%run)

        print "DONE"
                
        

    def calc_nPE_per_PMT(self, runAllRefRun=False, run=19359, writeNPEforAllRuns=False, year="2011-2012"):

        srcSn_nPE_Runs = [17238,17370,17521,17925,18361,18621,18749,19232,19359,19511,19857,19899,20520,20831,20905,21091,21315,21683,21918,22219,22298,22441,22771,22925] #These are runs for which the nPE per channel are calculated
        # May or may not match the reference runs

        srcRunPeriodRange = None
        runRange = None

        if year=="2011-2012":
            srcRunPeriodRange = [1,12]
            runRange = [16000,20000]

        elif year=="2012-2013":
            srcRunPeriodRange = [13,24]
            runRange = [20000,23173]

        if runAllRefRun:
            
            for srcRunPeriod in range(srcRunPeriodRange[0],srcRunPeriodRange[1]+1,1):
                r = srcSn_nPE_Runs[srcRunPeriod-1]
                os.system("cd ../calc_nPE/; ./calc_nPE.exe %i"%r)
                print "Ran calc_nPE.exe for run %i"%r
                #runs=[]
                #filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
                #infile = open(self.runListPath+filename,'r')
                #for line in infile:
                #    runs.append(int(line))
                #for rn in runs:
                #    os.system("cp %s/nPE_weights_%i.dat %s/nPE_weights_%i.dat"%(self.nPEcountPath,r,self.nPEcountPath,rn))
                #infile.close()

        elif writeNPEforAllRuns:
            for rn in range(runRange[0],runRange[1]+1,1):
                calPeriod = self.findCalibrationPeriod(rn)
                os.system("cp %s/nPE_weights_%i.dat %s/nPE_weights_%i.dat"%(self.nPEcountPath,srcSn_nPE_Runs[calPeriod-1],self.nPEcountPath,rn))
                os.system("cp %s/nPE_meanEtaVal_%i.dat %s/nPE_meanEtaVal_%i.dat"%(self.nPEcountPath,srcSn_nPE_Runs[calPeriod-1],self.nPEcountPath,rn))

        else:
            os.system("cd ../calc_nPE/; ./calc_nPE.exe %i"%run)
            print "Ran calc_nPE.exe for run %i"%run

        
        
        
        

    def fitSourcePositions(self,srcRunPeriod=1, overwrite=False):
        filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        infile = open(self.runListPath + filename,'r')
        runs = []
        for line in infile:
        #checking if positions have already been fit
            filepath = self.srcPositionsPath + "source_positions_%i.dat"%int(line)
            print filepath
            if not MButils.fileExistsAndNotEmpty(filepath) or overwrite:
                runs.append(int(line))

        for run in runs:
            os.system("root -b -q '../source_positions/fit_source_positions.C(\"%i\")'"%run)
            print "Running fit_source_positions.C on run %i"%run
        

    def fitSimSourcePeaks(self,srcRunPeriod=1, overwrite=True):
        filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        infile = open(self.runListPath+filename,'r')
        runs = []

        for line in infile:  
            if int(line) not in omittedRuns:
                #checking if peaks have already been fit
                filepath = self.srcPeakPath +"source_peaks_%i.dat"%int(line)
                print filepath
                if not MButils.fileExistsAndNotEmpty(filepath) or overwrite:
                    runs.append(int(line))
             
            

        for run in runs:
            
            os.system("cd ../source_peaks; ./sim_source_peak_fitter.exe %i"%run)
            #os.system("root -b -q '../source_peaks/plot_sim_source_peaks.C(\"%i\")'"%run)
            print "Ran sim_source_peak_fitter.exe on run %i"%run

            
    def fitSourcePeaks(self,srcRunPeriod=1, Simulation=False, overwrite=True):
        filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        infile = open(self.runListPath+filename,'r')
        runs = []
        for line in infile:       
            #checking if peaks have already been fit
            filepath = self.srcPeakPath +"source_peaks_%i.dat"%int(line)
            print filepath
            if not MButils.fileExistsAndNotEmpty(filepath) or overwrite:
                runs.append(int(line))

        for run in runs:
            if not Simulation:
                os.system("cd ../source_peaks; ./source_peaks.exe %i"%run)
                os.system("root -b -q '../source_peaks/plot_source_peaks.C(\"%i\")'"%run)
                print "Ran fit_source_peaks.C on run %i"%run
            else:
                os.system("cd ../source_peaks; ./sim_source_peaks.exe %i"%run)
                os.system("root -b -q '../source_peaks/plot_sim_source_peaks.C(\"%i\")'"%run)
                print "Ran fit_source_peaks.C on run %i"%run


    def fitSourcePeaksInEnergy(self,srcRunPeriod=1, Simulation=False, overwrite=True):
        filename = "Source_Calibration_Run_Period_%i.dat"%srcRunPeriod
        infile = open(self.runListPath+filename,'r')
        runs = []
        for line in infile:       
            #checking if peaks have already been fit
            filepath = self.srcPeakPath +"source_peaks_%i_Evis.dat"%int(line)
                
            print filepath
            if not MButils.fileExistsAndNotEmpty(filepath) or overwrite:
                runs.append(int(line))

        for run in runs:
            
            if Simulation:
                os.system("cd %s/source_peaks; ./sim_source_peaks_postCal.exe %i"%(os.getenv("ANALYSIS_CODE"),run))
                os.system("root -b -q '%s/source_peaks/plot_sim_source_peaks_postCal.C(\"%i\")'"%(os.getenv("ANALYSIS_CODE"),run))
            else:
                os.system("cd %s/source_peaks; ./source_peaks_postCal.exe %i"%(os.getenv("ANALYSIS_CODE"),run))
                os.system("root -b -q '%s/source_peaks/plot_source_peaks_postCal.C(\"%i\")'"%(os.getenv("ANALYSIS_CODE"),run))
                                
            print "Ran plot_source_peaks.C on run %i"%run


    def makePMTrunFile(self,CalibrationPeriod=1, master=False):
        if not master:
            outputFile = "%s/residuals/PMT_runQuality_SrcPeriod_%i.dat"%(os.getenv("ANALYSIS_CODE"),CalibrationPeriod)
            runList = []

            with open("%s/run_lists/Source_Calibration_Run_Period_%i.dat"%(os.getenv("ANALYSIS_CODE"),CalibrationPeriod)) as runlist:
                for run in runlist:
                    if os.path.isfile(self.srcListPath+"source_list_%i.dat"%int(run)) and int(run) not in omittedRuns:
                        srcList = open(self.srcListPath+"source_list_%i.dat"%int(run))
                        lines = []
                        for line in srcList:
                            lines.append(line)
                        if int(lines[0])>0:
                            runList.append(int(run))

            outfile = open(outputFile,'w')

            for run in runList:
                pmtList = [1,1,1,1,1,1,1,1]
                if run in EPMT1:
                    pmtList[0]=0
                if run in EPMT2:
                    pmtList[1]=0
                if run in EPMT3:
                    pmtList[2]=0
                if run in EPMT4:
                    pmtList[3]=0
                if run in WPMT1:
                    pmtList[4]=0
                if run in WPMT2:
                    pmtList[5]=0
                if run in WPMT3:
                    pmtList[6]=0
                if run in WPMT4:
                    pmtList[7]=0

                outfile.write("%i %i %i %i %i %i %i %i %i\n"%(run,pmtList[0],pmtList[1],pmtList[2],pmtList[3],
                                                          pmtList[4],pmtList[5],pmtList[6],pmtList[7]))

            outfile.close()
            print "Done writing PMT file for Source Period %i"%CalibrationPeriod

        #Update the master list of PMT quality
        if master:
            masterFile = open("%s/residuals/PMT_runQuality_master.dat"%(os.getenv("ANALYSIS_CODE")),'w')
            
            for run in range(16983,23174,1):
                pmtList = [1,1,1,1,1,1,1,1]
                if run in EPMT1:
                    pmtList[0]=0
                if run in EPMT2:
                    pmtList[1]=0
                if run in EPMT3:
                    pmtList[2]=0
                if run in EPMT4:
                    pmtList[3]=0
                if run in WPMT1:
                    pmtList[4]=0
                if run in WPMT2:
                    pmtList[5]=0
                if run in WPMT3:
                    pmtList[6]=0
                if run in WPMT4:
                    pmtList[7]=0

                masterFile.write("%i %i %i %i %i %i %i %i %i\n"%(run,pmtList[0],pmtList[1],pmtList[2],pmtList[3],
                                                          pmtList[4],pmtList[5],pmtList[6],pmtList[7]))
            masterFile.close()
            print "Updated master list of PMT run quality"



    def makeSourceCalibrationFile(self,CalibrationPeriod=1, Simulation=False, InEnergy=False):
        #This utilizes the omittedRuns and removes them from the calibration. Any time you make a change
        # to the runs which are to be omitted, you shoud rerun this!

        if not InEnergy:

            src_file_base = None
            peak_file = None
            peak_error_file = None
            fileEnding = None
            
            if Simulation:
                peak_file = "%s/residuals/SIM_source_runs_RunPeriod_%i.dat"%(os.getenv("ANALYSIS_CODE"),CalibrationPeriod)
                peak_error_file = "%s/residuals/SIM_source_runs_peakErrors_RunPeriod_%i.dat"%(os.getenv("ANALYSIS_CODE"),CalibrationPeriod)
                src_file_base = self.revCalSimPath + "/source_peaks/source_peaks_"
                fileEnding = "_etaEvis.dat"

            else:
                peak_file = "%s/residuals/source_runs_RunPeriod_%i.dat"%(os.getenv("ANALYSIS_CODE"),CalibrationPeriod)
                peak_error_file = "%s/residuals/source_runs_peakErrors_RunPeriod_%i.dat"%(os.getenv("ANALYSIS_CODE"),CalibrationPeriod)
                src_file_base = self.srcPeakPath + "source_peaks_"
                fileEnding = "_ADC.dat"
                
            runList = []

            with open("%s/run_lists/Source_Calibration_Run_Period_%i.dat"%(os.getenv("ANALYSIS_CODE"),CalibrationPeriod)) as runlist:
                for run in runlist:
                    if os.path.isfile(self.srcListPath+"source_list_%i.dat"%int(run)):
                        srcList = open(self.srcListPath+"source_list_%i.dat"%int(run))
                        lines = []
                        for line in srcList:
                            lines.append(line)
                        if int(lines[0])>0:
                            runList.append(int(run))

            #print runList

            outfile = open(peak_file,'w')
            outfileError = open(peak_error_file,'w')

            for run in runList:
                src_file = src_file_base+"%i"%run+fileEnding
                if MButils.fileExistsAndNotEmpty(src_file) and run not in omittedRuns:
                    infile = open(src_file,'r')
                    for line in infile:
                        outfile.write(line)
                    infile.close()

                src_file_err = src_file_base+"errors_%i"%run+fileEnding
                if MButils.fileExistsAndNotEmpty(src_file_err) and run not in omittedRuns:
                    infile = open(src_file_err,'r')
                    for line in infile:
                        outfileError.write(line)
                    infile.close()

            outfile.close()
            outfileError.close()

        else:

            EvisFile = None
            EvisWidthFile = None
            EreconFile = None
            src_file_base = None

            if Simulation: 
                src_file_base = self.revCalSimPath + "/source_peaks/"
                EvisFile = "%s/residuals/SIM_source_runs_Evis_RunPeriod_%i.dat"%(os.getenv("ANALYSIS_CODE"),CalibrationPeriod)
                EvisWidthFile = "%s/residuals/SIM_source_runs_EvisWidth_RunPeriod_%i.dat"%(os.getenv("ANALYSIS_CODE"),CalibrationPeriod)
                EreconFile = "%s/residuals/SIM_source_runs_Erecon_RunPeriod_%i.dat"%(os.getenv("ANALYSIS_CODE"),CalibrationPeriod)
            else:
                src_file_base = self.srcPeakPath
                EvisFile = "%s/residuals/source_runs_Evis_RunPeriod_%i.dat"%(os.getenv("ANALYSIS_CODE"),CalibrationPeriod)
                EvisWidthFile = "%s/residuals/source_runs_EvisWidth_RunPeriod_%i.dat"%(os.getenv("ANALYSIS_CODE"),CalibrationPeriod)
                EreconFile = "%s/residuals/source_runs_Erecon_RunPeriod_%i.dat"%(os.getenv("ANALYSIS_CODE"),CalibrationPeriod)
        

            runList = []

            with open("%s/run_lists/Source_Calibration_Run_Period_%i.dat"%(os.getenv("ANALYSIS_CODE"),CalibrationPeriod)) as runlist:
                for run in runlist:
                    if os.path.isfile(self.srcListPath+"source_list_%i.dat"%int(run)):
                        srcList = open(self.srcListPath+"source_list_%i.dat"%int(run))
                        lines = []
                        for line in srcList:
                            lines.append(line)
                        if int(lines[0])>0:
                            runList.append(int(run))

            #print runList

            EvisOutfile = open(EvisFile,'w')
            EvisWidthOutfile = open(EvisWidthFile,'w')
            EreconOutfile = open(EreconFile,'w')

            for run in runList:
                src_file = src_file_base+"/source_peaks_%i_Evis.dat"%run
                if MButils.fileExistsAndNotEmpty(src_file) and run not in omittedRuns:
                    infile = open(src_file,'r')
                    for line in infile:
                        EvisOutfile.write(line)
                    infile.close()

                src_file = src_file_base + "/source_widths_%i_Evis.dat"%run
                if MButils.fileExistsAndNotEmpty(src_file) and run not in omittedRuns:
                    infile = open(src_file,'r')
                    for line in infile:
                        EvisWidthOutfile.write(line)
                    infile.close()

                src_file = src_file_base + "/source_peaks_%i_EreconTot.dat"%run
                if MButils.fileExistsAndNotEmpty(src_file) and run not in omittedRuns:
                    infile = open(src_file,'r')
                    for line in infile:
                        EreconOutfile.write(line)
                    infile.close()

            EvisOutfile.close()
            EvisWidthOutfile.close()
            EreconOutfile.close()

        print "Made combined source peak file for Calibration Period %i"%CalibrationPeriod


    #### Calculate the linearity curves for each PMT from a source run period
    def LinearityCurves(self,CalibrationPeriod=1):
        filename = "%s/residuals/source_runs_RunPeriod_%i.dat"%(os.getenv("ANALYSIS_CODE"),CalibrationPeriod) #Files which hold fitted source peaks
        if os.path.isfile(filename):
            os.system("cd %s/linearity_curves/; root -b -q 'LinearityCurves.C (%i)'"%(os.getenv("ANALYSIS_CODE"),CalibrationPeriod))
            print "Calculated linearity curves for Calibration Period %i"%CalibrationPeriod
        else:
            print "No peak file to calculate linearity curves"
            sys.exit


            
    def calc_new_nPE_per_keV(self, CalibrationPeriod):
        print "Calculating new nPE per keV for Calibration Period %i"%CalibrationPeriod
        os.system("cd %s/source_peaks/; root -l -b -q 'width_fitter.C(%i)'"%(os.getenv("ANALYSIS_CODE"),CalibrationPeriod))
        print "DONE"



    #Calculates residuals for each PMT and as a whole for given run period
    def calculateResiduals(self, CalibrationPeriod=1):

        filename = "%s/residuals/source_runs_Erecon_RunPeriod_%i.dat"%(os.getenv("ANALYSIS_CODE"),CalibrationPeriod)

        if os.path.isfile(filename):
            
            os.system("root -b -q 'calc_residuals_postCal.C (%i)'"%CalibrationPeriod)
            
            print "Calculated residuals for Calibration Period %i"%CalibrationPeriod
        else:
            print "No peak file to calculate residuals"
            sys.exit



    #Combines residuals for the calPeriods in the list given, for the PMT chosen (0 is PMT as a whole), and for the side
    # which can be "East", "West", or "Both"
    def makeGlobalResiduals(self,CalPeriods=[1]):
        CalPeriods.sort()
        periodLow = CalPeriods[0]
        periodHigh = CalPeriods[len(CalPeriods)-1]

        outfile = open("%s/residuals/global_residuals_Erecon_runPeriods_%i-%i.dat"%(os.getenv("ANALYSIS_CODE"),periodLow,periodHigh),"w")

        for period in CalPeriods:
            filename = "%s/residuals/residuals_Erecon_runPeriod_%i.dat"%(os.getenv("ANALYSIS_CODE"),period)

            if os.path.isfile(filename):
                resid = open(filename)
                for line in resid:
                    outfile.write(line)

        outfile.close()

        print "Produced global residual file for weighted average of all PMTs, run periods %i-%i"%(periodLow,periodHigh)

        


    def plotErrorEnvelope(self, calPeriodLow=2, calPeriodHigh=10, PMT=0, InEnergy=False):
        ## This runs code which calculates the mean and RMS of the global residual file
        ## from Calibration run periods CalPeriodLow to CalPeriodHigh. It prints out the mean and RMS
        ## for each source. Later this will be input into code which actually plots the error envelope.
        
        filenameEast = None
        filenameWest = None

        if calPeriodLow!=calPeriodHigh and not PMT:
            filenameEast = "../residuals/residuals_global_East_periods_%i-%i.dat"%(calPeriodLow, calPeriodHigh)
            filenameWest = "../residuals/residuals_global_West_periods_%i-%i.dat"%(calPeriodLow, calPeriodHigh)
        elif calPeriodLow==calPeriodHigh and not PMT:
            filenameEast = "../residuals/residuals_East_runPeriod_%i.dat"%calPeriodLow
            filenameWest = "../residuals/residuals_West_runPeriod_%i.dat"%calPeriodLow
        elif calPeriodLow!=calPeriodHigh and PMT:
            filenameEast = "../residuals/residuals_global_East_periods_%i-%i_PMTE%i.dat"%(calPeriodLow, calPeriodHigh, PMT)
            filenameWest = "../residuals/residuals_global_West_periods_%i-%i_PMTW%i.dat"%(calPeriodLow, calPeriodHigh, PMT)
        elif calPeriodLow==calPeriodHigh and PMT:
            filenameEast = "../residuals/residuals_East_runPeriod_%i_PMTE%i.dat"%(calPeriodLow,PMT)
            filenameWest = "../residuals/residuals_West_runPeriod_%i_PMTW%i.dat"%(calPeriodLow,PMT)
        
        if MButils.fileExistsAndNotEmpty(filenameEast) and MButils.fileExistsAndNotEmpty(filenameWest):
            print "Making Error Envelope for Run Periods %i to %i"%(calPeriodLow,calPeriodHigh)
            os.system("root -l -b -q 'MB_errorEnvelope.C (%i,%i,%i)'"%(calPeriodLow,calPeriodHigh,PMT))



if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option("--fitSrcPeaks",dest="fitSrcPeaks",action="store_true",default=False,
                      help="Fit Source peaks for all runs unless otherwise stated")
    parser.add_option("--fitSrcPeaksInEnergy",dest="fitSrcPeaksInEnergy",action="store_true",default=False,
                      help="Fit Source peaks post replay_pass4 for all runs unless otherwise stated")
    parser.add_option("--fitSrcPositions",dest="fitSrcPositions",action="store_true",default=False,
                      help="Fit for Source Positions for all runs. Only to be done if new replay_pass3 has been done on sources!")
    parser.add_option("--makeAllCalFiles",dest="makeAllCalFiles",action="store_true",default=False,
                      help="Recombines all source peaks for each Calibration period, then calculates residuals, and finally makes global residual file")
    parser.add_option("--makePeakFiles",dest="makePeakFiles",action="store_true",default=False,
                      help="Make combined source peak files for individual source calibration periods.")
    parser.add_option("--calcResiduals",dest="calcResiduals",action="store_true",default=False,
                      help="Calculate the residuals for the given run periods")
    parser.add_option("--makeGlobalResiduals",dest="makeGlobalResiduals",action="store_true",default=False,
                      help="Combine all residuals into one file for drawing error envelope")
    parser.add_option("--ErrorEnvelope",dest="ErrorEnvelope",action="store_true",default=False,
                      help="Make error envelope and save mean and sigma to file in ../error_envelope.")
    parser.add_option("--makePMTrunFile",dest="makePMTrunFile",action="store_true",default=False,
                      help="Make file with booleans for whether to use each PMT for each run.")
    parser.add_option("--LinearityCurves",dest="LinearityCurves",action="store_true",default=False,
                      help="Run script which calculates linearity curves for PMTs.")
    parser.add_option("--makeDirectories",dest="makeDirectories",action="store_true",default=False,
                      help="Makes all the analysis directories.")
                      

    options, args = parser.parse_args()


    if options.makeDirectories:
        cal = CalReplayManager()
        cal.makeAllDirectories()


    ### This will fit all the source peaks for any runs in the runPeriods list below
    if options.fitSrcPeaks:
        
        runPeriods = [1,2,3,4,5,6,7,8,9,10,11,12]
        cal = CalibrationManager()
        for period in runPeriods:
            cal.fitSourcePeaks(period)

    if options.fitSrcPeaksInEnergy:
        
        runPeriods = [1,2,3,4,5,6,7,8,9,10,11,12]
        cal = CalibrationManager()
        for period in runPeriods:
            cal.fitSourcePeaksInEnergy(period)


    ### This will re-fit all the source positions. Note that by default, even with this option as true, 
    ### the overwrite option is False. These source positions have been checked by hand, and the only
    ### reason to redo them is if the position maps change.
    if options.fitSrcPositions:
        
        runPeriods = [1,2,3,4,5,6,7,8,9,10,11,12]
        cal = CalibrationManager()
        for period in runPeriods:
            cal.fitSourcePositions(period,False)
        
    ### Makes a file with each run followed by a boolean (0,1) for whether each PMT should be used or not
    if options.makePMTrunFile:
        runPeriods = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]#
        cal=CalibrationManager()
        for period in runPeriods:
            cal.makePMTrunFile(period)
            
        cal.makePMTrunFile(1,True) #Updates the master list of PMT quality over all runs


    ### If you have made changes to the runs which are to be ignored at the top of this script, you should run this 
    #def makeSourceCalibrationFile(self,CalibrationPeriod=1, PeaksInEnergy=False, PMTbyPMT=False, Simulation=False):
    if options.makeAllCalFiles:
    
        runPeriods = [1,2,3,4,5,6,7,8,9,10,11,12]
        cal = CalibrationManager()
        for period in runPeriods:
            #cal.makeSourceCalibrationFile(period)
            #cal.makePMTrunFile(period, master=False)
            cal.calculateResiduals(period)
            
        #cal.makeGlobalResiduals(runPeriods,PMT=1,Side="Both")


    ### useful if you are going to look at the linearity curves and residuals by eye separately
    if options.makePeakFiles:
        runPeriods = [13,14,16,17,18,19,20,21,22,23]#[1,2,3,4,5,6,7,8,9,10,11,12]
        cal=CalibrationManager()
        for period in runPeriods:
            cal.makeSourceCalibrationFile(period, False)


    ### Does the linearity curves for all the source calibration periods
    if options.LinearityCurves:
        runPeriods =[1,2,3,4,5,6,7,8,9,10,11,12]
        cal = CalibrationManager()
        for runPeriod in runPeriods:
            cal.LinearityCurves(runPeriod)


    ### Saves the results of calculating the mean and RMS of all the global residuals for the 
    ### given combination of calibration periods and for a certain PMT (PMT=0 is for all 4 combined)
    if options.ErrorEnvelope:
        cal = CalibrationManager()
        cal.plotErrorEnvelope(calPeriodLow=1,calPeriodHigh=12,PMT=0)

    

    ## Makes file holding all the residuals for each PMT for each run which is to be used
    if options.makeGlobalResiduals:
        cal = CalibrationManager()
        runPeriods = [16,17,18,19,20,21,22,23,24]#,6,7,8,9,10,11,12]#[16,17,18,19,20,21,22,23,24]#,#[[1],[2],[3],[4],[5],[6],[7],[8],[9],[10],[11],[12]]
        
        cal.makeGlobalResiduals(runPeriods)


    #### All the steps for completely replaying runs (without doing a new calibration or new position maps along the way)
    if 0:
        rep = CalReplayManager()
        cal = CalibrationManager()
        runPeriods = [21,22]#,17,18,19,20,21,22,23,24]#,16,19,20,21,22,23,24]#,16,17,18,19,20,21,22,23,24]#[11,12]#,4,5,6,7,8,9,10,11,12]#[13,14,16,17,18,19,20,21,22,23,24]# 
        for runPeriod in runPeriods:
            #rep.makeBasicHistograms(runPeriod, sourceORxenon="source")
           
            rep.findPedestals(runPeriod)
            rep.runReplayPass1(runPeriod)
            rep.runGainBismuth(runPeriod)
            rep.findTriggerFunctions(runPeriod)
            rep.runReplayPass2(runPeriod)
            #cal.fitSourcePositions(runPeriod)
            
        
    
    ### Source Run Calibration Steps...
    ### 13,14,15 all bad!
    if 0: 
        runPeriods = [24,23]#[1,2,3,4,5,6,7,8,9,10,11,12]#[16,20,21,22,24,23]#[16,17,18,19,20,21,22,23,24]#[1,12]##[13,14,16,17,18,19,20,21,22,23,24]#
        rep = CalReplayManager()
        cal = CalibrationManager()

        iterations = 1 # number of times to run through the calibration

        for i in range(0,iterations,1):
        

            for runPeriod in runPeriods:

                # Calculate new linearity curves and nPE/keV values from previous iterations peaks
                if i>0:
                    cal.calc_new_nPE_per_keV(runPeriod) # compare widths of simulated peaks and data peaks to make new alphas
                    cal.LinearityCurves(runPeriod) # Calculate new Linearity Curves using new peak values
            

                # Data Stuff

                cal.runSourceCalReplayPeakFitter(runPeriod,doOnlyCeSnBi=True);
                cal.makeSourceCalibrationFile(runPeriod, Simulation=False, InEnergy=False) # gather source peak information in ADC
                cal.makeSourceCalibrationFile(runPeriod, Simulation=False, InEnergy=True) # gather source peak information in Energy

                #Simulation Stuff

                rep.runReverseCalibration(runPeriod) #Apply detector response model to simulation
                cal.fitSimSourcePeaks(runPeriod) #fit the source peaks in eta*Evis
                cal.makeSourceCalibrationFile(runPeriod, Simulation=True, InEnergy=False) #gather source peak information in eta*Evis
                cal.makeSourceCalibrationFile(runPeriod, Simulation=True, InEnergy=True)  #gather source peak information in Energy
                
                cal.calculateResiduals(runPeriod) # compare data peaks to simulated peaks


                # Calculate new linearity curves and nPE/keV values from previous iterations peaks
                if 1:#i<(iterations-1):
                    cal.calc_new_nPE_per_keV(runPeriod) # compare widths of simulated peaks and data peaks to make new alphas
                    cal.LinearityCurves(runPeriod) # Calculate new Linearity Curves using new peak values
            


    ### Replaying Xe Runs. Note that the position maps are calculated post replayPass2 and only need to
    ### be done once unless fundamental changes to the code are made upstream
    if 0: 
        runPeriods = [9]#,3,4,5,7]#[2,3,4,5,7,8,9,10]#,3,4,5,7] #[8,9,10]##### 1-7 are from 2011/2012, while 8-10 are from 2012/2013
        rep = CalReplayManager()
        cal = CalibrationManager()
        #cal.calc_nPE_per_PMT(runAllRefRun=False,writeNPEforAllRuns=True)
        for runPeriod in runPeriods:    
            #rep.makeBasicHistograms(runPeriod, sourceORxenon="xenon")
            #rep.findPedestals(runPeriod, sourceORxenon="xenon")
            rep.runReplayPass1(runPeriod, sourceORxenon="xenon")
            #rep.runGainBismuth(runPeriod, sourceORxenon="xenon")
            rep.runReplayPass2(runPeriod, sourceORxenon="xenon")
            #rep.runReplayPass3(runPeriod, sourceORxenon="xenon")
            #rep.runReplayPass4(runPeriod, sourceORxenon="xenon")

            
    
