#!/usr/bin/python

import sys
import os
from math import *
#sys.path.append(".")
from EnergyUncertainty.EnergyErrors import *

# I need to read in statistical uncertainties to store in the 



def getEnergyBin(energy):
    return int(energy)/10
    
def getBinEnergyMid(b):
    return b*10. + 5.

def getBinEnergyLowEdge(b):
    return b*10.

def getBinEnergyUpperEdge(b):
    return (b+1.)*10.

def weightRealStats(data,stats,e0,e1):
        """Weight the given data array by the real experimental statistics"""

        numer = sum([(1./stats[i])**2 * data[i] 
                     for i in range(getEnergyBin(e0),getEnergyBin(e1)) if e0<=getBinEnergyMid(i)<=e1])
        denom = sum([(1./stats[i])**2  
                     for i in range(getEnergyBin(e0),getEnergyBin(e1)) if e0<=getBinEnergyMid(i)<=e1])
        
        return numer/denom

def totalStatErr(data,dataErr,e0,e1):

    if len(data)!=len(dataErr):
        return 0

    numer = sum([(1./dataErr[i])**2 * data[i] 
                 for i in range(getEnergyBin(e0),getEnergyBin(e1)) if e0<=getBinEnergyMid(i)<=e1])
    weightsum = sum([(1./dataErr[i])**2 
                     for i in range(getEnergyBin(e0),getEnergyBin(e1)) if e0<=getBinEnergyMid(i)<=e1])
    A = numer/weightsum if weightsum>0. else 1.
    return fabs( ( 1./sqrt(weightsum) ) / A ) if weightsum>0. and fabs(A)>0. else 1.

def sumErrors(err):
    return sqrt(sum([er**2 for er in err]))

def getAsymmetry(anaCh,octLow,octHigh,elow,ehigh,sim=True):
    os.system("./MC_Corrections/AsymmCalculator.exe %s %i %i %f %f %i false > log.txt"%(anaCh,octLow, octHigh, elow,ehigh,int(sim)))
    infile = open("asymm_hold.txt",'r')
    asymm=[0,0]
    if infile:
        for line in infile:
            entry = line.split()
            asymm[0] = float(entry[1])
            asymm[1] = float(entry[2])
                #print("type0sim %f"%type0sim)
        infile.close()
    return asymm


class uncertaintyHandler:
    
    def __init__(self,year,anaChoice="C"):
        self.anaChoice = anaChoice
        self.stats = []
        self.A0 = 0.
        self.statUncert = 0.
        self.year = 0
        self.octLow = 0
        self.octHigh = 0
        self.realA = []
        self.realAerr = []
        self.stat_percent_err = []
        self.energy_err = []
        self.theory_err = []
        self.theory_corr = []
        self.systA = []
        self.syst_err = []
        self.syst_corr = []

        self.year = year
        
        if self.year == 2011:
            self.octLow = 0
            self.octHigh = 59
            
        elif self.year == 2012:
            self.octLow = 60
            self.octHigh = 121
        else:
            print('bad year given!!')
            exit(0)

    def statUncertainties(self):
        """Reads in the Statistical Uncertainty of every bin, calculates and 
        stores the total fitted asymmetry and the total percent error
        from statistics alone"""
        
        infile = open( os.environ["ANALYSIS_RESULTS"]+"/Asymmetries/"+
                       "UnCorr_OctetAsymmetries_AnaCh%s_Octets_%i-%i_BinByBin.txt"
                       %(self.anaChoice,self.octLow,self.octHigh),'r' )
        A_en = []
        A = []
        Aerr = []
        percErr = []

        if infile:
            for line in infile:
                l = line.split()
                if float(l[0])<1000.:
                    A_en.append(float(l[0]))
                    A.append(float(l[1]))
                    Aerr.append(float(l[2]))
                    if float(l[1])!=0.: 
                        percErr.append( fabs( float(l[2])/float(l[1]) ) )
                    
                    else:
                        percErr.append( 1. )

        else: 
            print("Couldn't open file for statistical uncertainties")
            exit(0)
        
        self.realA = A
        self.realAerr = Aerr 
        self.stat_percent_err = percErr


    def readMCSystematicCorrections(self):
        """Read in the systematic correction bin-by-bin for final
        asymmetries, calculates delta_A/A"""

        A_corr = []

        # begin with making sure uncorrected asymm is loaded
        if len(self.stat_percent_err)==0: 
            self.statUncertainties() 

         
        # Now read in corrected asymmetry
        infile = open( os.environ["ANALYSIS_RESULTS"]+"/Asymmetries/"+
                       "DeltaExpOnly_OctetAsymmetries_AnaCh%s_Octets_%i-%i_BinByBin.txt"
                       %(self.anaChoice,self.octLow,self.octHigh),'r' )

        if infile:
            for line in infile:
                l = line.split()
                if float(l[0])<1000.:
                    A_corr.append(float(l[1]))
                    
        else: 
            print("Couldn't open file for statistical uncertainties")
            exit(0)

        infile.close()

        
        # Calculate the percent error on the correction
        corr = []
        percErr = []
        
        for i in range( 0, len(A_corr) ):
            correction = 100.
            if fabs(self.realA[i])>0.:
                correction = A_corr[i]/self.realA[i] - 1. 
            
            corr.append(correction)
            percErr.append( fabs(correction*0.2) )
            
        self.systA = A_corr
        self.syst_err = percErr
        self.syst_corr = corr

    def calcMCSystematicError(self,elow,ehigh,BStypes=[]):

        percErrType0Sim = 0

        t0errs = []
       
        infile = None

        #read in Type 0 MC error from Type 1 shuffle
        os.system("./MC_Corrections/EventShuffle.exe %i %i %f %f -0.1 0.0 0.0 0.0 > log.txt"%(self.octLow, self.octHigh, elow,ehigh))
        infile = open("resultHolderEventShuffle.txt",'r')
    
        if infile:
            for line in infile:
                if line[0]=="p":
                    entry = line.split()
                    t0errs.append(float(entry[1]))
                    #percErrType0Sim = float(entry[1])
                    #print("perc Err on T0 %f"%percErrType0Sim)

        infile.close()
        
        #read in Type 0 MC error from Type 2 shuffle
        os.system("./MC_Corrections/EventShuffle.exe %i %i %f %f 0.0 -0.25 0.0 0.0 > log.txt"%(self.octLow, self.octHigh, elow,ehigh))
        infile = open("resultHolderEventShuffle.txt",'r')
    
        if infile:
            for line in infile:
                if line[0]=="p":
                    entry = line.split()
                    t0errs.append(float(entry[1]))
                    #percErrType0Sim = float(entry[1])
                    #print("perc Err on T0 %f"%percErrType0Sim)

        infile.close()

        #read in Type 0 MC error from Type 3 shuffle
        os.system("./MC_Corrections/EventShuffle.exe %i %i %f %f 0.0 0.0 -0.05 0.0 > log.txt"%(self.octLow, self.octHigh, elow,ehigh))
        infile = open("resultHolderEventShuffle.txt",'r')
    
        if infile:
            for line in infile:
                if line[0]=="p":
                    entry = line.split()
                    t0errs.append(float(entry[1]))
                    #percErrType0Sim = float(entry[1])
                    #print("perc Err on T0 %f"%percErrType0Sim)

        infile.close()
        
        #read in Type 0 MC error from Type 4 shuffle
        os.system("./MC_Corrections/EventShuffle.exe %i %i %f %f 0.0 0.0 0.0 0.25 > log.txt"%(self.octLow, self.octHigh, elow,ehigh))
        infile = open("resultHolderEventShuffle.txt",'r')
    
        if infile:
            for line in infile:
                if line[0]=="p":
                    entry = line.split()
                    t0errs.append(float(entry[1]))
                    #percErrType0Sim = float(entry[1])
                    #print("perc Err on T0 %f"%percErrType0Sim)

        infile.close()
        
        
        percErrType0Sim = sqrt( sum(err*err for err in t0errs) )
        #print(percErrType0Sim)
       
        asymm = []
        A = []
        deltaA_MC = [] 
        deltaA_Stat = []

        asymm = getAsymmetry("D",self.octLow,self.octHigh,elow,ehigh,sim=False)
       
        A.append( asymm[0] )
        deltaA_MC.append( fabs(asymm[0]*percErrType0Sim) )
        deltaA_Stat.append( fabs(asymm[1]) )
        print("D asymm: %f +/- %f +/- %f"%(A[0],asymm[1],deltaA_MC[0]))

        for bs in BStypes:
            asymm = getAsymmetry(bs,self.octLow,self.octHigh,elow,ehigh,sim=False)
            print("%s asymm: %f +/- %f +/- %f(%f)"%(bs,asymm[0],asymm[1],fabs(asymm[0]-A[0]),
                                                    sqrt(asymm[1]**2+deltaA_Stat[0]**2))
                  )
            A.append( asymm[0] )
            deltaA_MC.append( fabs(asymm[0]-A[0]) )
            deltaA_Stat.append( fabs(asymm[1]) )
            
        

        deltaA_MCFinal = sqrt(sum( (deltaA_MC[i]/deltaA_Stat[i]**2)**2 for i in range(0,len(A)) ) ) / sum( 1/err**2 for err in deltaA_Stat )
        A_final = sum( A[i]/deltaA_Stat[i]**2 for i in range(0,len(A)) ) / sum( 1/err**2 for err in deltaA_Stat )

        return fabs(deltaA_MCFinal/A_final)
        

        


    def readEnergyUncertainties(self):
        """Read in the Energy uncertainties and store the percent 
        error from every bin."""
        
        infile = open( os.environ["ANALYSIS_RESULTS"]+"/Corrections/"+
                       "EnergyLinearityUncertainty_%i.txt"%self.year)

        percErr = []

        if infile:
            for line in infile:
                if line[0]!="#":
                    l = line.split()
                    percErr.append( fabs(float(l[3])) )

        else: 
            print("Couldn't open file for Energy uncertainties")
            exit(0)
        

        self.energy_err = percErr

    def makeSystematicCorrections(self,errFac=0.2):
        """Read in the systematic correction bin-by-bin for final
        asymmetries, calculates delta_A/A, and then multiplies by the 
        error factor before filling the percent error for every bin."""

        A_uncorr = []
        A_corr = []

        # begin with uncorrected asymmetry
        infile = open( os.environ["ANALYSIS_RESULTS"]+"/Asymmetries/"+
                       "UnCorr_OctetAsymmetries_AnaCh%s_Octets_%i-%i_BinByBin.txt"
                       %(self.anaChoice,self.octLow,self.octHigh),'r' )

        if infile:
            for line in infile:
                l = line.split()
                if float(l[0])<1000.:
                    A_uncorr.append(float(l[1]))

        else: 
            print("Couldn't open file for statistical uncertainties")
            exit(0)

         

        infile.close()

        # Now read in corrected asymmetry
        infile = open( os.environ["ANALYSIS_RESULTS"]+"/Asymmetries/"+
                       "DeltaExpOnly_OctetAsymmetries_AnaCh%s_Octets_%i-%i_BinByBin.txt"
                       %(self.anaChoice,self.octLow,self.octHigh),'r' )

        if infile:
            for line in infile:
                l = line.split()
                if float(l[0])<1000.:
                    A_corr.append(float(l[1]))

        else: 
            print("Couldn't open file for statistical uncertainties")
            exit(0)

        infile.close()

        #for x in A_corr:
        #    print(x)

        # Calculate the percent error on the correction
        corr = []
        percErr = []
        
        for i in range( 0, len(A_corr) ):
            correction = 100.
            if fabs(A_uncorr[i])>0.:
                correction = A_corr[i]/A_uncorr[i] - 1. 
            
            corr.append(correction)
            percErr.append( fabs(correction*errFac) )

        self.syst_err = percErr
        self.syst_corr = corr


    


    def makeTheoryUncertainties(self,lambda_err=0.):
        """Calculate the theory uncertainties given the error in lambda provided"""
        percErr = []
        
    def calcEnergyUncert(self,emin,emax):

        if len(self.energy_err)==0: 
            self.readEnergyUncertainties() 
            
        if len(self.stat_percent_err)==0: 
            self.statUncertainties() 
            
        return weightRealStats(self.energy_err,self.stat_percent_err,emin,emax)

    def calcStatUncert(self,emin,emax):

        #if len(self.stat_percent_err)==0: 
        #    self.statUncertainties() 
            
        #return totalStatErr(self.realA,self.realAerr,emin,emax) 
        asymm = getAsymmetry(self.anaChoice,self.octLow,self.octHigh,emin,
                             emax,False)
        self.A0 = asymm[0]
        return fabs( asymm[1]/asymm[0] )

    def calcSystematicUncert(self,emin,emax):

        if len(self.syst_err)==0: 
            self.makeSystematicCorrections() 
            
        if len(self.stat_percent_err)==0: 
            self.statUncertainties() 
            
        return weightRealStats(self.syst_err,self.stat_percent_err,emin,emax)

    def gainUncert(self,emin,emax):

        if len(self.stat_percent_err)==0: 
            self.statUncertainties() 

        infile = open( os.environ["ANALYSIS_RESULTS"]+"/Corrections/"+
                       "gainUncertainty_%i.txt"%self.year)

        gainErr = []

        if infile:
            for line in infile:
                if line[0]!="#":
                    l = line.split()
                    gainErr.append( fabs(float(l[3])) )
                    #print(float(l[3]))

        else: 
            print("Couldn't open file for gain uncertainties")
            exit(0)
        
            
        return weightRealStats(gainErr,self.stat_percent_err,emin,emax)



    def minimizer(self):
        """Minimizes the error and outputs all errors in order or largest to smallest"""
        
        BStypes = []

        if self.anaChoice=="A":
            BStypes = ["F","G"]
        elif self.anaChoice=="B":
            BStypes = ["F"]
        elif self.anaChoice=="C":
            BStypes = ["F","J","K"]

        self.statUncertainties()
        self.readEnergyUncertainties()
        #self.makeSystematicCorrections()
        #self.readMCSystematicCorrections()

        ofile = open("%i_Uncertainties.txt"%self.year, 'w')
        

        minErr = 100.
        minEnergyErr = 0.
        minStatisticalErr = 0.
        minSystematicErr = 0.
        enBinLow = None
        enBinHigh = None

        minNumBins = 40
        lowbound = 10
        upperbound = 76
        for lowBin in range(lowbound,upperbound+1-minNumBins):
            for highBin in range(lowBin+minNumBins,upperbound+1):

                Errors = []

                Errors.append( self.calcEnergyUncert(     getBinEnergyMid(lowBin),getBinEnergyMid(highBin) ) )
                Errors.append( self.calcStatUncert(       getBinEnergyLowEdge(lowBin),getBinEnergyUpperEdge(highBin) ) )
                
                ##### Systematic error handling

                # calcSystematicUncert applies a uniform percent uncertainty to every bin correction
                #Errors.append( self.calcSystematicUncert( getBinEnergyMid(lowBin),getBinEnergyMid(highBin) ) )
                
                # Ignores systematics in minimization
                #Errors.append(0.)

                # Multiplies final correction by some fractional percent uncertainty ( 20% in 2010 )
                #Errors.append( fabs( weightRealStats(self.syst_corr,self.stat_percent_err,getBinEnergyMid(lowBin),getBinEnergyMid(highBin)) ) * 0.2)

                Errors.append( self.calcMCSystematicError(getBinEnergyLowEdge(lowBin),
                                                          getBinEnergyUpperEdge(highBin),BStypes) );

                errSum = sumErrors(Errors)
                


                print("****************************")
                print("%i-%i keV total Error = %f\n"%(getBinEnergyLowEdge(lowBin),getBinEnergyUpperEdge(highBin),errSum))
                print("Asymm = %f\n"%self.A0)
                print("Weighted Energy Error: %f"%Errors[0])
                print("Total Statistical Error: %f"%Errors[1])
                print("Total Systematic Error: %f"%Errors[2])

                
                ofile.write("****************************\n")
                ofile.write("Asymm= %f\n"%self.A0)
                ofile.write("%i-%i keV_total_Error= %f\n\n"%(getBinEnergyLowEdge(lowBin),getBinEnergyUpperEdge(highBin),errSum))
                ofile.write("Weighted_Energy_Error: %f\n"%Errors[0])
                ofile.write("Total_Statistical_Error: %f\n"%Errors[1])
                ofile.write("Total_Systematic_Error: %f\n"%Errors[2])

                
                if errSum < minErr:
                    minErr = errSum
                    enBinLow = lowBin
                    enBinHigh = highBin
                    minEnergyErr = Errors[0]
                    minStatisticalErr = Errors[1]
                    minSystematicErr = Errors[2]


        ofile.close()

       # minSystCorrTotal = weightRealStats(self.syst_corr,self.stat_percent_err,getBinEnergyMid(enBinLow),getBinEnergyMid(enBinHigh))
        print
        #print("Systematic Correction: %f"%minSystCorrTotal)
        print
        print("Weighted Energy Error: %f"%minEnergyErr)
        print("Total Statistical Error: %f"%minStatisticalErr)
        print("Total Systematic Error: %f"%minSystematicErr)
        print("Min total Err = %f"%minErr)
        print("Energy Range = %f-%f"%(getBinEnergyLowEdge(enBinLow),getBinEnergyUpperEdge(enBinHigh)))

if __name__ == "__main__":
    
    uncert = uncertaintyHandler(2011,"C")
    #uncert.minimizer()

    low = 230
    high = 750
    print(uncert.calcEnergyUncert(low,high))
    print(uncert.calcStatUncert(low,high))
    #print(uncert.calcSystematicUncert(low,high))
    print(uncert.calcMCSystematicError(low,high,["F","J","K"]) )
    asymm = getAsymmetry(uncert.anaChoice,uncert.octLow,uncert.octHigh,low,high,False)
    print(sumErrors([uncert.calcEnergyUncert(low,high),fabs( asymm[1]/asymm[0] ),uncert.calcMCSystematicError(low,high,["F","J","K"])]))
    print("A: %f +/- %f"%(asymm[0],asymm[1]))

    print("Gain Err: +/- %f"%uncert.gainUncert(low,high))
    #220-750 for 2011 : 0.00612 
    
    #290-750 for 2012 : 0.01114

    
