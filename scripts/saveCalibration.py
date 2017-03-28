#!/usr/bin/python

import os

calibrationName = "03272017_5umSrcFoil_it1"#"07312016_3PoissonMethod_with_Trigger_noBi2_linear_it4"
year = "2012-2013" #"2011-2012" or "2012-2013"
calibrationName = calibrationName+"_"+year
location = "/extern/UCNA/goodCalibrations/" + calibrationName + "/"

os.system("mkdir -p %s"%(location))

calPeriods = None

if year=="2011-2012":
    calPeriods = [1,2,3,4,5,6,7,8,9,10,11,12]
else:
    calPeriods = [13,14,15,16,17,18,19,20,21,22,23,24]

for per in calPeriods:

    os.system("cp %s/linearity_curves/linCurves_SrcPeriod_%i.pdf %s"%(os.getenv("ANALYSIS_CODE"),per,location))
    os.system("cp %s/lin_curves_srcCal_Period_%i.dat %s"%(os.getenv("LINEARITY_CURVES"),per,location))
    os.system("cp %s/simulation_comparison/nPE_per_keV/nPE_per_keV_%i.dat %s"%(os.getenv("ANALYSIS_CODE"),per,location))
    os.system("cp %s/simulation_comparison/nPE_per_keV/width_comp_%i.pdf %s"%(os.getenv("ANALYSIS_CODE"),per,location))
    

if year=="2011-2012":
    os.system("cp %s/error_envelope/error_envelope_calPeriods_1-12.dat %s"%(os.getenv("ANALYSIS_CODE"),location))
    os.system("cp %s/residuals/global_residuals_Erecon_runPeriods_1-12.dat %s"%(os.getenv("ANALYSIS_CODE"),location))
 
if year=="2012-2013":
    os.system("cp %s/error_envelope/error_envelope_calPeriods_16-24.dat %s"%(os.getenv("ANALYSIS_CODE"),location))
    os.system("cp %s/residuals/global_residuals_Erecon_runPeriods_16-24.dat %s"%(os.getenv("ANALYSIS_CODE"),location))
    
os.system("cp final_residuals.pdf %s"%location)


