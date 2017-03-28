#!/usr/bin/python

import os
import shutil

group = ["Octet"]#, "Quartet", "Pair"]
anaChoices = [1,2,3,4,5,6,7,8,9]
for octet in range(60,121+1,1):
    for g in group:
        for anaCh in anaChoices:
            #filenameRM = "lowStats_corrections/PureOverProc_%s-%i_Analysis-%i.txt"%(g,octet,anaCh)
            filename = "DeltaExp_OctetByOctetCorrections/ThOverProc_%s-%i_Analysis-%i.txt"%(g,octet,anaCh)
            #filename = "lowStats_corrections/ThOverProc_%s-%i_Analysis-%i.txt"%(g,octet,anaCh)
            
            #print filename
            if os.path.isfile(filename):
                outdir = "%s/Octet_%i/%sAsymmetry/Systematics/"%(os.getenv("ANALYSIS_RESULTS"),octet,g)
                #os.system("rm %s/%s"%(outdir,filenameRM))
                shutil.copy(filename,outdir)
                
                print "%s to %s"%(filename,outdir)



