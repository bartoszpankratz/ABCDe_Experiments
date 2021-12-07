#!/usr/bin/env python
# coding: utf-8
# gen_LFR.py n gamma min_deg degs beta coms mu xi path no_iters
import sys
import os.path

n = int(sys.argv[1])
gamma = float(sys.argv[2])
min_deg = int(sys.argv[3])
degs = sys.argv[4]
degs = open(degs,'r').read().splitlines()
degs = [int(d) for d in degs]
beta = float(sys.argv[5])
coms = sys.argv[6]
coms = open(coms,'r').read().splitlines()
coms = [int(c) for c in coms]
mu = float(sys.argv[7])
xi = float(sys.argv[8])
res_dir = sys.argv[9]
iters = int(sys.argv[10])

import random
import logging
import networkit as nk

if mu > 1 or mu < 0:
    logging.basicConfig(filename = logname)
    logging.warning("mu equal to: " + str(mu))
    logging.warning("for parameters: " + str(n) + ", " + str(xi) + ", " + str(beta) + ", "  + str(gamma) + ", " + str(min_deg))
    raise Exception("mu equal to: " + str(mu)) 
                  
for iter in range(iters):
    random.seed(iter + 22)
    logname = "logs/LFR_" +str(iter + 1) + "_" + str(n) + "_" + str(xi) + "_" + str(beta) + "_" + str(gamma) + "_" + str(min_deg) + "_.log"
    LFRname = res_dir + "/LFR_" + str(iter + 1) + "_" + str(n) + "_" + str(xi) + "_" + str(beta) + "_"  + str(gamma) + "_" + str(min_deg) + "_.dat"
    comname = res_dir + "/LFR_clusters_" + str(iter + 1) + "_" + str(n) + "_" + str(xi) + "_" + str(beta) + "_"  + str(gamma) + "_" + str(min_deg) + "_.dat"
    try:
        lfr = nk.generators.LFRGenerator(n)
        lfr.setDegreeSequence(degs)
        lfr.setCommunitySizeSequence(coms)
        lfr.setMu(mu)
        lfrG = lfr.generate()

        with open(LFRname, 'w') as f:
            for ine,oute in lfrG.iterEdges():
                print(ine+1,"\t",oute+1, file = f)
        f.close()
        
        with open(comname, 'w') as f:
            for i,c in enumerate(lfr.getPartition().getVector()):
                print(i+1,"\t",c+1, file = f)
        f.close()
        
    except Exception as err:
        logging.basicConfig(filename = logname)
        logging.warning(err)
        logging.warning("for parameters: " + str(n) + ", " + str(xi) + ", " + str(beta) + ", "  + str(gamma) + ", " + str(min_deg))
