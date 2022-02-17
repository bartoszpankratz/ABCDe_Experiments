#!/usr/bin/env python
# coding: utf-8
# gen_LFR.py threads n gamma min_deg E_deg degs beta coms mu xi path seed
import sys
import os.path

threads = int(sys.argv[1])
n = int(sys.argv[2])
gamma = float(sys.argv[3])
min_deg = int(sys.argv[4])
E_deg = float(sys.argv[5])
degs = sys.argv[6]
degs = open(degs,'r').read().splitlines()
degs = [int(d) for d in degs]
beta = float(sys.argv[7])
coms = sys.argv[8]
coms = open(coms,'r').read().splitlines()
coms = [int(c) for c in coms]
mu = float(sys.argv[9])
xi = float(sys.argv[10])
res_dir = sys.argv[11]

logname = "logs/LFR_" + str(threads) + "_" + str(n) + "_" + str(xi) + "_" + str(beta) + "_" + str(gamma) + "_" + str(min_deg) + ".log"
LFRname = res_dir + "/LFR_" + str(threads) + "_" + str(n) + "_" + str(xi) + "_" + str(beta) + "_"  + str(gamma) + "_" + str(min_deg) + ".dat"

if os.path.isfile(LFRname) or os.path.isfile(logname):
    sys.exit() 																	
import random
import logging
import networkit as nk
from timeit import default_timer as timer

random.seed(int(sys.argv[12]))

nk.setNumberOfThreads(threads)

if mu > 1 or mu < 0:
    logging.basicConfig(filename = logname)
    logging.warning("mu equal to: " + str(mu))
    logging.warning("for parameters: " + str(threads) + ", " + str(n) + ", " + str(xi) + ", " + str(beta) + ", "  + str(gamma) + ", " + str(min_deg))
    raise Exception("mu equal to: " + str(mu)) 
                  

try:
    start = timer()
    lfr = nk.generators.LFRGenerator(n)
    #lfr.generatePowerlawDegreeSequence(min_deg, int(n**0.5), -gamma)
    lfr.setDegreeSequence(degs)
    #lfr.generatePowerlawCommunitySizeSequence(int(0.005*n), int(0.2*n), -beta)
    lfr.setCommunitySizeSequence(coms)
    lfr.setMu(mu)
    lfrG = lfr.generate()
    time = timer() - start
    n_edges = lfrG.numberOfEdges()
    res = [n, xi, mu, beta, gamma, min_deg, E_deg, n_edges, "NaN", threads, time]
    with open(LFRname, 'w') as f:
        for r in res:
            print(r, file = f)
    f.close()
except Exception as err:
    logging.basicConfig(filename = logname)
    logging.warning(err)
    logging.warning("for parameters: " + str(threads) + ", " + str(n) + ", " + str(xi) + ", " + str(beta) + ", "  + str(gamma) + ", " + str(min_deg))
