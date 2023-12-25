#!/usr/bin/python

import sys
import math
from gurobipy import *
from datetime import date
import time

if len(sys.argv) < 2:
    print('Usage: myguropi.py lpfilename')
    exit(0)

logfile = 'mygurobi.log'
if len(sys.argv) >= 3:
    logfile = sys.argv[2]
print('logfile:', logfile)

log = open(logfile,"w")

today = date.today()


# Read and solve model

model = read(sys.argv[1])

log.write('Solving %s\n' % sys.argv[1])
log.write("variables = " + str(model.NumVars) + "\n")
log.write("constraints = " + str(model.NumConstrs) + "\n")

model.params.Cuts = 0
#model.params.presolve = 0
#model.params.NonConvex = 2

time1 = time.time()


model.optimize()

runtime = time.time() - time1


success = False

if model.status == GRB.status.INF_OR_UNBD:
    #logger.joint('->LP infeasible or unbounded\n')
    model.Params.DualReductions = 0
    t0p = time.time()
    model.optimize()
    t1p = time.time()

    runtime += t1p - t0p
if model.status == GRB.status.INFEASIBLE:
    #logger.joint('->LP infeasible')
    stringvalue = 'INFEASIBLE'
    model.computeIIS()
    model.write("model.ilp")
    success = False
        
elif model.status == GRB.status.UNBOUNDED:
    #logger.joint('->LP unbounded')
    stringvalue = 'UNBOUNDED'
elif model.status == GRB.OPTIMAL:
    #logger.joint(' OPTIMAL for {:s}, optimal value: {:.4e} in time {:.4e}'.format(signature, model.objval, runtime))
    success = True
    stringvalue = 'OPTIMAL'

elif model.status == GRB.SUBOPTIMAL:
    #logger.joint(' SUBOPTIMAL for {:s}, (sub)optimal value: {:.4e} in time {:.4e}'.format(signature, model.objval, runtime))
    success = True
    stringvalue = 'SUBOPTIMAL'
else:
    #logger.joint(' unexpected status ' + str(model.status))
    stringvalue = 'UNEXPECTED'

if success:
    log.write('Optimal objective = %g\n' % model.objVal)

    model.printQuality()    


    count = 0
    thresh = -1
    for v in model.getVars():
        if thresh < 0 or math.fabs(v.x) > thresh:
            count += 1

    print("\n" + str(count) + " nonzero variables at threshold " +str(thresh) +" in solution\n")        
    log.write("\n"+str(count) + " nonzero variables at threshold " +str(thresh) +" in solution\n")
    for v in model.getVars():
        if thresh < 0 or math.fabs(v.x) > min(thresh,1e-16):
            print( '%s = %.16e'%(v.varname,v.x))
            log.write( '%s = %.16e\n'%(v.varname,v.x))        

log.write("bye.\n")
log.close()
