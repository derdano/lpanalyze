# -*- coding: utf-8 -*-
import gurobipy as gp
import numpy as np
from gurobipy import *
from numpy import linalg as LA
import sys
import struct
from log import danoLogger
from myutils import breakexit, askfordouble
import time

def simplebreak():
    stuff = input('break > ')
    if stuff == 'x' or stuff == 'q':
        sys.exit("bye")


def ladoit(alldata):
    log = alldata['log']

    global MYDEBUG
    MYDEBUG = False

    tstart = time.time()

    log.joint('\nAnalyzing %s.\n' % alldata['LPFILE'])

    model = alldata['model']

    
    log.joint("Q constraints = " + str(model.NumQConstrs)+'.\n')    
    Qconstrs = model.getQConstrs()
    LOUD = False

    code = 0

    #first check
    for Qconstr in Qconstrs:
        therow = model.getQCRow(Qconstr)
        if LOUD:
            log.joint(Qconstr.QCName+ ' ' + Qconstr.QCSense + ' ' + str(Qconstr.QCRHS) + '\n')
            #print('therow', therow)

        linQc  = therow.getLinExpr()
        rhs = Qconstr.qcrhs

        if rhs != 0 or Qconstr.QCSense != '<':
            log.joint(Qconstr.QCName + ' fails 1.\n')
            code = 1
            break

    log.joint('First test passed.\n')

 

    code = 0
    varlindeg = alldata['varlindeg']
    varlindeg_bin = {}
    seen = {}


    distcount = 0  # number of distance sets    
    distvarset = {}   # variables in current distance set
    distsize = {}  # number of vars in current distance set
    distvar_owner = {}  # binary variable owner for current distance set
    distvar_partner = {}  # continuous variable multiplied times owner in current distance set
    distconstr = {} # distance constraint as we index it

    for var in model.getVars():
        if var.vtype == GRB.BINARY:
            varlindeg_bin[var.varname] = 0
        seen[var.varname] = 0
    

    quantsquares = 0

    LOUD = False

    for Qconstr in Qconstrs:
        therow = model.getQCRow(Qconstr)
        if LOUD:
            log.joint(Qconstr.QCName+ ' ' + Qconstr.QCSense + ' ' + str(Qconstr.QCRHS) + '\n')
            print('therow', therow)

        linQc  = therow.getLinExpr()
        rhs = Qconstr.qcrhs

        #print(linQc)

        if linQc.size() > 0:
            log.joint(Qconstr.QCName + ' fails 2.\n')
            code = 1
            break

        #simplebreak()

        distconstr[distcount] = Qconstr
        distvarset[distcount] = []

        numbilin = 0
        for j in range(therow.size()):
            v1 = therow.getVar1(j)
            v2 = therow.getVar2(j)


            #print('>>>','v1',v1.varname,'v2',v2.varname)
            coeff = therow.getCoeff(j)
            if v1.varname != v2.varname:
                if numbilin > 0:
                    log.joint(Qconstr.QCName+ ' fails bilin.\n')
                    print(v1.varname, v2.varname)
                    simplebreak()
                    code = 1
                    break
                numbilin += 1                
                if coeff != -1.0:
                    log.joint(Qconstr.QCName+ ' fails coeff.\n')                
                    code = 1
                    break
                if varlindeg[v1.varname] ==0 and v1.vtype != GRB.BINARY and v2.vtype == GRB.BINARY and varlindeg_bin[v2.varname]==0:
                    # good case
                    varlindeg[v1.varname] = 1
                    varlindeg_bin[v2.varname] = 1

                    distvar_owner[distcount] = v2
                    distvar_partner[distcount] = v1
                    
                elif varlindeg_bin[v1.varname] == 0 and v1.vtype == GRB.BINARY and varlindeg[v2.varname]== 0 and v2.vtype != GRB.BINARY:
                    # reverse good case
                    varlindeg_bin[v1.varname] = 1
                    varlindeg[v2.varname] = 1

                    distvar_owner[distcount] = v1
                    distvar_partner[distcount] = v2
                   
                else:
                    log.joint(Qconstr.QCName+ ' fails bindeg.\n')                
                    print(v1.varname, varlindeg[v1.varname], v2.varname, varlindeg[v2.varname])
                    code = 1
                    simplebreak()                    
                    break
            else: #v1 = v2
                distvarset[distcount].append(v1)
                
                if coeff != 1.0:
                    log.joint(Qconstr.QCName+ ' fails coeff=1.\n')                
                    print(v1.varname)
                    code = 1
                    simplebreak()                    
                    break
                if seen[v1.varname]:
                    log.joint(Qconstr.QCName+ ' fails seen.\n')                
                    print(v1.varname)
                    code = 1
                    simplebreak()                    
                    break
                seen[v1.varname] = 1
        distsize[distcount] = len(distvarset[distcount])

        if False:
            log.joint('distance set %d of size %d.\n'%(distcount, distsize[distcount]))
            for j in range(distsize[distcount]):
                log.joint(' %s'%(distvarset[distcount][j].varname))
            log.joint('\n')
        #print(distcount, distsize[distcount], distvarset[distcount])
        distcount += 1
            
        if code == 1:
            log.joint(Qconstr.QCName+ ' fails.\n')
            simplebreak()
            break
        else:
            quantsquares += therow.size()-1

        #simplebreak()

    log.joint('Final code: ' + str(code) + '.\n')
    log.joint('Quantsquares: %d, plus bilinears %d.  Total n = %d.\n'%(quantsquares, quantsquares + 2*(len(Qconstrs)), model.NumVars))


    alldata['distcount'] = distcount
    alldata['distsize'] = distsize
    alldata['distvarset'] = distvarset
    alldata['distvar_owner'] = distvar_owner
    alldata['distvar_partner'] = distvar_partner
    alldata['distconstr'] = distconstr

    tend = time.time()

    log.joint('Done analyzing in time %g\n'%(tend - tstart))
    

    laresetLP(alldata)

    lasolveLP(alldata)

    return code


def laresetLP(alldata):
    log = alldata['log']

    log.joint('Now building new LP.\n')

    distcount = alldata['distcount']
    distsize = alldata['distsize']
    distvarset = alldata['distvarset']
    
    model = alldata['model']

    #set up new objective
    distobj = 0
    for i in range(distcount):
        for j in range(distsize[i]):
            dvarij = distvarset[i][j]
            distobj += dvarij*dvarij
        '''
        log.joint('distance set %d of size %d.\n'%(i, distsize[i]))
        for j in range(distsize[i]):
            log.joint(' %s'%(distvarset[i][j].varname))
        log.joint('\n')
        '''
        #print('***',i,distobj)
        #simplebreak()
        
    model.setObjective(distobj,GRB.MINIMIZE)

    log.joint('Reset objective.\n')

    #now update distance-defining constraints

    #simplebreak()

    LOUD = False

    distvar_owner = alldata['distvar_owner'] 
    distvar_partner = alldata['distvar_partner']
    distconstr = alldata['distconstr']

    cheatcoeff = alldata['cheatcoeff'] = 125

    USELAZY = True

    if USELAZY == False:
        boundvar = {}
        ubval = cheatcoeff**.5

    for i in range(distcount):
        Qconstr = distconstr[i]
        owner = distvar_owner[i]
        partner = distvar_partner[i]
        therow = model.getQCRow(Qconstr)
    
        if LOUD:
            log.joint(Qconstr.QCName+ ' ' + Qconstr.QCSense + ' ' + str(Qconstr.QCRHS) + '\n')
            #print('therow', therow)
            #print('newrow', newrow)

        rhs = Qconstr.qcrhs
        model.remove(Qconstr)

        if USELAZY:
            newrow = therow - cheatcoeff*owner + partner*owner
            model.addConstr(newrow <= 0, name = Qconstr.QCName)

        else:
            # first, add new variable
            boundvar[i] = model.addVar(obj = 0.0, lb = 0, ub = cheatcoeff, name = "up_" + str(i))
            newrow = therow - boundvar[i] + partner*owner
            model.addConstr(newrow <= 0, name = Qconstr.QCName)

            model.addConstr(boundvar[i] <= ubval*owner, name = 'vplus'+str(i))


        #finally

        owner.vtype = GRB.CONTINUOUS

    if USELAZY:
        model.write('lazy.lp')
    else:
        model.write('notlazy.lp')
        simplebreak()


def lasolveLP(alldata):
    log = alldata['log']

    log.joint('Now solving reset LP.\n')

    distcount = alldata['distcount']
    distsize = alldata['distsize']
    distvarset = alldata['distvarset']
    distvar_owner = alldata['distvar_owner']
    
    model = alldata['model']

    model.optimize()

    log.joint('Solving %s' % sys.argv[1])
    log.joint("variables = " + str(model.NumVars))
    log.joint("constraints = " + str(model.NumConstrs))

    if model.status == GRB.status.INF_OR_UNBD:
        log.joint('->LP infeasible or unbounded')
        model.Params.DualReductions = 0
        model.optimize()
    if model.status == GRB.status.INFEASIBLE:
        log.joint('->LP infeasible')                                            
        model.computeIIS()
        model.write("model.ilp")
        sys.exit(0)

    elif model.status == GRB.status.UNBOUNDED:
        log.joint('->LP unbounded')                                             
        sys.exit(0)
    
    elif model.status == GRB.OPTIMAL:
        log.joint(' ->OPTIMAL\n')

    log.joint('Optimal objective = %g' % model.objVal)

    model.printQuality()

    # Now, get solution stats.

    cheatcoeff = alldata['cheatcoeff']
    distvalue2 = 0
    setvalue2 = np.zeros(distcount)

    numpositive = 0
    TOL = 1e-10

    LOUD = False
    for i in range(distcount):
        owneri = distvar_owner[i]
        #ub_sum += cheatcoeff*owneri.x

        
        setvalue2[i] = 0
        for j in range(distsize[i]):
            dvarij = distvarset[i][j]
            xvalue = dvarij.x
            if abs(xvalue) > TOL:
                setvalue2[i] += dvarij.x*dvarij.x
                if LOUD:
                    log.joint('Set %d, elt %d (%s) at value %g\n'%(i,j,dvarij.varname, dvarij.x))

        distvalue2 += setvalue2[i]
        numpositive += setvalue2[i] > TOL
        log.joint('Set %d sum-of-squares %g; ub %g\n'%(i, setvalue2[i], cheatcoeff*owneri.x))
        #simplebreak()

    log.joint('Positives: %d; sum total %g\n'%(numpositive, distvalue2))

    ind = np.argsort(-setvalue2)

    ordered = setvalue2[ind]

    '''
    print(ordered)
    print(np.sum(ordered[75:]))
    '''

    sumsmallest2 = 0
    cardthresh = 75
    for k in range(cardthresh, distcount):
        i = ind[k]
        owneri = distvar_owner[i]

        sumsmallest2 += ordered[k]

        log.joint('Ordered set %d is %d (%s, %g) at %g\n'%(k,i, owneri.varname, owneri.x, ordered[k]))


    log.joint('Sum of %d smallest: %g\n'%(distcount - cardthresh, sumsmallest2))
    log.joint('So overall lower bound: %g\n'%(model.objVal + sumsmallest2))
        

