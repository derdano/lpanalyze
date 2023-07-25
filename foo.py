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


def foo(alldata):
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

    #hard-coded, boo
    cheatcoeff = alldata['cheatcoeff'] = 125
    alldata['notmaxcard'] = 25

    
    #keep old model
    alldata['originalmodel_relax'] = alldata['model'].relax()
    alldata['originalmodel_relax'].write('foo.lp')
    #simplebreak()


    alldata['keyorder'] = np.zeros(distcount, dtype = int)



    
    laresetLPtoSumOfSquares(alldata)
    simplebreak()

    lasolveLP(alldata,'ss')
    simplebreak()

    alldata['original_ss'] = alldata['model'].objval

    doparametric = False

    if doparametric:
        laresetLPtoParametric(alldata, alldata['original_ss'])

        code = {}
        pushvalue = np.zeros(distcount)

        for i in range(distcount):    #(2):
            code[i], pushvalue[i] = lapushdownParametric(alldata, i)


        #sort

        log.joint('Sorted pushes:\n')
        ind = np.argsort(pushvalue)
        ordered = pushvalue[ind]
        for k in range(distcount):
            i = ind[k]
            log.joint('k %d i %d value %g\n'%(k,i, pushvalue[i]))
    

    return code


def laresetLPtoSumOfSquares(alldata):
    log = alldata['log']

    log.joint('Now building new LP with sum-of-squares objective.\n')

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

    log.joint('Reset objective to sum-of-squares.\n')

    #now update distance-defining constraints

    #simplebreak()

    LOUD = False

    distvar_owner = alldata['distvar_owner'] 
    distvar_partner = alldata['distvar_partner']
    distconstr = alldata['distconstr']
    cheatcoeff = alldata['cheatcoeff']


    USELAZY = False

    if USELAZY == False:
        #boundvar = {}
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
            for j in range(distsize[i]):
                dvarij = distvarset[i][j]

                model.addConstr(dvarij <= ubval*owner, name = 'vplus'+str(i)+','+str(j))
                model.addConstr(dvarij >= -ubval*owner, name = 'vminus'+str(i)+','+str(j))                
        
        '''
        else:
            # first, add new variable
            boundvar[i] = model.addVar(obj = 0.0, lb = 0, ub = cheatcoeff, name = "up_" + str(i))
            newrow = therow - boundvar[i] + partner*owner
            model.addConstr(newrow <= 0, name = Qconstr.QCName)

            model.addConstr(boundvar[i] <= ubval*owner, name = 'vplus'+str(i))
        '''

        #finally

        owner.vtype = GRB.CONTINUOUS

    model.reset()
    
    if USELAZY:
        model.write('lazy.lp')
    else:
        model.write('notlazy.lp')
        #simplebreak()


def lasolveLP(alldata, header):
    log = alldata['log']
    # Initialize a list to store x^*
    x_sol = [] 

    log.joint('Now solving LP. ')
    if header:
        log.joint('%s'%(header))
    log.joint('.\n')

    distcount = alldata['distcount']
    distsize = alldata['distsize']
    distvarset = alldata['distvarset']
    distvar_owner = alldata['distvar_owner']
    
    model = alldata['model']

    log.joint('Solving %s' % sys.argv[1])
    log.joint("variables = " + str(model.NumVars))
    log.joint("constraints = " + str(model.NumConstrs))

    model.optimize()

    code = 0

    if model.status == GRB.status.INF_OR_UNBD:
        log.joint('->LP infeasible or unbounded')
        model.Params.DualReductions = 0
        model.optimize()
    if model.status == GRB.status.INFEASIBLE:
        log.joint('->LP infeasible\n')                                            
        model.computeIIS()
        model.write("model.ilp")
        code = 1
        return code, 1e10  #bad: hard-coded
    elif model.status == GRB.status.UNBOUNDED:
        log.joint('->LP unbounded\n')                                             
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

    # new 
    for i in range(distcount):
        x_i = []
        for j in range(distsize[i]):
            dvarij = distvarset[i][j]
            xvalue = dvarij.x

            # if xvalue too small just append 0
            if abs(xvalue) <= TOL: x_i.append(0)
            if abs(xvalue) > TOL:
                # new
                x_i.append(xvalue)
        # new
        for j in range(len(x_i)):
            log.joint('%g\n'%(x_i[j]))
    log.joint('\n \n')

    # new
    log.joint('Minimize \n obj: \n  ')

    for i in range(distcount):
        owneri = distvar_owner[i]

        # new
        x_i = []
        # varij_name = []
        # for j in range(distsize[i]):
        #     x_i.append(distvarset[i][j])
        #ub_sum += cheatcoeff*owneri.x

        
        setvalue2[i] = 0
        for j in range(distsize[i]):
            dvarij = distvarset[i][j]
            xvalue = dvarij.x

            # new 
            # varij_name.append(dvarij.varname)

            # if xvalue too small just append 0
            if abs(xvalue) <= TOL: x_i.append(0)


            if abs(xvalue) > TOL:
                # new
                x_i.append(xvalue)

                setvalue2[i] += dvarij.x*dvarij.x
                if LOUD:
                    if header:
                        log.joint('%s '%(header))
                    log.joint('Set %d, elt %d (%s) at value %g\n'%(i,j,dvarij.varname, dvarij.x))

        distvalue2 += setvalue2[i]
        numpositive += setvalue2[i] > TOL
        # log.joint('Set %d sum-of-squares %g; binary owner value %g\n'%(i, setvalue2[i], owneri.x))

        # new: objective
        for j in range(len(x_i)):
            # ind_j = j + i*12
            if x_i[j] >= 0: 
                # log.joint('( x%d - %g ) ^ 2 + '%(ind_j, x_i[j]))
                # log.joint('[ x%d ^ 2 - %g x%d + %g ] + '%(ind_j, (2*x_i[j]), ind_j, (x_i[j])**2))
                # log.joint('- %g x%d '%((2*x_i[j]), ind_j))
                log.joint('- %g %s '%((2*x_i[j]), distvarset[i][j].varname))
            if x_i[j] < 0: 
                # log.joint('( x%d + %g ) ^ 2 + '%(ind_j, abs(x_i[j])))
                # log.joint('[ x%d ^ 2 + %g x%d + %g ] + '%(ind_j, (2*abs(x_i[j])), ind_j, (x_i[j])**2))
                # log.joint('+ %g x%d '%((2*abs(x_i[j])), ind_j))
                log.joint('+ %g %s '%((2*abs(x_i[j])), distvarset[i][j].varname))
        # log.joint('\n')
        #simplebreak()
    # log.joint('\n') # new

    log.joint('+ [ 2 %s ^ 2 '%(distvarset[0][0].varname))

    for j in range(1, len(x_i)):
            # ind_j = j + i*12
            # log.joint('+ 2 x%d ^ 2 '%(ind_j))
            log.joint('+ 2 %s ^ 2 '%(distvarset[0][j].varname))

    for i in range(1, distcount):
        x_i = []
        # varij_name = []

        for j in range(distsize[i]):
            dvarij = distvarset[i][j]
            xvalue = dvarij.x

            # varij_name.append(dvarij.varname)

            # if xvalue too small just append 0
            if abs(xvalue) <= TOL: x_i.append(0)
            if abs(xvalue) > TOL:
                # new
                x_i.append(xvalue)
        # new
        for j in range(len(x_i)):
            # ind_j = j + i*12
            # log.joint('+ 2 x%d ^ 2 '%(ind_j))
            log.joint('+ 2 %s ^ 2 '%(distvarset[i][j].varname))
    # log.joint('\n')

    log.joint('] / 2 \n \n')
    log.joint('Subject To \n  c0: ')

    # new: constraints 

    for i in range(distcount):
        x_i = []
        # varij_name = []

        for j in range(distsize[i]):
            dvarij = distvarset[i][j]
            xvalue = dvarij.x

            # if xvalue too small just append 0
            if abs(xvalue) <= TOL: x_i.append(0)
            if abs(xvalue) > TOL:
                # new
                x_i.append(xvalue)
        # new
        for j in range(len(x_i)):
            # ind_j = j + i*12
            if x_i[j] >= 0: 
                # log.joint('+ %g * ( x%d - %g ) '%(x_i[j], ind_j, x_i[j]))
                # log.joint('+ [ %g x%d - %g ] '%(x_i[j], ind_j, (x_i[j])**2))
                # log.joint('+ %g x%d '%(x_i[j], ind_j))
                log.joint('+ %g %s '%(x_i[j], distvarset[i][j].varname))
            if x_i[j] < 0: 
                # log.joint('- %g * ( x%d + %g ) '%(abs(x_i[j]), ind_j, abs(x_i[j])))
                # log.joint('- [ %g x%d + %g ] '%(abs(x_i[j]), ind_j, (x_i[j])**2))
                log.joint('- %g %s '%(abs(x_i[j]), distvarset[i][j].varname))
    # log.joint('\n')

    sum_of_constants = 0

    for i in range(distcount):
        x_i = []
        for j in range(distsize[i]):
            dvarij = distvarset[i][j]
            xvalue = dvarij.x

            # if xvalue too small just append 0
            if abs(xvalue) <= TOL: x_i.append(0)
            if abs(xvalue) > TOL:
                # new
                x_i.append(xvalue)
        # new
        for j in range(len(x_i)):
            ind_j = j + i*12
            sum_of_constants += - (x_i[j])**2
    # log.joint('sum of constants is: %g '%(sum_of_constants))
    # log.joint('\n')
    log.joint('= %g \n'%(abs(sum_of_constants)))

    for i in range(distcount):
        for j in range(distsize[i]): 
            ind_j = j + i*12
            # log.joint('  c%d: -12 y%d - x%d <= 0 \n'%((2*ind_j+1), i, ind_j))
            log.joint('  c%d: -12 b%d - %s <= 0 \n'%((2*ind_j+1), (i+102), distvarset[i][j].varname))

            # log.joint('  c%d: x%d - 12 y%d <= 0 \n'%((2*ind_j+2), ind_j, i))
            log.joint('  c%d: %s - 12 b%d <= 0 \n'%((2*ind_j+2), distvarset[i][j].varname, (i+102)))

    log.joint('  c2401:')
    for i in range(distcount-1):
        log.joint(' b%d +'%(i+102))
    log.joint(' b%d <= 75 '%(distcount-1+102))

    log.joint('\n \nBounds \n')
    for i in range(distcount): 
        for j in range(distsize[i]): 
            log.joint('  %s Free \n'%(distvarset[i][j].varname))

    log.joint('\nBinary \n  ')

    for i in range(distcount):
        log.joint('b%d '%(i+102))
    log.joint('\n \nEnd\n \n')

    # if header:
        # log.joint('%s '%(header))
        # log.joint('')
    # log.joint('Positives: %d; sum total %g\n'%(numpositive, distvalue2))

    ind = np.argsort(-setvalue2)
    ordered = setvalue2[ind]
    
    '''
    print(ordered)
    print(np.sum(ordered[75:]))
    '''


    for k in range(distcount):
        i = ind[k]
        alldata['keyorder'][k] = i   
    

    sumsmallest2 = 0
    cardthresh = distcount - alldata['notmaxcard']

    rhssum = 0
    for k in range(cardthresh, distcount):
        i = ind[k]
        owneri = distvar_owner[i]

        sumsmallest2 += ordered[k]

        rhssum += owneri.x*(ordered[ind[cardthresh-1]] - ordered[k])

        # if header:
        #     log.joint('')
        # log.joint('Ordered set %d is %d (%s, %g) at %g\n'%(k,i, owneri.varname, owneri.x, ordered[k]))
        

    # log.joint('rhssum: %g\n'%(rhssum))
    # if header:
        # log.joint('%s '%(header))
    # log.joint('Sum of %d smallest: %g\n'%(distcount - cardthresh, sumsmallest2))
    # if header:
        # log.joint('%s '%(header))
    # log.joint('So overall lower bound: %g\n'%(distvalue2 + sumsmallest2))


    return code, distvalue2 + sumsmallest2
        

def laresetLPtoParametric(alldata, ss):
    log = alldata['log']

    log.joint('Now building parametric RHS LP using ss %g.\n'%(ss))

    distcount = alldata['distcount']
    distsize = alldata['distsize']
    distvarset = alldata['distvarset']
    
    model = alldata['model']

    lambdavar = alldata['lambdavar'] = model.addVar(obj = 1.0, lb = 0, name = "lambda")

    #set up new objective
    lambdaobj = lambdavar
        
    model.setObjective(lambdaobj,GRB.MINIMIZE)

    log.joint('Reset objective to lambda.\n')



    #simplebreak()

    LOUD = False

    #now create parametric constraint

    distvar_owner = alldata['distvar_owner'] 
    distvar_partner = alldata['distvar_partner']
    distconstr = alldata['distconstr']

    cheatcoeff = alldata['cheatcoeff'] # not needed

    USELAZY = False

    if USELAZY == False:
        #boundvar = {}
        ubval = cheatcoeff**.5

    lhs = 0        
    for i in range(distcount):
        for j in range(distsize[i]):
            dvarij = distvarset[i][j]
            lhs += dvarij*dvarij
    
    model.addConstr(lhs <= ss*lambdavar, name = 'lambda_up')

    model.reset()

    if USELAZY:
        model.write('lambdalazy.lp')
    else:
        model.write('lambdanotlazy.lp')
    simplebreak()



def lapushdownParametric(alldata, keyi):
    log = alldata['log']

    log.joint('Now pushing down parametric RHS LP using keyi %d.\n'%(keyi))

    distcount = alldata['distcount']
    distsize = alldata['distsize']
    distvarset = alldata['distvarset']
    
    model = alldata['model']

    LOUD = False

    distvar_owner = alldata['distvar_owner'] 
    distvar_partner = alldata['distvar_partner']
    distconstr = alldata['distconstr']

    USELAZY = False


    ub = {}
    lb = {}
    for j in range(distsize[keyi]):
            dvarij = distvarset[keyi][j]
            ub[j] = dvarij.ub
            lb[j] = dvarij.lb
            dvarij.lb = dvarij.ub = 0
    
    model.reset()

    if USELAZY:
        model.write('lambdalazypush'+str(keyi)+'.lp')
    else:
        model.write('lambdanotlazypush'+str(keyi)+'.lp')


    code, pushvalue = lasolveLP(alldata,'push'+str(keyi))

    #restore

    for j in range(distsize[keyi]):
            dvarij = distvarset[keyi][j]
            dvarij.ub = ub[j]
            dvarij.lb = lb[j]

    log.joint('key%d done with code %d val %g\n'%(keyi, code, pushvalue))
    #simplebreak()

    return code, pushvalue

    
