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


def golp(alldata):
    log = alldata['log']

    global MYDEBUG
    MYDEBUG = False

    tstart = time.time()

    getlp(alldata)
    
    n = alldata['n']

 
def getlp(alldata):
    log = alldata['log']

    log.joint('Handling %s.\n' % alldata['LPFILE'])

    model = read(alldata['LPFILE'])

    alldata['model'] = model

    myvars = {}
    revvars = {}
    revvarsbyname = {}


    reversepairs = {}
    pairs = {}
    qobjpairs = {}    

    log.joint('Scanning model.\n')    

    log.joint("Number of variables = " + str(model.NumVars)+'.\n')

    alldata['n'] = model.NumVars
    alldata['internal_varvalue_holder'] = np.zeros(alldata['n'])
    
    numLunbounded = numUunbounded = numBinary = numInteger = 0
    varcount = 0

    LOUD = False

    alldata['varlindeg'] = varlindeg = {}

    for var in model.getVars():
        varcount += 1
        myvars[varcount] = var  #base1, not zero
        varlindeg[var.varname] = 0
        revvars[var] = varcount
        revvarsbyname[var.varname] = varcount
        Lunbounded = (var.lb <= -GRB.INFINITY)
        numLunbounded += Lunbounded
        Uunbounded = (var.ub >= GRB.INFINITY)
        numUunbounded += Uunbounded
        numBinary += (var.vtype == GRB.BINARY)
        numInteger += (var.vtype == GRB.INTEGER)
        if LOUD:
            log.joint(var.varname +' type ' + var.vtype + ' lb = ' + str(var.lb) + ' ' + str(Lunbounded) +' ub = ' +  str(var.ub) + ' '  + str(Uunbounded)+'\n')
    log.joint(' Lunbounded: '+ str(numLunbounded) +', Uunbounded: '+ str(numUunbounded) + ', Binary: '+ str(numBinary) + ', Integer: '+ str(numInteger)+'.\n')
    
    #simplebreak()
    
    pairscount = 0
    obj = model.getObjective()
    objisquad = False

    internal_quadsize = 0
    internal_quad = {} #both for now

    internal_pair_varindices1 = {}

    if LOUD:
        log.joint('obj = '+ str(obj))
        
    if isinstance(obj, QuadExpr):
        log.joint('\n->quadratic objective!\n')
        objisquad = True
        lin_obj = obj.getLinExpr()
        q_obj = obj - lin_obj
        if LOUD:
            print('q_obj ='+ str(q_obj))

        internal_quadsize = q_obj.size()
        for j in range(q_obj.size()):
            v1 = q_obj.getVar1(j)
            v2 = q_obj.getVar2(j)
            val = q_obj.getCoeff(j)
            if ((v1,v2) not in reversepairs) and ((v2,v1) not in reversepairs) :
                pairscount += 1                
                reversepairs[(v1,v2)] = pairscount
                pairs[pairscount] = (v1,v2)
                qobjpairs[pairscount] = q_obj.getCoeff(j)

                k1 = revvars[v1]  #recall that revvars is the variable index, counted from zero
                k2 = revvars[v2]
                internal_pair_varindices1[pairscount] = (k1, k2)
                internal_quad[j] = (k1,k2, v1,v2, pairscount, q_obj.getCoeff(j))  #a bit redundant
                
                #print(k1, k2)
                #breakexit('')
            if LOUD:
                log.joint( str(pairscount) + '. '+ str(q_obj.getCoeff(j)) + ' ' + str(q_obj.getVar1(j))+ ' ' + str(q_obj.getVar2(j)) + '\n')
                
    else:
        log.joint('\nLinear objective.\n')        
        lin_obj = obj
        if LOUD:
            for j in range(lin_obj.size()):
                log.joint(' '+ str(lin_obj.getCoeff(j))+ ' ' + str(lin_obj.getVar(j))+'\n')
    log.joint('Quad count, after obj: %d.\n'  %(pairscount))
    alldata['internal_quadsize'] = internal_quadsize
    alldata['internal_quad'] = internal_quad
    
    log.joint('Constant '+ str(lin_obj.getConstant()) + '.\n')
    alldata['internal_constant'] = lin_obj.getConstant()


    internal_linsize = lin_obj.size()
    internal_lin = {}
    for j in range(lin_obj.size()):
        k = revvars[lin_obj.getVar(j)]
        internal_lin[j] = (k,lin_obj.getVar(j),lin_obj.getCoeff(j)) #not sure all are needed

    alldata['internal_linsize'] = internal_linsize
    alldata['internal_lin'] = internal_lin
    
    #log.joint("\nQ constraints = " + str(model.NumQConstrs)+'\n')    
    Qconstrs = model.getQConstrs()
    numnonzerorhsQconstrs = 0

    LOUD = False
    for Qconstr in Qconstrs:
        therow = model.getQCRow(Qconstr)
        if LOUD:
            log.joint(Qconstr.QCName+ ' ' + Qconstr.QCSense + ' ' + str(Qconstr.QCRHS) + '\n')
            #print('therow', therow)

        linQc  = therow.getLinExpr()
        rhs = Qconstr.qcrhs

        numnonzerorhsQconstrs += (rhs != 0)
        quadQc = therow - linQc

        if LOUD:
            log.joint('quadQc size '+ str(quadQc.size()) + '.\n')

            simplebreak()

        for j in range(quadQc.size()):
            v1 = quadQc.getVar1(j)
            v2 = quadQc.getVar2(j)
            if ((v1,v2) not in reversepairs) and ((v2,v1) not in reversepairs) :
                pairscount += 1                
                reversepairs[(v1,v2)] = pairscount
                pairs[pairscount] = (v1,v2)

                k1 = revvars[v1]  #recall that revvars is the variable index, counted from zero
                k2 = revvars[v2]
                internal_pair_varindices1[pairscount] = (k1, k2)

            if LOUD:
                log.joint( str(pairscount) + '. '+ str(quadQc.getCoeff(j))+ ' ' + str(quadQc.getVar1(j)) + ' '  + str(quadQc.getVar2(j))+'\n')

    log.joint('Number of quadratic constraints: %d; nonhomogeneous: %d.\n' %(len(Qconstrs), numnonzerorhsQconstrs))                

    alldata['internal_pair_varindices1'] = internal_pair_varindices1
    alldata['internal_pairvalue_holder'] = np.zeros(pairscount)
                
    log.joint('Quadcount, after constraints: %d.\n' %(pairscount))
    #simplebreak()
    #log.joint("\nlinear constraints = " + str(model.NumConstrs) + '\n')
    constrs = model.getConstrs()
    numnonzerorhsconstrs = 0

    for constr in constrs:
        numnonzerorhsconstrs += (constr.rhs != 0)
        lhs = model.getRow(constr)
        for j in range(lhs.size()):
            varlindeg[lhs.getVar(j).varname] += 1

        if LOUD:
            log.joint(constr.ConstrName + ' ' + constr.Sense + ' ' + str(constr.RHS) + '\n')            
            log.joint('constr size '+ str(lhs.size()) + '\n')
            for j in range(lhs.size()):
                print(lhs.getCoeff(j),lhs.getVar(j))
            log.joint('constant '+ str(lhs.getConstant()) + '\n')

    log.joint('Number of linear constraints: %d; nonhomogeneous: %d.\n' %(len(constrs), numnonzerorhsconstrs))

    constrs = model.getConstrs()

    numOfLinConstrs = len(constrs)
    alldata['numOfLinConstrs'] = numOfLinConstrs

    LOUD = False 
    isallbinary = np.zeros(len(constrs),dtype=int)
    ishomogeneous = np.zeros(len(constrs),dtype=int)    
    i = 0
    for constr in constrs:
        lhs = model.getRow(constr)
        if LOUD:
            log.joint(constr.ConstrName + ' ' + constr.Sense + ' ' + str(constr.RHS) + '\n')

        isallbinary[i] = 1
        for j in range(lhs.size()):
            #print(i,j, lhs.getCoeff(j),lhs.getVar(j).varname, lhs.getVar(j).vtype)
            if lhs.getVar(j).vtype == 'C':
                isallbinary[i] = 0
                #print(' --> not all binary!')
                break


        ishomogeneous[i] = (constr.RHS == 0)
        if ishomogeneous[i] == False and isallbinary[i] == 0:
            log.joint(constr.Constrname + ' nonhomogeneous, RHS '+ str(constr.RHS) + '\n')
        i += 1

    log.joint('Number of linear constraints = %d\n'%(len(constrs)))
    log.joint('Number of all-binary linear constraints = %d\n'%(np.sum(isallbinary)))
    log.joint('Number of homogeneous linear constraints = %d\n'%(np.sum(ishomogeneous)))

    alldata['isallbinary'] = isallbinary
    alldata['ishomogeneous'] = ishomogeneous
    breakexit('bar')

    if LOUD:
        log.joint('Bilinear/quadratic terms:')
        for k in range(1,1+pairscount):
            log.joint('%d. %s %s -> %d %d\n' %(k, pairs[k][0].varname, pairs[k][1].varname, revvars[pairs[k][0]], revvars[pairs[k][1]]) + '\n')


    alldata['pairscount'] = pairscount
    alldata['pairs'] = pairs
    alldata['myvars'] = myvars
    alldata['revvarsbyname'] = revvarsbyname
    alldata['revvars'] = revvars    


    n = alldata['n']
    #simplebreak()





