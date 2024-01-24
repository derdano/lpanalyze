# -*- coding: utf-8 -*-
import gurobipy as gp
import numpy as np
from gurobipy import *
from numpy import linalg as LA
import sys
import struct
from log import danoLogger
from myutils import breakexit, askfordouble
from lalift import *
import time



def laCuttingPlane(alldata):

    log = alldata['log']
    model = alldata['model']

    # 0. First, we solve the initial relaxation (OR READ from file)

    distcount = alldata['distcount']
    distsize = alldata['distsize']
    distvarset = alldata['distvarset']
    distvar_owner = alldata['distvar_owner'] 
    distvar_partner = alldata['distvar_partner']
    distconstr = alldata['distconstr']
    varSet = alldata['varSet']
    distset_byname = alldata['distset_byname']
    
    # Change variable type for owner variables

    for i in range(distcount):
        owner = distvar_owner[i]
        owner.vtype = GRB.CONTINUOUS

    alldata['scores'] = scores = np.zeros(distcount)    
    # Then, iterate T times (at the beginning, T should be small, like 1

    # 1. Pick the K "most fractional" binary variables.  Most fractional means closest to 1/2.  And K should be small at the beginning, like K = 1.  Later K = 10 (maybe -- need to experiment with this).

    # 2. For each of the identified K binary variables, we call lalift using that variable as the first argument, but always using the same solution to the relaxation in the vectordictionary

    sorted_order = laSelectMostFractional(alldata)

    worst_ind = sorted_order[0]

    mostFractionalVariableName = distvar_owner[worst_ind].VarName

    log.joint('Most fractional value is %f by %s\n'%(scores[worst_ind], mostFractionalVariableName))

    lalift(alldata, liftingvariablename = mostFractionalVariableName, vectordictionary = alldata['vectordictionary'])

    # 3. After each call to lalift, store the computed cut (if a cut was successfuly computed) in some data structure.

    # 4. Once the K calls are completed, add the cuts to the relaxation.

    # 5. Solve the updated relaxation.  Go to 1.

    log.joint('Now solving updated model\n')

    model.optimize()


def laSelectMostFractional(alldata):

    log = alldata['log']
    vectordictionary = alldata['vectordictionary']

    distcount = alldata['distcount']
    distsize = alldata['distsize']
    distvarset = alldata['distvarset']
    distvar_owner = alldata['distvar_owner'] 
    distvar_partner = alldata['distvar_partner']
    distconstr = alldata['distconstr']
    varSet = alldata['varSet']
    distset_byname = alldata['distset_byname']
    scores = alldata['scores']

    for i in range(distcount):
        owner = distvar_owner[i]
        varname = owner.VarName
        value = abs(vectordictionary[varname] - 0.5)
        scores[i] = value
        print(i, varname, value)

    sorted_order = np.argsort(scores)
    #print('scores',scores)
    #print('sorted_order', sorted_order)

    #print(sorted_order)

    return sorted_order
