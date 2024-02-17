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

    # vectordictionary will be updated in laCuttingPlaneEachRound each time it is called. 

    log.joint('Starting cutting plane algorithm.\n')



    distcount = alldata['distcount']
    distsize = alldata['distsize']
    distvarset = alldata['distvarset']
    distvar_owner = alldata['distvar_owner'] 
    distvar_partner = alldata['distvar_partner']
    distconstr = alldata['distconstr']
    varSet = alldata['varSet']
    distset_byname = alldata['distset_byname']
    
    # Change variable type for owner variables
    log.joint('Switching to continuous relaxation of original model.\n')
    for i in range(distcount):
        owner = distvar_owner[i]
        owner.vtype = GRB.CONTINUOUS

    # 0. First, we solve the initial relaxation (OR READ from file)        

    # if solution not available: solve continuous relaxation
    if 'vectordictionary' not in alldata.keys():
        log.joint('Solving continuous relaxation of original model\n')
        model.optimize()
        alldata['vectordictionary'] = vectordictionary = {}
        for v in model.getVars():
            vectordictionary[v.varname] = v.x
        log.joint('Solved, and mapped solution to vectordictionary.\n')
        #sys.exit('pooooooo')

    vectordictionary = alldata['vectordictionary']        
   
    # Then, iterate T times (at the beginning, T should be small, like 1
    # T = 5 for now, we might change later. 

    T = 5 # Change this later to be a parameter rather than hard-coded. 

    for i in range(T): 
        log.joint("\n*****This is the start of round %d\n"%(i+1))

        # Record the index of the current round. 
        alldata['round_ind'] = round_ind = i

        # Initialize scores in each round. 
        # alldata['scores'] = scores = np.zeros(distcount) 

        # Sort the binaries based on how close the values are to 1/2. 
        # alldata['sorted_order'] = sorted_order = laSelectMostFractional(alldata)

        # Add the cuts and solve them. 
        laCuttingPlaneEachRound(alldata)

def laCuttingPlaneEachRound(alldata):

    # 1. Pick the K "most fractional" binary variables.  Most fractional means closest to 1/2.  
    # And K should be small at the beginning, like K = 1.  
    # Later K = 10 (maybe -- need to experiment with this).

    K = 25 # Change this later to be a parameter rather than hard-coded. 

    log = alldata['log']
    model = alldata['model']


    distcount = alldata['distcount']
    if K > distcount: K = distcount
    distsize = alldata['distsize']
    distvarset = alldata['distvarset']
    distvar_owner = alldata['distvar_owner'] 
    distvar_partner = alldata['distvar_partner']
    distconstr = alldata['distconstr']
    varSet = alldata['varSet']
    distset_byname = alldata['distset_byname']

    # Get sorted_order and scores to pick the K most fractional binary variables. 

    # sorted_order = alldata['sorted_order']
    # scores = alldata['scores']

    # vectordictionary records the solution to the model in the last round, 
    # and will be updated after this round. 

    vectordictionary = alldata['vectordictionary']

    # 2. For each of the identified K binary variables, we call lalift using that variable as the first argument, 
    # but always using the same solution to the relaxation in the vectordictionary. 

    # Import index for current round to name the constraints. 
    round_ind = alldata['round_ind']

    list_of_expr = []
    list_of_cutrhs = []
    list_of_diff = []
    list_of_constr_name = []

    for i in range(distcount): 
        # curr_ind = sorted_order[i]
        curr_ind = i

        # Get variable name corresponding to current index. 
        currIndVariableName = distvar_owner[curr_ind].VarName

        # The i^th most fractional value is... 
        # log.joint('The %d most fractional value is %f by %s with difference %f\n'%(i+1, vectordictionary[currIndVariableName], currIndVariableName, scores[curr_ind]))
        
        # Get output of lalift. 
        expr, cutrhs, diff = lalift(alldata, liftingvariablename = currIndVariableName, vectordictionary = vectordictionary)

        list_of_expr.append(expr)
        list_of_cutrhs.append(cutrhs)

        # This is negative for the sake of argsort (argsort sorts in ascending order). 
        list_of_diff.append(-diff)

        # 3. After each call to lalift, store the computed cut (if a cut was successfuly computed) in some data structure.
        # This was done in lalift, which returns a linear expression and a value for the rhs. 
        
        # Add constraint from the cut outputted by lalift. 
        # Add a name for the constraint 
        constr_name = 'cut_round_' + str(round_ind) + '_var_' + currIndVariableName
        list_of_constr_name.append(constr_name)

        # model.addConstr(expr >= cutrhs, name=constr_name)

    # 4. Once all calls are completed, add the K most violated cuts to the relaxation.
    sorted_order_violation = np.argsort(list_of_diff)
    diff_threshold = 1e-4
    print("Diff_threshold: %.3e"%diff_threshold)
    for i in range(K): 
        i_most_violated_ind = sorted_order_violation[i]

        diff_i = -list_of_diff[i_most_violated_ind] 


        if diff_i > diff_threshold: 
            violatedIndVariableName = distvar_owner[i_most_violated_ind].VarName

            expr_i = list_of_expr[i_most_violated_ind]
            cutrhs_i = list_of_cutrhs[i_most_violated_ind]
            constr_name_i = list_of_constr_name[i_most_violated_ind]

            model.addConstr(expr_i >= cutrhs_i, name=constr_name_i)

            log.joint('The %d most violated cut is obtained by lifting %s with violation %f\n'%(i+1, violatedIndVariableName, diff_i))

    # Write model with K cuts. 
    model.write('Dcut'+str(alldata['round_ind'])+'.lp')

    # 5. Solve the updated relaxation.  Go to 1.

    log.joint('Now solving updated model\n')

    model.optimize()

    # Print optimal objective to log. 
    log.joint('Optimal objective = %g\n' % model.objVal)

    # newsolution records the optimal solution. 
    newsolution = {}

    for v in model.getVars():
        if v.varname in vectordictionary: 
            newsolution[v.varname] = v.x

    alldata['vectordictionary'] = newsolution


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
        # print(i, varname, value)

    sorted_order = np.argsort(scores)

    return sorted_order
