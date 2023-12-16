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



def lalift(alldata, distset_index):
    log = alldata['log']

    global MYDEBUG
    MYDEBUG = False

    # For now, in each iteration we will formulate and solve a new problem
    # The problem will consider a given binary variable: z_j
    # Given (x*, y*, z*) the problem will compute the distance from (x*,y*,z*) to the
    # projection, in (x,y,z) space, of the z_j - convexified lifted formulation of the problem,
    # That is to say, the convex hull of feasible solutions where z_j = 0 or 1.

    # The formulation in the lifted space will have four parts:
    # (a) A copy of the formulation multiplied by z_j and linearized; corresponding to z_j = 1.
    # (b) A copy of the formulation multiplied by 1-z_j and linearized; corresponding to z_j = 0.
    # (c) A constraint saying that (x,y,z) = (x1, y1, z1) + (x0, y0, z0), where the "1" or "0".
    #     refers to the vector in (a) or (b), respectively.
    # (d) Have to handle individual variable bounds.
    # (e) Constraints and variables used to model the distance between (x,y,z) and (x*,y*,z*).

    

    # For (a) and (b) we will use the same code, where the "1" and "0" case are distinguished
    # using a flag.  We need logic so as to

    #

    # (1) in the zero case, we need to 
    #   Keep track of the variables that are in set S_j (lhs coefficient becomes zero).
    #   Appropriately handle rhs (e.g., >= 5 becomes  5 z_j >= 5).
    #   Appropriately handle cardinality constraint.

    # (2) In the one case, 
    #   Appropriately handle rhs (e.g., >= 5 becomes -5 z_j >= 0).
    #   Appropriately handle cardinality constraint.

    #

    # To help with this, varclass is used to keep track of the nature of variables
    # 'indicator' (binary), 'distance_sum' (y variable) and 'groundset' (x variable).


    distcount = alldata['distcount']
    distsize = alldata['distsize']
    distvarset = alldata['distvarset']
    distvar_owner = alldata['distvar_owner'] 
    distvar_partner = alldata['distvar_partner']
    distconstr = alldata['distconstr']
    

    liftingvariable = distvar_owner[distset_index]
    log.joint('Lifting using set %d variable %s\n' %(distset_index, liftingvariable.varname))

    model = alldata['model']

    Dmodel = Model('disjunctmodel')

    # First we will create the variables used in the disjunction: "_0" and "_1".
    #
    # We run through all variables in the original model and we create the _0 and _1 lifted variables,
    #   as well as a copy of the variable itself (which we do first first)

    for var in model.getVars():
        Dmodel.addVar(obj = 0.0, lb = var.lb, ub = var.ub, name = var.varname)
    for DisjunctionCase in range(2):
        for var in model.getVars():
            newname = var.varname + '_' + str(DisjunctionCase)
            #if newname[0] == 'b': print(newname)
            if var.varname == liftingvariable.varname: print(newname)            
            Dmodel.addVar(obj = 0.0, lb = var.lb, ub = var.ub, name = newname)

    # there is some waste in this, in particular in terms of the lifting variable, but it's OK for now

    #
    Dmodel.update()
     
    constrs = model.getConstrs()
    for DisjunctionCase in range(2):
        Dliftingvariable = Dmodel.getVarByName(liftingvariable.Varname+'_'+str(DisjunctionCase))
        for constr in constrs:
            expr = LinExpr() #the new constraint
            
            lhs = model.getRow(constr)
            exprSize = lhs.size()
            for k in range(exprSize): 
                actualVar = lhs.getVar(k)
                #print(actualVar.Varname)
                liftedvar = Dmodel.getVarByName(actualVar.Varname+'_'+str(DisjunctionCase))
                
                coefficient = lhs.getCoeff(k)

                if actualVar.Varname == liftingvariable.Varname and DisjunctionCase == 0:
                    coefficient = 0
            
                expr += coefficient*liftedvar

            rhsval = constr.RHS
            expr -= rhsval*Dliftingvariable

            print(constr.ConstrName)
            print(expr)
            if constr.Sense == '=': 
                Dmodel.addConstr(expr == 0, name = constr.ConstrName + '_' + str(DisjunctionCase))
            elif constr.Sense == '>': 
                Dmodel.addConstr(expr >= 0, name = constr.ConstrName + '_' + str(DisjunctionCase))
            elif constr.Sense == '<': 
                Dmodel.addConstr(expr <= 0, name = constr.ConstrName + '_' + str(DisjunctionCase))

    Dmodel.update()

    Dmodel.write('D.lp')


    breakexit('lalifted')
