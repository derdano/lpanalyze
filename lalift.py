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



def lalift(alldata):
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
