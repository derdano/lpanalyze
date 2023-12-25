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

    # 0. First, we solve the initial relaxation

    # Then, iterate T times (at the beginning, T should be small, like 1

    # 1. Pick the K "most fractional" binary variables.  Most fractional means closest to 1/2.  And K should be small at the beginning, like K = 1.  Later K = 10 (maybe -- need to experiment with this).

    # 2. For each of the identified K binary variables, we call lalift using that variable as the first argument, but always using the same solution to the relaxation in the vectordictionary

    # 3. After each call to lalift, store the computed cut (if a cut was successfuly computed) in some data structure.

    # 4. Once the K calls are completed, add the cuts to the relaxation.

    # 5. Solve the updated relaxation.  Go to 1.
