# -*- coding: utf-8 -*-
import csv
import math
import time
import datetime
import traceback
import os
import sys


from log import danoLogger
from utils import breakexit, returnbreakexit
from versioner import *
from lalp import golp
from ladoit import ladoit, lascan
from lalift import lalift
from lacut import *

def laread_solution(alldata, filename):

    log = alldata['log']

    log.joint("reading solution file " + filename + "\n")

    try:
        f = open(filename, "r")
        lines = f.readlines()
        f.close()
    except:
        log.stateandquit("cannot open file " + filename + "\n")

    vectordictionary = {}
    model = alldata['model']
    seen = {}    
    for var in model.getVars():
        seen[var.varname] = 1


    linenum = 0

    while linenum < len(lines):
     thisline = lines[linenum].split()

     if len(thisline) > 0 and thisline[0][0] != '#':
         if thisline[0] != 'END':         
            varname = thisline[0]
            varvalue = float(thisline[2])
            #print(varname, varvalue)
            if varname not in seen:
                 print(varname, varvalue)                 
                 breakexit('hey')
            vectordictionary[varname] = varvalue
         elif thisline[0] == 'END':
             break
     linenum += 1
    
    return vectordictionary
    
def read_config(log, filename):

    log.joint("reading config file " + filename + "\n")

    try:
        f = open(filename, "r")
        lines = f.readlines()
        f.close()
    except:
        log.stateandquit("cannot open file " + filename + "\n")

    alldata = {}
    lpfile = 'NONE'
    solutionfile = None
    linenum = 0
    ignorevarlindeg = False

    MAXCARD = 75
    UB = 125
    VERSION = 'SumOfSquares'
    ACTION = 'Standard'
    

    while linenum < len(lines):
     thisline = lines[linenum].split()

     if len(thisline) > 0 and thisline[0][0] != '#':
       word = thisline[0]
       if thisline[0] == 'LPFILE':
         lpfile = thisline[1]
       elif thisline[0] == 'MAXCARD':
         MAXCARD = int(thisline[1])
       elif thisline[0] == 'UB':
         UB = float(thisline[1])
       elif thisline[0] == 'VERSION':
         VERSION = thisline[1]
       elif thisline[0] == 'ACTION':
         ACTION = thisline[1:3]
       elif thisline[0] == 'SOLUTIONFILE':
         solutionfile = thisline[1]
       elif thisline[0] == 'IGNOREVARLINDEG':
         ignorevarlindeg = True
       elif thisline[0] == 'END':
         break
       else:
         log.stateandquit("illegal input " + thisline[0] + "\n")

     linenum += 1

    for x in [('LPFILE',lpfile), ('MAXCARD',MAXCARD), ('UB', UB), ('VERSION',VERSION), ('ACTION', ACTION), ('SOLUTIONFILE', solutionfile), ('IGNOREVARLINDEG', ignorevarlindeg)]:
      if x[1] == 'NONE':
        log.stateandquit(' no ' + x[0] + ' input'+'\n')
      alldata[x[0]] = x[1]
      log.joint(x[0] + ' ' + str(x[1])+ '\n')

    alldata['log'] = log
 
    return alldata

if __name__ == "__main__":

    logfile = 'la.log'

    if len(sys.argv) > 3:
        sys.exit('usage: la0.py configfile [logfile]')
    if len(sys.argv) == 1:
        configfile = 'la.conf'
    if len(sys.argv) == 2:
        configfile = sys.argv[1]
    if len(sys.argv) == 3:
        logfile = sys.argv[2]
    log = danoLogger(logfile)
    stateversion(log)

    alldata = read_config(log, configfile)
    alldata['log'] = log


    if alldata['LPFILE'] != 'NONE':
        golp(alldata) # creates the initial LP (by reading the .lp file that has the LP)

        lascan(alldata) # checks all the variables and creates some of the variable dictionaries

        if alldata['SOLUTIONFILE'] != 'None':
            VECDIC= laread_solution(alldata, alldata['SOLUTIONFILE'])
            alldata['vectordictionary'] = VECDIC

        if alldata['ACTION'][0] == 'Lift_and_Project':
            laCuttingPlane(alldata)
            #lalift(alldata, liftingvariablename = alldata['ACTION'][1], vectordictionary = VECDIC)
        else:
            ladoit(alldata)

    log.closelog()    
