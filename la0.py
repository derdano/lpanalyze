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
from ladoit import ladoit

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
    linenum = 0

    while linenum < len(lines):
     thisline = lines[linenum].split()

     if len(thisline) > 0 and thisline[0][0] != '#':
       word = thisline[0]
       if thisline[0] == 'LPFILE':
         lpfile = thisline[1]
       elif thisline[0] == 'END':
         break
       else:
         log.stateandquit("illegal input " + thisline[0] + "\n")

     linenum += 1

    for x in [('LPFILE',lpfile)]:
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
        golp(alldata)

        ladoit(alldata)

    log.closelog()    
