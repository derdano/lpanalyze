import sys
import os.path

def simplebreak():
    stuff = input('break > ')
    if stuff == 'x' or stuff == 'q':
        sys.exit("bye")

def breakexit(foo):
    stuff = input("("+foo+") break> ")
    if stuff == 'x' or stuff == 'q':
        sys.exit("bye")

def askfordouble(foo):
    stuff = input("("+foo+") input value> ")
    value = float(stuff)
    return value

def checkstop(log,name):
    if os.path.isfile(name) :
        breakexit("stop?")

def myprintfile(log, mfilename, casefilelines):
    try:
        f = open(mfilename, "w")
    except:
        log.stateandquit("cannot open file " + mfilename + "\n")

    for line in casefilelines.values():
        f.write(line)

    f.close()

    return 1

def myreadfile(log, filename):
    code = 0
    try:
        f = open(filename, "r")
        lines = f.readlines()
        f.close()
    except:
        log.joint("cannot open file " + filename + "\n")
        code = 1

    return code, lines
