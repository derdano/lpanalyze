import gurobipy as gp
import numpy as np
from gurobipy import *
from numpy import linalg as LA
import sys
import struct
from log import danoLogger
from myutils import breakexit, askfordouble
import time

def check_equality():

	# file1 is the txt file that contains x^*
	file1 = open('hyperperspcheck-little_brother_test.sol', 'r')

	array1 = []
	for line in file1: array1.append([x.strip() for x in line.split()])

	p_star = []
	p_star_varname = []

	for i in range(931, len(array1)): 
		item = float(array1[i][1])
		varname = str(array1[i][0])
		if varname[0] == 'p': 
			p_star.append(item)
			p_star_varname.append(varname)

	print("The length of p is: " + str(len(p_star)))

	for i in range(len(p_star)): 
		print(p_star_varname[i] + ": " + str(p_star[i]))
check_equality()