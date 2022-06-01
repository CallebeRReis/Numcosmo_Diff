import numpy as np
import matplotlib.pyplot as plt
# from Derivative_Numcosmo.derivative import rsCoef
import derivative as d
import diff_table
import copy 
def fowardDiff(f,x,h):
    return (f(x+h)-f(x))/h
def func(x):
    return np.exp((x**2))
def dFunc(x):
    return 2*x*np.exp(x**2.0)

def parameter(index):
    return np.power(2.0,index)
 


size = 5


a = d.Derivative(func, 0.1, size)
b = diff_table.Diff_Table(a,0.1,parameter,1.0)
b.build_Table_Backward()
b.printTable()

