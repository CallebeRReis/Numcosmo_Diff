import numpy as np 
import matplotlib.pyplot as plt
from scipy.fftpack import diff 
import diff_table

class Derivative:
    def __init__(self, func,initial_step, max_order = 1, diff_func = None) :
        self.func = func 
        self.initial_Step = initial_step
        self.max_Order = max_order
        self.diff_Func = diff_func

    def foward_Diff(self, x, h):
        return (self.func(x+h)-self.func(x))/h
    def backward_Diff(self, x, h):
        return (-self.func(x-h)+self.func(x))/h
    def central_Diff(self, x, h):
        return((self.func(x+h)-self.func(x-h))/(2.0*h))


    def get_func(self):
        return self.func
    def get_diff_Table(self):
        return self.diff_Table
    def get_max_Order(self):
        return self.max_Order
    def get_diff_Func(self):
        return self.diff_Func
    

    def set_func(self,func):
        if callable(func):
             self.func = func
        else:
            print("Cannot take the derivative, not a function!")
            return 
    def set_initial_Step(self,new_Step):
        if not (new_Step.isnumeric() and new_Step > 0):
            print("Step is not an number!")
        else: 
            self.set_initial_Step = new_Step
    def set_max_Order(self,new_Max_Order):
        if new_Max_Order.isnumeric() and new_Max_Order >= 0:
            self.max_order = new_Max_Order
    def set_diff_Func(self, new_Diff_func):
        if callable(new_Diff_func):
            self.diff_Func = new_Diff_func
        else:
            print("The derivative is not a function!")
            return 
        
        



