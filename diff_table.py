import derivative
import numpy as np

class Diff_Table: 
    def __init__(self, derivative, initial_step,parameter_func,point):
        self.diff = derivative
        self.max_Order = derivative.max_Order
        self.delta_Forward = [] 
        self.delta_Central =  [] 
        self.delta_Backward = [] 
        self.parameter_func = parameter_func
        self.elements = []
        self.steps = []
        self.rsCoef = []
        self.size = len(self.steps)
        self.initial_step = initial_step
        self.point = point
 
    def build_delta_Foward(self):
        self.delta_Forward = np.zeros(self.size) 
        for i in range(self.size):
            self.delta_Forward[i] = self.diff.foward_Diff(self.point,  self.steps[i])
    
    def build_delta_Central(self):
        self. delta_Central = np.zeros(self.size) 
        for i in range(self.size):
            self.delta_Central[i] = self.diff.central_Diff(self.point,  self.steps[i])
    
    def build_delta_Backward(self):
        self. delta_Backward = np.zeros(self.size) 
        for i in range(self.size):
            self.delta_Backward[i] = self.diff.backward_Diff(self.point, self.steps[i])

    def build_Table_Central(self):
        
        b = self.size
        for k in range(b,self.max_Order):
            n = k + 2

            if k == 0:
                g0 = self.parameter_func(0.0) 
                g1 = self.parameter_func(1.0)
                self.steps.append(1.0/g0*self.initial_step)
                self.steps.append(1.0/g1*self.initial_step)
                self.rsCoef.append(1.0/(1-g1/g0))
                self.rsCoef.append(1.0/(1-g0/g1))
                self.size += 2
                self.richardsonStepCentral()         
            
            else:
                l = n-1 
                
                gl = self.parameter_func(l)
                dh = np.zeros(n)
                lhlambda = self.rsCoef
                dlambda = np.zeros(n)
                lambda_l = 1.0

                for i in range(l):
                    gi = self.parameter_func(i)
                    lambda_i =lhlambda[i]/ (1.0 - gl / gi)
                    dlambda[i] = lambda_i
                    lambda_l *= 1.0 / (1.0 - gi / gl)
                    dh[i] = self.steps[i]
                    
                dh[l] = 1.0/gl
                self.steps.append(dh[-1]*self.initial_step)
                dlambda[l] = lambda_l
                self.rsCoef= dlambda
                self.size += 1
                self.richardsonStepCentral()         
    def build_Table_Foward(self):
        
        b = self.size
        for k in range(b,self.max_Order):
            n = k + 2

            if k == 0:
                g0 = self.parameter_func(0.0) 
                g1 = self.parameter_func(1.0)
                self.steps.append(1.0/g0*self.initial_step)
                self.steps.append(1.0/g1*self.initial_step)
                self.rsCoef.append(1.0/(1-g1/g0))
                self.rsCoef.append(1.0/(1-g0/g1))
                self.size += 2
                self.richardsonStepCentral()         
            
            else:
                l = n-1 
                
                gl = self.parameter_func(l)
                dh = np.zeros(n)
                lhlambda = self.rsCoef
                dlambda = np.zeros(n)
                lambda_l = 1.0

                for i in range(l):
                    gi = self.parameter_func(i)
                    lambda_i =lhlambda[i]/ (1.0 - gl / gi)
                    dlambda[i] = lambda_i
                    lambda_l *= 1.0 / (1.0 - gi / gl)
                    dh[i] = self.steps[i]
                    
                dh[l] = 1.0/gl
                self.steps.append(dh[-1]*self.initial_step)
                dlambda[l] = lambda_l
                self.rsCoef= dlambda
                self.size += 1
                self.richardsonStepFoward()         
    def build_Table_Backward(self):
        
        b = self.size
        for k in range(b,self.max_Order):
            n = k + 2

            if k == 0:
                g0 = self.parameter_func(0.0) 
                g1 = self.parameter_func(1.0)
                self.steps.append(1.0/g0*self.initial_step)
                self.steps.append(1.0/g1*self.initial_step)
                self.rsCoef.append(1.0/(1-g1/g0))
                self.rsCoef.append(1.0/(1-g0/g1))
                self.size += 2
                self.richardsonStepCentral()         
            
            else:
                l = n-1 
                
                gl = self.parameter_func(l)
                dh = np.zeros(n)
                lhlambda = self.rsCoef
                dlambda = np.zeros(n)
                lambda_l = 1.0

                for i in range(l):
                    gi = self.parameter_func(i)
                    lambda_i =lhlambda[i]/ (1.0 - gl / gi)
                    dlambda[i] = lambda_i
                    lambda_l *= 1.0 / (1.0 - gi / gl)
                    dh[i] = self.steps[i]
                    
                dh[l] = 1.0/gl
                self.steps.append(dh[-1]*self.initial_step)
                dlambda[l] = lambda_l
                self.rsCoef= dlambda
                self.size += 1
                self.richardsonStepBackward()         


    def richardsonStepCentral(self):
        self.elements = np.append(self.elements,0)
        self.build_delta_Central()
        for i in range(self.size):
            self.elements[-1] += self.rsCoef[i]*self.delta_Central[i] 
    def richardsonStepFoward(self):
        self.elements = np.append(self.elements,0)
        self.build_delta_Foward()
        for i in range(self.size):
            self.elements[-1] += self.rsCoef[i]*self.delta_Forward[i] 
    def richardsonStepBackward(self):
        self.elements = np.append(self.elements,0)
        self.build_delta_Backward()
        for i in range(self.size):
            self.elements[-1] += self.rsCoef[i]*self.delta_Backward[i]    
            
    def printTable(table):
        for i in range(7*table.max_Order):
            print("-",end='-')
        print()
        print("Steps:                  |",end=" ")
        for i in range(table.max_Order):    
            print(format (table.steps[i],".6f"),end= ' | ')
        print()
        for i in range(7*table.max_Order):
            print("-",end='-')
        print()
        print("Richardson Coeficients: |",end=" ")
        for i in range(table.max_Order):    
            print(format(table.rsCoef[i],".6f"),end= ' | ')
        print()
        for i in range(7*table.max_Order):
            print("-",end='-')
        print()
        print("Derivatives:            |",end=" ")
        for i in range(table.max_Order):    
            print(format(table.elements[i],".6f"),end= ' | ')   
        print()
        for i in range(7*table.max_Order):
            print("-",end='-')