import numpy as np
import matplotlib.pyplot as plt

class Problem():
    def __init__(self):
        
        self.m = 1      # mass [kg]
        self.c = 0.8      # damping constant [N/(m/s)]
        self.k = 5      # spring constant [N/m]

        self.x0 = 1     # initial position [m]
        self.v0 = 0     # initial velocity [m]

        self.t0 = 0     # initial time [s]
        self.tf = 20    # final time [s]
        self.dt = 0.01  # timestep [s]
        self.N = int((self.tf-self.t0)/self.dt + 1)

        self.t_array = np.linspace(self.t0,self.tf,self.N)     # t array
        self.x_array = np.zeros_like(self.t_array)              # x array
        self.v_array = np.zeros_like(self.t_array)              # v array

        self.x_array[0] = self.x0
        self.v_array[0] = self.v0

        self.solve()
        self.plot_solution()

    def force(self,t):
        # self.Amplitude = self.k*self.x0
        # self.Amplitude = 0
        # self.omega = np.sqrt(self.k/self.m+0.01)
        # self.F = self.Amplitude*np.sin(self.omega*t)
        self.F = 0
        return self.F

    def solve(self):
        #   m*a + c*v + k*x = F 
        #   <=> a = - k/m*x - c/m*v + F/m
        #   <=> {   dx/dt = + 0*x + 1*v + 0
        #       {   dv/dt = - k/m*x - c/m*v + F/m
        #   <=> X = [[x],
        #            [v]]
        #       dX/dt = [[dx/dt],
        #                [dv/dt]]
        #       A = [[0,1],
        #            [-k/m,-c/m]]
        #       b = [[0],
        #            [F/m]]
        #       dX/dt = A*X + b

        self.A = np.array([[0,1],
                           [-self.k/self.m,-self.c/self.m]])
        self.X = np.array([[self.x0],
                           [self.v0]])
        for i, t in enumerate(self.t_array[1:]):
            self.b = np.array([[self.force(t)/self.m],
                               [0]])
            self.dXdt = np.matmul(self.A,self.X)+self.b
            self.X = self.integrateFirstOrderEuler(self.X,self.dXdt,self.dt)
            #print(f"t = {'%.4f' %t}")
            #print(np.array2string(self.X, formatter={'float_kind':lambda x: "%.4f" % x}))
            self.x_array[i] = self.X[0,0]
            self.v_array[i] = self.X[1,0]
    
    def integrateFirstOrderEuler(self,X,dXdt,dt):
        # x[n+1] = x[n] + dx/dt[n]*dt
        return X+dXdt*dt

    def plot_solution(self):
        plt.figure(1)
        plt.plot(self.t_array, self.x_array, label = "x")
        plt.plot(self.t_array, self.v_array, label = "v")
        plt.legend()
        plt.xlabel("time [s]")
        plt.ylabel("position [m] and velocity [m/s]")
        plt.show()

p=Problem()