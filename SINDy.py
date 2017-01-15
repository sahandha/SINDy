import warnings
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
from progressbar import ProgressBar
from itertools import combinations


class SINDy:
    def __init__(self,data,polyorder=2,usesine=False,cutoff=0.025):

        self._t         = data[:,0]
        self._x         = data[:,1:]
        self._dx        = []
        self._dim       = np.shape(data)[1]-1 # take out the time column
        self._polyorder = polyorder
        self._usesine   = usesine
        self._cutoff    = cutoff
        self._theta     = []
        self._xi        = []
        self._rx        = []

        print "Initiated a class for Sparse Identification from Numerical Dynamics"

    def ComputeDerivatives(self):
        if self._dx == []:
            print "Write numerical computation of derivatives and set it to self._dx"

            print "no derivative data provided..."
            print "computing derivatives using ... method..."

    def ComputeTheta(self):
        self._theta = self.PoolData(self._x)
        print "**** Candidate functions library has been created ****"

    def LeastSquare(self, A, b):
        num_vars = A.shape[1]
        rank = A.shape[0] #np.linalg.matrix_rank(A)
        nonzeros = num_vars*np.shape(b)[1]
        finalsol = 0
        print "Rank: ", rank
        print "Number of variables", num_vars
        print A.shape[0]
        if rank == num_vars:
            print "there is a unique solution"
            sol = np.linalg.lstsq(A, b)[0]    # not under-determined
            finalsol = sol
        elif rank > num_vars:
            print "Over-determined system"
            for nz in combinations(range(rank), num_vars):
                try:
                    sol = np.zeros((num_vars, np.shape(b)[1]))
                    sol = np.asarray(np.linalg.solve(A[nz, :], b[nz,:]))
                    if np.count_nonzero(sol)<nonzeros:
                        finalsol = sol
                        nonzeros = np.count_nonzero(sol)
                except np.linalg.LinAlgError:
                    pass                    # picked bad variables, can't solve
        else:
            print "under-determined system"
            for nz in combinations(range(num_vars), rank):    # the variables not set to zero
                try:
                    sol = np.zeros((num_vars, np.shape(b)[1]))
                    sol[nz, :] = np.asarray(np.linalg.solve(A[:, nz], b))
                    if np.count_nonzero(sol)<nonzeros:
                        finalsol = sol
                        nonzeros = np.count_nonzero(sol)
                except np.linalg.LinAlgError:
                    pass                    # picked bad variables, can't solve
        return finalsol

    def PolyConvolve(self, data1, data2, initialrun=False):

        if (np.shape(data1)[0] != np.shape(data2)[0]):
            warnings.warn("Data dimension mismatch.")

        if not initialrun:
            data1 = data1[:,1:]

        row = []
        res = []

        for n in xrange(np.shape(data1)[0]):
            a, b = data1[n], data2[n]
            for i in xrange(len(a)):
                for j in xrange(min(i,len(b)-1),len(b)):
                    row.append(a[i]*b[j])
            res.append(row)
            row = []

        return np.array(res)

    def PoolData(self,x):
        yout = np.ones((np.shape(x)[0],1))
        initialrun = True
        for i in xrange(self._polyorder):
            yout = np.append(yout,self.PolyConvolve(yout,x,initialrun=initialrun),axis=1)
            initialrun = False
        if self._usesine:
            for k in xrange(1,10):
                yout = np.append(np.append(yout, np.sin(k*x),axis=1),np.cos(k*x),axis=1)
        return yout

    def RunSINDy(self,simulate=True):

        self.ComputeDerivatives()
        self.ComputeTheta()
        self.SparsifyDynamics()
        if simulate:
            self.SimulateSINDy()

    def SetDerivative(self,derdata):
        self._dx = derdata[:,1:]
        print "**** Derivative Set ****"

    def SimulateSINDy(self):
        print "**** Identification is complete. We now use it to simulate the system. ****"
        dt = np.diff(self._t[0:2])[0]
        t0 = self._t[0]
        te = self._t[-1]
        r = ode(self.SparseGalerkin).set_integrator('dopri5', nsteps=1000)
        r.set_initial_value(self._x[0], t0)
        pbar = ProgressBar()
        #while r.successful() and r.t < te:
        for t in pbar(np.arange(t0,te,dt)):
            self._rx.append(np.insert(r.integrate(r.t+dt),0,r.t+dt))

        self._rx = np.array(self._rx)

    def SINDyPlot(self, statesymbols=[],datacolors=[],simcolors=[]):

        if len(statesymbols) != self._dim:
            warnings.warn("The number of state symbols don't match the state dimension.")
            statesymbols = np.arange(self._dim)+1
        if len(datacolors) != self._dim:
            warnings.warn("The number of color specs don't match the state dimension.")
            datacolors = ['b','r','g']

        plt.figure
        for i in xrange(self._dim):
            ps = plt.plot(self._t,self._x[:,i],label="{}".format(statesymbols[i]))
            plt.setp(ps, 'Color', datacolors[i], 'linewidth', 3)
            if self._rx != []:
                ps = plt.plot(self._rx[:,0],self._rx[:,i+1],'--')
                plt.setp(ps, 'Color', simcolors[i], 'linewidth', 3,)

        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=3, mode="expand", borderaxespad=0.)
        plt.grid(True)
        plt.show()

    def SparseGalerkin(self,t,y):
        yy = np.array(y).reshape((1,self._dim))
        yPool = self.PoolData(yy)
        return np.dot(yPool,self._xi)

    def SparsifyDynamics(self):
        print "**** Performing regression. Please wait... ****"
        #self._xi = np.dot(np.linalg.pinv(self._theta),self._dx)
        self._xi = np.linalg.lstsq(self._theta,self._dx)[0]
        print "Iteration in progress: ",
        for k in xrange(10):
            print "{},".format(k+1),
            smallinds = np.abs(self._xi)<self._cutoff
            self._xi[smallinds] = 0
            for ind in range(self._dim):
                biginds = ~smallinds[:,ind]
                #self._xi[biginds,ind] = np.dot(np.linalg.pinv(self._theta[:,biginds]),self._dx[:,ind])
                self._xi[biginds,ind] = np.linalg.lstsq(self._theta[:,biginds],self._dx[:,ind])[0]
        print ""

if __name__ == "__main__":

    from SIR import *

    '''
        Simulate Data
    '''
    sir = SIR(tstart=0.001, tend=10, dt=.01, beta=3, gamma=2, N=1)
    sir.Initialize(S0=0.9, I0=0.1, R0=0);
    sir.Simulate();


    '''
        Add noise to my data
    '''
    import random as rn
    eps = 0.025

    data  = np.transpose(np.array([sir._Time,
                               sir._SS + eps*np.random.randn(sir._SS.shape[0]),
                               sir._II + eps*np.random.randn(sir._SS.shape[0]),
                               sir._RR + eps*np.random.randn(sir._SS.shape[0])]))
    ddata = np.transpose(np.array([sir._Time,
                               sir._dSS + eps*np.random.randn(sir._SS.shape[0]),
                               sir._dII + eps*np.random.randn(sir._SS.shape[0]),
                               sir._dRR + eps*np.random.randn(sir._SS.shape[0])]))


    '''
        Run SINDy
    '''
    sin = SINDy(data=data,polyorder=2,usesine=False)
    sin.SetDerivative(ddata)
    sin.RunSINDy(simulate=True)
    sin.SINDyPlot(statesymbols=["Susceptible","Infected","Recovered"],
                datacolors=[[0.7, 0.7, 1.0],[1.0, 0.7, 0.7],[0.7, 1.0, 0.7]],
                simcolors =[[0.0, 0.0, 1.0],[1.0, 0.0, 0.0],[0.0, 1.0, 0.0]])
