import warnings
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt

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

    def setDerivative(self,derdata):
        self._dx = derdata[:,1:]
        print "**** Derivative Set ****"

    def PoolData(self,x):

        yout = np.ones((np.shape(x)[0],1))

        for i in xrange(self._polyorder):
            yout = np.append(yout,self.PolyConvolve(yout,x),axis=1)

        if self._usesine:
            for k in xrange(1,10):
                yout = np.append(np.append(yout, np.sin(k*x),axis=1),np.cos(k*x),axis=1)
        return yout

    def SparsifyDynamics(self):
        print "**** Performing regression. Please wait... ****"
        self._xi = np.linalg.lstsq(self._theta,self._dx)

        for k in xrange(10):
            print "\tIteration {} in progress...".format(k)
            self._x[np.abs(self._x)<self._cutoff] = 0
            self._xi = np.linalg.lstsq(self._theta,self._dx)

    def SparseGalerkin(self,t,y):
        yy = np.array(y).reshape((1,self._dim))
        yPool = self.PoolData(yy)
        return np.dot(yPool,self._xi[0])

    def PolyConvolve(self, data1, data2):

        if (np.shape(data1)[0] != np.shape(data2)[0]):
            warnings.warn("Data dimension mismatch.")

        row = np.zeros(np.shape(data1)[1]*np.shape(data2)[1])
        res = [];

        for n in xrange(np.shape(data1)[0]):
            a, b = data1[n], data2[n]
            ind = 0
            for i in xrange(len(a)):
                for j in xrange(len(b)):
                    row[ind] = a[i]*b[j]
                    ind += 1
            res.append(list(row))
        return np.array(res)

    def ComputeTheta(self):
        self._theta = self.PoolData(self._x)
        print "**** Candidate functions library has been created ****"

    def ComputeDerivatives(self):
        if self._dx == []:
            print "Write numerical computation of derivatives and set it to self._dx"

            print "no derivative data provided..."
            print "computing derivatives using ... method..."

    def SimulateSINDy(self):
        print "**** Identification is complete. We now use it to simulate the system. ****"
        dt = np.diff(self._t[0:2])[0]
        te = self._t[-1]
        r = ode(self.SparseGalerkin).set_integrator('dopri5', nsteps=3000)
        r.set_initial_value(self._x[0], self._t[0])

        while r.successful() and r.t < te:

            self._rx.append(np.insert(r.integrate(r.t+dt),0,r.t+dt))
            #print "results are in: ", (r.t+dt, r.integrate(r.t+dt))
        self._rx = np.array(self._rx)
    def RunSINDy(self):

        self.ComputeDerivatives()
        self.ComputeTheta()
        self.SparsifyDynamics()
        self.SimulateSINDy()

    def SINDyPlot(self, statesymbols=[]):

        if len(statesymbols) != self._dim:
            warnings.warn("The number of state symbols don't match the state dimension.")
            statesymbols = np.arange(self._dim)+1

        plt.figure
        for i in xrange(self._dim):
            plt.plot(self._t,self._x[:,i],label="{}".format(statesymbols[i]))

            plt.plot(self._rx[:,0],self._rx[:,i+1],label="{}".format(statesymbols[i]))
            plt.grid()
            plt.legend(loc=0)
        plt.show()



if __name__ == "__main__":

    from SIR import *

    sir = SIR(tstart=0, tend=5, dt=.1, beta=3, gamma=2, N=1)
    sir.Initialize(S0=0.9, I0=0.1, R0=0);
    sir.Simulate();


    data  = np.transpose(np.array([sir._Time, sir._SS, sir._II, sir._RR]))
    ddata = np.transpose(np.array([sir._Time, sir._dSS, sir._dII, sir._dRR]))

    sin = SINDy(data=data,polyorder=2,usesine=False)
    sin.setDerivative(ddata)
    sin.RunSINDy()
    sin.SINDyPlot(["S","I","R"])
