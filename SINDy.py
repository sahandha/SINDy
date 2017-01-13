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

    def setDerivative(self,derdata):
        self._dx = derdata[:,1:];

    def PoolData(self,x):

        yout = np.ones((np.shape(x)[0],1))

        for i in xrange(self._polyorder):
            yout = np.append(yout,self.PolyConvolve(yout,x),axis=1)

        if self._usesine:
            for k in xrange(1,10):
                yout = np.append(np.append(yout, np.sin(k*x),axis=1),np.cos(k*x),axis=1)
        return yout

    def SparsifyDynamics(self):
        self._xi = np.linalg.lstsq(self._theta,self._dx)

        for k in xrange(10):
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

    def ComputeDerivatives(self):
        if self._dx == []:
            print "Write numerical computation of derivatives and set it to self._dx"

            print "no derivative data provided..."
            print "computing derivatives using ... method..."


    def SimulateSINDy(self):
        dt = np.diff(self._t[0:2])[0]
        te = self._t[-1]
        r = ode(self.SparseGalerkin).set_integrator('dopri5')
        r.set_initial_value(self._x[0], self._t[0])

        while r.successful() and r.t < te:
            print (r.t+dt, r.integrate(r.t+dt))


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
            plt.grid()
            plt.legend(loc=0)
        plt.show()



if __name__ == "__main__":

    data  = np.array([[1,.2,3,.3],[2.1,4,5,2],[3,4,6,.6]])
    ddata = np.array([[2,2,3,4],[3,4,5,3],[3,4,5,1.2]])

    sindy = SINDy(data=data,polyorder=1,usesine=False,cutoff=.5)
    sindy.setDerivative(ddata)
    sindy.RunSINDy()
    #sindy.SINDyPlot(statesymbols=["S","I","R"])
