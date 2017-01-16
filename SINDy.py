import warnings
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
from collections import Counter

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
        r.set_initial_value(np.mean(self._x[0:int(0.1/dt)],axis=0), t0)

        while r.successful() and r.t < te-2*dt: #TODO: figure out a better way to do this.
            self._rx.append(np.insert(r.integrate(r.t+dt),0,r.t+dt))

        self._rx = np.array(self._rx)

    def SINDyPlot(self, fignum=1, statesymbols=[],datacolors=[],simcolors=[]):

        if len(statesymbols) != self._dim:
            warnings.warn("The number of state symbols don't match the state dimension.")
            statesymbols = np.arange(self._dim)+1
        if len(datacolors) != self._dim:
            warnings.warn("The number of color specs don't match the state dimension.")
            datacolors = ['b','r','g']

        plt.figure(fignum)
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

    def StringConvolve(self,l1,l2):
        res = []
        for i in xrange(len(l1)):
            for j in xrange(min(i,len(l2)-1),len(l2)):
                if l1[i] == "1":
                    res.append(l2[j])
                else:
                    res.append(l1[i]+l2[j])
        return res

    def StringMultFormat(self,str):
        sortedStr = sorted(str)
        availablesVars = list(set(sortedStr))
        counts = Counter(sortedStr)

        rs = ''

        for i in range(len(availablesVars)):
            power = counts[availablesVars[i]]
            if power == 1:
                rs = rs + availablesVars[i] + " "
            else:
                rs = rs + availablesVars[i]+"^{} ".format(power)
        return rs

    def StringTerms(self,strlist):
        terms = []
        s1 = "1"
        for i in self._polyorder:
            terms.join(self.StringConvolve(s1,strlist))
            s1 = terms

    def StringModelView(self,StateVariables=[]):
        terms = []
        s1 = ["1"]
        for i in range(self._polyorder):
            terms += self.StringConvolve(s1,StateVariables)
            s1 = terms
        terms = ["1"] + terms
        for i in xrange(len(StateVariables)):
            if self._xi[0,i] > 1e-3:
                row = "d"+StateVariables[i]+"/dt = " + "{: 2.3f}".format(self._xi[0,i])
                sp = " + "
            else:
                row = "d"+StateVariables[i]+"/dt = "
                sp = ""

            for j in xrange(1,len(self._xi)):
                ss = "{: 2.2f}".format(self._xi[j,i])
                if self._xi[j,i]>1e-3:
                    row = row + sp + "{: 2.3f}".format(self._xi[j,i])+" "+self.StringMultFormat(terms[j])
                    sp = " + "
            print row


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
    eps    = 0.05
    noise  = eps*np.random.randn(3,sir._Time.shape[0])
    dnoise = eps*np.random.randn(3,sir._Time.shape[0])

    data  = np.transpose(np.insert(np.array([sir._SS , sir._II , sir._RR]) + noise,0,sir._Time,axis=0))
    ddata = np.transpose(np.insert(np.array([sir._dSS, sir._dII, sir._dRR])+dnoise,0,sir._Time,axis=0))


    '''
        Run SINDy
    '''
    sin = SINDy(data=data,polyorder=2,usesine=False,cutoff=0.1)
    sin.SetDerivative(ddata)
    sin.RunSINDy(simulate=True)
    sin.SINDyPlot(statesymbols=["Susceptible","Infected","Recovered"],
              datacolors=[[0.8, 0.8, 1.0],[1.0, 0.8, 0.8],[0.8, 1.0, 0.8]],
              simcolors =[[0.0, 0.0, 1.0],[1.0, 0.0, 0.0],[0.0, 1.0, 0.0]])

    sin.StringModelView(["S","I","R"])
