import warnings
import numpy as np

class SINDy:
    def __init__(self,data,polyorder=2,usesine=False,cutoff=0.025):

        self._x         = data
        self._dx        = []
        self._dim       = np.shape(data)[1]
        self._polyorder = polyorder
        self._usesine   = usesine
        self._cutoff    = cutoff
        self._theta     = []
        self._xi        = []

    def setddata(self,derdata):
        self._dx = derdata;

    def PoolData(self):

        yout = np.ones((np.shape(self._x)[0],1))

        for i in xrange(self._polyorder):
            yout = np.append(yout,self.PolyConvolve(yout,self._x),axis=1)

        if self._usesine:
            for k in xrange(1,10):
                yout = np.append(np.append(yout, np.sin(k*self._x),axis=1),np.cos(k*self._x),axis=1)

        return yout

    def SparsifyDynamics(self):
        self._xi = np.linalg.lstsq(self._theta,self._dx)

        for k in xrange(10):
            smallinds = np.abs(self._x)<self._cutoff
            self._x[smallinds] = 0

            for ind in xrang(self._dim):
                biginds = !smallinds(:,ind)
                self._xi[biginds,ind] = np.linalg.lstsq(self._theta(:,biginds),self._dx(:,inds))


    #def SparseGalerkin(self):

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
        self._theta = self.PoolData()

    def ComputeDerivatives(self):
        if self._dx == []:
            print "Write numerical computation of derivatives and set it to self._dx"

            print "no derivative data provided..."
            print "computing derivatives using ... method..."

    def RunSINDy(self):

        self.ComputeDerivatives()
        self.ComputeTheta()
        self.SparsifyDynamics()




if __name__ == "__main__":
    data = np.array([[1,2,3],[3,4,5]])
    yout = np.ones((2,1))

    sindy = SINDy(data=data,polyorder=1,usesine=False)
    res = sindy.PoolData()

    print data
    print res
    print sindy._cutoff
