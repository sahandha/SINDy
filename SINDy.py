import warnings
import numpy as np

class SINDy:
    def __init__(self, data, polyorder ):

        self.data      = data
        self.dim       = np.shape(data)[1]
        self.polyorder = polyorder
        self.Theta     = []


    def PoolData(self):

        yout = np.ones((np.shape(self.data)[0],1))

        for i in xrange(self.polyorder):
            print self.PolyConvolve(yout,self.data)
            yout = np.append(yout,self.PolyConvolve(yout,self.data),axis=1)

        return yout

    #def SparsifyDynamics(self):

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


if __name__ == "__main__":
    data = np.array([[1,2,3],[3,4,5]])
    yout = np.ones((2,1))

    sindy = SINDy(data,1)
    res = sindy.PoolData()

    print data
    print res
