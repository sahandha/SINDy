import warnings
import numpy as np
from scipy.signal import filtfilt, butter
import matplotlib.pyplot as plt
import mysql.connector

gr ,  re,  bl,  yl,  bk = '#32AF4B', '#AF324B', '#323BAF', '#AFAF4B', '#000000'
lgr, lre, lbl, lyl, lbk = '#ccffcc', '#ffcccc', '#ccccff', '#AFAFBB', '#333333'

class DataBaseObject:
    def __init__(self,usrnm='root',psswd='',hst='localhost',db=''):
        self._Connection = mysql.connector.connect(user=usrnm, password=psswd,host=hst,database=db)
        self._Cursor = self._Connection.cursor()
        self.getHistoryData()

    def getHistoryData(self):
        sqlquery = 'SELECT * FROM HistoryData';
        self._Cursor.execute(sqlquery)
        data = self._Cursor.fetchall()
        self._T = [row[1] for row in data]
        self._S = [row[2] for row in data]
        self._I = [row[3] for row in data]
        self._P = [row[4] for row in data]
        self._R = [row[5] for row in data]
        self._D = [row[6] for row in data]
        self._B = [row[7] for row in data]
        self._N = [row[8] for row in data]

    def getSINDyReadyData(self):
        return np.array([
                self._T,
                self._S,
                [x+y for x,y, in zip(self._P,self._I)],
                self._R]).transpose()

    def DBPlot(self):
        p = plt.plot(self._T,self._S,'-o',label="Susceptible")
        plt.setp(p, 'Color', 'b', 'linewidth', 3)

        p = plt.plot(self._T,np.array(self._P)+np.array(self._I),'-o',label="Infected")
        plt.setp(p, 'Color', 'r', 'linewidth', 3)

        p = plt.plot(self._T,self._R,'-o',label="Recovered")
        plt.setp(p, 'Color', 'g', 'linewidth', 3)

        p = plt.plot(self._T,np.array(self._S)+np.array(self._P)+np.array(self._I)+np.array(self._R),label="Total")
        plt.setp(p, 'Color', 'k', 'linewidth', 3)
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=3, mode="expand", borderaxespad=0.)
class Differentiate:
    def __init__(self,data,order=1,method='CentralDifference',filter=None):

        self._t          = data[:,0]
        self._x          = data[:,1:]
        self._dt         = np.diff(self._t[0:2])[0]
        self._dim        = data.shape[1]-1
        self._method     = method
        self._filter     = filter
        print "filter: ", self._filter
        print "method: ", self._method

        self.FilterChoice()

    def CentralDifference(self):
        self._dx = np.zeros_like(self._x)
        self._dx[0   ,:] = (self._x[1   ,:]-self._x[0,:])/self._dt
        self._dx[1:-2,:] = (self._x[2:-1,:]-self._x[0:-3,:])/(2*self._dt)
        self._dx[-1  ,:] = (self._x[-1  ,:]-self._x[-2,:])/self._dt

    def Derive(self):
        if self._method=="CentralDifference":
            self.CentralDifference()
        elif self._method=="TotalVariationRegularized":
            self.TotalVariationRegularized()
        else:
            print "Please provide a valid method."

    def FilterChoice(self):
        if self._filter==None:
            self.Derive()
        elif self._filter=="IIRFilter":
            self.IIRFilter()
        elif self._filter=="RunningMean":
            self.RunningMean()
        else:
            print "please provide a valid filter method"

    def RunningMean(self, N=10):
        for i in range(self._x.shape[1]):
            self._x[:,i] = np.convolve(self._x[:,i], np.ones((N,))/N)[(N-1):]
        self.Derive()

    def TotalVariationRegularized(self):
        print "To be implemented"

    def IIRFilter(self,order=3, cutoff=0.02):

        b, a = butter(order, cutoff)
        for i in range(self._x.shape[1]):
            self._x[:,i] = filtfilt(b,a,self._x[:,i])
        self.Derive()

    def DerivativePlot(self,fignum=1,labels=[],colors=[],**kwargs):
        plt.figure(fignum)
        plt.suptitle(kwargs['title'])
        plt.subplots_adjust(top=0.85)
        for i in range(self._dim):
            plt.subplot(2,1,1)
            ps = plt.plot(self._t,self._x[:,i],label=labels[i]);
            plt.setp(ps, 'Color', colors[i], 'linewidth', 3)
            plt.subplot(2,1,2)
            ps = plt.plot(self._t,self._dx[:,i],label=labels[i]);
            plt.setp(ps, 'Color', colors[i], 'linewidth', 3)
        plt.subplot(2,1,1)
        plt.ylabel('Data')
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=3, mode="expand", borderaxespad=0.)
        plt.grid(True)
        plt.subplot(2,1,2)
        plt.ylabel('Derivative Data')
        plt.xlabel('Time')
        plt.grid(True)


if __name__ == "__main__":
    db = DataBaseObject(db='sim_v1_1_16_2017')

    db.DBPlot()
