import numpy as np
import math
import matplotlib.pylab as plt

class SIR:
    def __init__(self, tstart, tend, dt, beta, gamma, N):

        self._tstart = tstart;
        self._dt     = dt;
        self._tend   = tend;
        self._beta   = beta;
        self._gamma  = gamma;
        self._N      = N;
        self._dS     = 0;
        self._dI     = 0;
        self._dR     = 0;


    def Flow(self):
        self._dS = -self._beta*self._I*self._S/self._N;
        self._dI =  self._beta*self._I*self._S/self._N - self._gamma*self._I;
        self._dR =  self._gamma*self._I;

    def Update(self):
        self.Flow();
        self._S = self._S + self._dt*self._dS;
        self._I = self._I + self._dt*self._dI;
        self._R = self._R + self._dt*self._dR;

    def UpdateRK(self):
        # produce k1
        self.Flow();
        ys, yi, yr = self._S, self._I, self._R;
        ks, ki, kr = self._dS, self._dI, self._dR;
        # produce k2
        self._S = ys + self._dS*self._dt/2.0;
        self._I = yi + self._dI*self._dt/2.0;
        self._R = yr + self._dR*self._dt/2.0;
        ks = ks + 2*self._dS;
        ki = ki + 2*self._dI;
        kr = kr + 2*self._dR;
        self.Flow();
        # produce k3
        self._S = ys + self._dS*self._dt/2.0;
        self._I = yi + self._dI*self._dt/2.0;
        self._R = yr + self._dR*self._dt/2.0;
        ks = ks + 2*self._dS;
        ki = ki + 2*self._dI;
        kr = kr + 2*self._dR;
        self.Flow();
        # produce k4
        self._S = ys + self._dS*self._dt;
        self._I = yi + self._dI*self._dt;
        self._R = yr + self._dR*self._dt;
        ks = ks + self._dS;
        ki = ki + self._dI;
        kr = kr + self._dR;
        # advance
        self._S = ys + ks*self._dt/6.0;
        self._I = yi + ki*self._dt/6.0;
        self._R = yr + kr*self._dt/6.0;

    def Initialize(self,S0,I0,R0):
        self._S    = S0;
        self._I    = I0;
        self._R    = R0;
        self.Flow()

        self._Time = np.arange(self._tstart,self._tend,self._dt);

        self._SS   = np.zeros(len(self._Time));
        self._II   = np.zeros(len(self._Time));
        self._RR   = np.zeros(len(self._Time));

        self._dSS  = np.zeros(len(self._Time));
        self._dII  = np.zeros(len(self._Time));
        self._dRR  = np.zeros(len(self._Time));

    def Simulate(self):
        for ii in range(len(self._Time)):
            self._SS[ii] = self._S;
            self._II[ii] = self._I;
            self._RR[ii] = self._R;

            self._dSS[ii] = self._dS;
            self._dII[ii] = self._dI;
            self._dRR[ii] = self._dR;

            self.Update();

    def PlotSIR(self, num):
        plt.figure(num)
        ps = plt.plot(self._Time,self._SS,label="Susceptible");
        plt.setp(ps, 'Color', 'b', 'linewidth', 3)
        pi = plt.plot(self._Time,self._II,label="Infected");
        plt.setp(pi, 'Color', 'r', 'linewidth', 3)
        pr = plt.plot(self._Time,self._RR,label="Recovered");
        plt.setp(pr, 'Color', 'g', 'linewidth', 3)
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=3, mode="expand", borderaxespad=0.)
        plt.grid(True)

    def PlotDSIR(self, num):
        plt.figure(num)
        ps = plt.plot(self._Time,self._dSS,label="Susceptible");
        plt.setp(ps, 'Color', 'b', 'linewidth', 3)
        pi = plt.plot(self._Time,self._dII,label="Infected");
        plt.setp(pi, 'Color', 'r', 'linewidth', 3)
        pr = plt.plot(self._Time,self._dRR,label="Recovered");
        plt.setp(pr, 'Color', 'g', 'linewidth', 3)
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=3, mode="expand", borderaxespad=0.)
        plt.grid(True)


if __name__ == "__main__":
    sir = SIR(0,8,.01,.8,1.2,100);
    sir.Initialize(999,1,0);
    sir.Simulate();
    sir.PlotSIR(1)
    plt.show()
