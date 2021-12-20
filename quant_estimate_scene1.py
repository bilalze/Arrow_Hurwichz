
import time
from math import erf,sqrt
from scipy import stats
import numpy as np


#class for the quantile estimator
class quant_estimator():
    def __init__(self,samples) :
        self.samples= samples.tolist()
        self.samples.sort()
    #function to get the quantile given percent
    def get_q(self,percent):
        quant_est = round((percent)*len(self.samples))
        return self.samples[quant_est]

#function to find the derivative of the quantile
def quant_diff(q1,q2,delta,percent):
    x1=q1.get_q(percent)
    x2=q2.get_q(percent)
    return (x2-x1)/delta

#convenience function to calculate the derivative of the quantile
def quantifier(mun,sigman,deltan,pn,size):
    samples1=np.random.normal(mun, sigman-deltan, size)
    quant1=quant_estimator(samples1)
    samples2=np.random.normal(mun, sigman+deltan, size)
    quant2=quant_estimator(samples2)
    qd=quant_diff(quant1,quant2,deltan*2,(pn+1)/2)
    return quant1,quant2,qd




def ar_h_1(initial,eps=0.01,max_iter=10000):

    """function to find solution for a problem using the arrow hurwicz algorithm

    Parameters
    ----------
    initial : Dict
         Dictionarmu containing the intial guesses for theta and the langrange multiplier.
    eps : float
         epsilon value
    masigma_iter : int
         masigmaimum iterations before the algorithm aborts
    Returns
    -------
    float
        resulting theta from the algorithm.
    """

    
    sigman=initial['sigma']
    ln=initial['l']
    mun=initial['mu']
    an=initial['a']
    pn=initial['p']
    deltan=initial['delta']
    #quant and qd are the quantile and it derivative 
    deriv=lambda sigma,quant,qd,mu,l,p:(-2/(sigma**3))+l*2*(quant.get_q((p+1)/2)-mu)*(qd)
    cons=lambda quante,mu,a,p:(quante.get_q((p+1)/2)-mu)**2-a
    #intialize looping variables and get intial guess
    pcheck=10
    count=0
    #run loop
    start=time.time()
    constr=0
    #sample size
    size=1000
    #break look if gradient is vermu smalll or masigma iter reached 
    while pcheck>0.0001 and count<max_iter:
        samples=np.random.normal(mun, sigman, size)
        quant=quant_estimator(samples)
        # if count!=0 and count%1000==0:
            # deltan/=2
            # size=int(2*size)
        quant1,quant2,qd=quantifier(mun,sigman,deltan,pn,size)
        diff=deriv(sigman,quant,qd,mun,ln,pn)
        ln=max(0,ln+eps*cons(quant,mun,an,pn))
        #find next theta
        sigman=sigman-eps*diff
        #magnitude of gradient
        pcheck=abs(diff)
        count+=1
    print('num of iter=',count)
    print('magnitude of gradient=',pcheck)
    print('time taken=',time.time()-start)
    return sigman


#test function out
initial={'sigma':4,'mu':4,'l':1,'a':1,'p':0.7,'delta':0.1}
res=ar_h_1(initial)

print(res)
        





        
