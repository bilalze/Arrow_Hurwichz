from sympy import Symbol, Derivative,sympify,lambdify
import time
from math import erf,sqrt
from scipy import stats




def ar_h_1(initial,eps=0.01,max_iter=10000):

    """function to find solution for a problem using the arrow hurwicz algorithm

    Parameters
    ----------
    initial : Dict
         Dictionary containing the intial guesses for theta and the langrange multiplier.
    eps : float
         epsilon value
    max_iter : int
         maximum iterations before the algorithm aborts
    Returns
    -------
    float
        resulting theta from the algorithm.
    """


    xn=initial['x']
    ln=initial['l']
    yn=initial['y']
    an=initial['a']
    pn=initial['p']
    deriv=lambda x,y,l,p:(-2/(x**3))+l*(-2*(stats.norm(y,x).ppf((p+1)/2)-y)*((stats.norm(y,x).ppf((p+1)/2)-y)/(x**2)))
    cons=lambda x,y,a,p:(stats.norm(y,x).ppf((p+1)/2)-y)**2-a
    #intialize looping variables and get intial guess
    pcheck=10
    count=0
    #run loop
    start=time.time()
    #break look if gradient is very or max iter reached 
    while pcheck>0.00001 and count<max_iter:
        diff=deriv(xn,yn,ln,pn)
        #find next lambda
        # aa=(-2/(xn**3))
        # bb=-2*(stats.norm(yn,xn).ppf((pn+1)/2)-yn)
        # cc=(stats.norm(yn,xn).ppf((pn+1)/2)-yn)/(xn**2)
        # dd=bb*cc
        ln=max(0,ln+eps*cons(xn,yn,an,pn))
        #find next theta
        xn=xn-eps*diff
        #magnitude of gradient
        pcheck=abs(diff)
        count+=1
    print('num of iter=',count)
    print('magnitude of gradient=',p)
    print('time taken=',time.time()-start)
    return xn


#test function out
initial={'x':4,'y':4,'l':1,'a':1,'p':0.7}
res=ar_h_1(initial)

print(res)
        

