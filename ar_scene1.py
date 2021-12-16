from sympy import Symbol, Derivative,sympify,lambdify
import time
from math import erf,sqrt
from scipy import stats

def phi(x, mu, sigma):
     return (1 + erf((x - mu) / sigma / sqrt(2))) / 2


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
    deriv=lambda x,y,l,a:(-2/(x**3))+l*((sqrt(a)/x**2)*stats.norm(y,x).pdf(y+sqrt(a))-(-sqrt(a)/x**2)*stats.norm(y,x).pdf(y-sqrt(a)))
    cons=lambda x,y,a,p:p-phi(y+sqrt(a),y,x) + phi(y-sqrt(a),y,x)
    #intialize looping variables and get intial guess
    p=10
    count=0
    #run loop
    start=time.time()
    #break look if gradient is very or max iter reached 
    while p>0.00001 and count<max_iter:
        diff=deriv(xn,yn,ln,an)
        #find next lambda
        ln=max(0,ln+eps*cons(xn,yn,an,pn))
        #find next theta
        xn=xn-eps*diff
        #magnitude of gradient
        p=abs(diff)
        count+=1
    print('num of iter=',count)
    print('magnitude of gradient=',p)
    print('time taken=',time.time()-start)
    return xn


#test function out
initial={'x':4,'y':4,'l':1,'a':1,'p':0.7}
res=ar_h_1(initial)

print(res)
        

