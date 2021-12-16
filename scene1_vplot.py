from sympy import Symbol, Derivative,sympify,lambdify
import time
from math import erf,sqrt
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt

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


 
    yn=initial['y']
    an=initial['a']
    pn=initial['p']
    deriv=lambda x,y,l,a:(-2/(x**3))+l*((sqrt(a)/x**2)*stats.norm(y,x).pdf(y+sqrt(a))-(-sqrt(a)/x**2)*stats.norm(y,x).pdf(y-sqrt(a)))
    cons=lambda x,y,a,p:p-phi(y+sqrt(a),y,x) + phi(y-sqrt(a),y,x)
    #intialize looping variables and get intial guess
    x=np.arange(0.5,11,0.2)
    y=np.arange(0.5,11,0.2)
    derivs=[]
    conss=[]
    for xer in x:
        for yer in y:
            derivs.append(-deriv(xer,yn,yer,an))
            conss.append(cons(xer,yn,an,pn))
    return x,y,derivs,conss

   


#test function out
initial={'y':4,'a':1,'p':0.7}
x,y,dx,dy=ar_h_1(initial)

Y,X=np.meshgrid(y,x)
dx=np.array(dx)
dy=np.array(dy)
n=2
color = np.sqrt(((dx+n)/2)*2 + ((dy+n)/2)*2)
dx/=np.linalg.norm(dx)
dy/=np.linalg.norm(dy)
# dx=dx*2
# dy=dy*2

# print(X)
# print(Y)

# n = 2
fig,ax=plt.subplots(figsize=(16, 12), dpi=80)

# Creating plot
ax.quiver(X, Y, dx, dy,color,angles='xy')
# ax.xaxis.set_ticks([])
# ax.yaxis.set_ticks([])
# ax.set_aspect('equal')
ax.set_title('gradient')
 
 

# show figure
plt.tight_layout()
plt.show()

        

