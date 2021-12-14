from sympy import Symbol, Derivative,sympify,lambdify
import time


def ar_h(func,initial,cons,eps=0.01,max_iter=10000):

    """function to find solution for a problem using the arrow hurwicz algorithm

    Parameters
    ----------
    func : String
         function we want to minimize.
    initial : Dict
         Dictionary containing the intial guesses for theta and the langrange multiplier.
    cons : String
         constraint function
    eps : float
         epsilon value
    max_iter : int
         maximum iterations before the algorithm aborts
    Returns
    -------
    float
        resulting theta from the algorithm.
    """


    #convert string to sympy expression
    func=sympify(func)
    cons=sympify(cons)
    #create symbol for x and l
    x= Symbol('x')
    l=Symbol('l')
    # create lagarngian
    lang=func+l*cons
    #derivate lagrangian
    deriv= Derivative(lang, x).doit()
    #intialize looping variables and get intial guess
    p=10
    count=0
    xn=initial['x']
    ln=initial['l']
    #convert sympy expression to lambda function for faster execution
    deriv=lambdify( [x,l], deriv, "math" )
    cons=lambdify( [x], cons, "math" )
    #run loop
    start=time.time()
    #break look if gradient is very or max iter reached 
    while p>0.0001 and count<max_iter:
        diff=deriv(xn,ln)
        #find next lambda
        ln=max(0,ln+eps*cons(xn))
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
func='0.5*x**2'
initial={'x':2,'l':1}
cons='5-x'

res=ar_h(func,initial,cons)

print(res)
        

