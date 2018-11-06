%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt

def func(x):
    a = -2
    b = 10
    return np.exp(a*x) * np.cos(b*x)

def func_integral(x):
    a = -2
    b = 10
    c = 5
    d = 52
    return np.exp(a*x) * (c* np.sin(b*x) - np.cos(b*x)) * (1/d)








def trapezoid_core(f,x,h):
    return 0.5 * h * (f(x + h) + f(x))

def trapezoid_method(f,a,b,N):
    x = np.linspace(a,b,N)
    h = x[1] - x[0]
    
    fint = 0.0
    
    for i in range(0,len(x) - 1, 1):
        fint += trapezoid_core(f,x[i],h)
    
    print("Iterations = " + str(i))
    return(fint)
    
    

    
def simpson_core(f,x,h):
    return (h/3)* (f(x) + 4*f(x+h) + f(x+2*h))

def simpsons_method(f,a,b,N):
    x = np.linspace(a,b,N)
    h = x[1] - x[0]
    
    Fint = 0.0
    
    for i in range (0, len(x) - 2, 2):
        Fint += simpson_core(f,x[i],h)
        
    if (N%2 == 0):
        Fint += simpson_core(f,x[-2], 0.5 * h)
        
    print("Iterations = " + str(i))   
    return Fint




def romberg_core(f,a,b,i):
    h = b - a
    dh = h/2.**i
    K = h/2.**(i+1)
    
    M = 0.0
    for j in range(2**i):
        M += f(a + 0.5*dh + j*dh)
        
    return(K*M)


def romberg_integration(f,a,b,tol):
    
    i = 0 
    
    imax = 1000
    
    delta = 100.0*np.fabs(tol)
    
    I = np.zeros(imax,dtype=float)
    
    I[0] = 0.5*(b-a)*(f(a) + f(b))

    i +=1
    
    while(delta>tol):
        
        I[i] = 0.5*I[i-1] + romberg_core(f,a,b,i)
        
        delta = np.fabs( (I[i] - I[i-1]) / I[i])
        
        print(i,I[i],I[i-1],delta)
        
        if (delta>tol):
            i += 1
            
            if (i > imax):
                print('Max iterations reached.')
                raiseStopIteration('Stopping iterations after ', i)
        
        
    print("Iterations = " + str(i))
    return(I[i])







print("function = e^(-2x) * cos(10x)")
answer = func_integral(np.pi) - func_integral(0)
print("Integral = " + str(answer))
print(" ")
print('Trapezoid method')
print(trapezoid_method(func,0,np.pi,10))
print(" ")
print("Simpson's method")
print(simpsons_method(func,0,np.pi,10))
print(" ")
print("Romberg's method")
tolerance = 1e-6
RI = romberg_integration(func,0,np.pi,tolerance)
print(RI, (RI - answer)/answer, tolerance