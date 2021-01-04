# -*- coding: utf-8 -*-
'''
@author : Alexandre Castel
'''
import numpy as np
import matplotlib.pyplot as plt
a,b=0,10000 #Oil Dimension
alpha=1.833333 #Water flow
N=100 #number of points
h=(b-a)/N
X=np.linspace(a,b,N) #Generate N points between a and b
T=50 #Time
k=h/(3*alpha) 
r=k/h
S=np.zeros((len(X),T)) #saturation vector
S[0,:]=1 #boundary conditions

#Definitions
def fwater(s):
    return s**2
def foil(s):
    return ((1-s)**2)/4

def f(s):
    return alpha*fwater(s)/(fwater(s)+foil(s))
for n in range(0,T-1):
    for i in range(1,len(S)-1):
       S[i,n+1] = S[i,n] - r*(f(S[i,n])-f(S[i-1,n]))
K=41*(10**(-12)) #soil permeability
g=9.81 
S1=np.zeros((len(X),T)) 
S1[0,:]=1 #boundary conditions
rho_w,rho_o = 1000,892  #Density in kg/m^3
phi=0.40
beta = (rho_w - rho_o)*(K/phi)*g
def G(a,b):
    if -alpha+beta*fwater(a)<=0:
        G=(fwater(a)*(alpha+beta*foil(a)))/(fwater(a)+foil(a))
    else:
        G=(fwater(a)*(alpha+beta*foil(b)))/(fwater(a)+foil(b))
    return(G)
for n in range(0,T-1):
    for i in range(1,len(S)-1):
       a=S1[i,n]
       b=(k/h)
       c=G(S1[i,n],S1[i+1,n])-G(S1[i-1,n],S1[i,n])
       S1[i,n+1] = a - (b*c)
       
#Plots
plt.clf()
plt.plot(X,S[:,-1])
plt.plot(X,S1[:,-1])
plt.plot
plt.xlabel('lenght in meters')
plt.ylabel('Water saturation')
plt.legend(["without gravity",'with gravity'])
plt.title('Upwind method')
plt.show()

