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
S1=np.zeros((len(X),T)) #saturation vector 
S2=S1 #boundary conditions
S3=S1 
S1[0,:]=1 
S2[0,:]=S1[0,:]
S3[0,:]=S1[0,:]
U=np.zeros(N) #Godunov
U[0]=1 #boundary conditions

#Definitions
def fw(s):
    return s**2
def fo(s):
    return ((1-s)**2)/4
def f(s):
    return alpha*fw(s)/(fw(s)+fo(s))

#values of saturation :
#Upwind :
for n in range(0,T-1):
    for i in range(1,len(S1)-1):
       S1[i,n+1] = S1[i,n] - r*(f(S1[i,n])-f(S1[i-1,n]))

#Lax-Friederich:
for n in range(0,T-1):
    for i in range(1,len(S2)-1):
       S2[i,n+1] = (0.5*(S2[i+1,n]+S2[i-1,n]))-(0.5*r*(f(S2[i+1,n])-f(S2[i-1,n])))

#Godunov:
for n in range(0,T-1):
    for i in range(1,len(S3)-1):
        U[i]=S3[i,n]-r*(f(S3[i,n])-f(S3[i-1,n]))
        S3[i,n+1] = 0.5*(S3[i,n]+U[i]-r*(f(U[i])-f(U[i-1])))

# Plotting  
plt.clf()
plt.plot(X,S1[:,-1],label='Upwind method')
plt.plot(X,S2[:,-1],label='Lax-Friederich method')
plt.plot(X,S3[:,-1],label='Godunov method')
plt.xlabel('length in meter')
plt.ylabel('water saturation')
plt.title('Comparison of differents methods for oil extraction without gravity')
plt.legend()
plt.show()