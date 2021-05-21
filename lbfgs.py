#Libraries
import math
import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal
import scipy
from scipy.optimize import minimize

import periodic
import initialise
import update
import energy as en

Lz=30
Lt=30
conv = 1e-6
iters = 1e7
splitt = 4.0
splitz = 5.0
timesplit = int(Lt / splitt)
timesplitf = int((splitt-1)*(Lt / splitt))
lowz = int(Lz/splitz)
highz = int((splitz-1))*int(Lz/splitz)
dz = 1.0
coza,cozb = 0,0
#time energy
sig = 4.0

#Bulk constants
a=300.0
b=200.0
c=100.0
s = (b + math.sqrt(b**2 + 24*a*c))/(4.0*c)
print(s)
#Elastic constants

ks = 1e-12
kt = 1e-12
#initial and final windings for initial conditions

w0 = 1.0
w1 = 3.0
#Chirality
wideal=(w0+w1)/2.0
q0 = wideal*math.pi/(dz*(Lz-1))

fq5 = np.zeros([Lz,Lt])
fq3 = np.zeros([Lz,Lt])
fgamma = np.zeros([Lz,Lt])
falpha = np.zeros([Lz,Lt])
fbeta = np.zeros([Lz,Lt])
alpha = np.zeros([Lz,Lt])
beta = np.zeros([Lz,Lt])
gamma = np.zeros([Lz,Lt])
alphatemp = np.zeros([Lz,Lt])
betatemp = np.zeros([Lz,Lt])
gammatemp = np.zeros([Lz,Lt])
Q3 = np.zeros([Lz,Lt])
Q5 = np.zeros([Lz,Lt])
Q3temp = np.zeros([Lz,Lt])
Q5temp = np.zeros([Lz,Lt])
Q1 = np.zeros([Lz,Lt])
Q2 = np.zeros([Lz,Lt])
Q4 = np.zeros([Lz,Lt])

(alpha,beta,gamma,Q3,Q5,Lz,Lt,w0,w1,a,b,c) = initialise.initialise(splitt,alpha,beta,gamma,Q3,Q5,Lz,Lt,w0,w1,a,b,c)
Q1 = alpha + beta
Q4 = alpha - beta
#Q5 = -Q1-Q4
Q2 = gamma
file = open("initialguess.dat", 'w')
for t in range(0,Lt):
#    for x in range(0,10):
    for i in range(0,Lz):
        file.write("%f\n" % Q1[i,t])
        file.write("%f\n" % Q2[i,t])
        file.write("%f\n" % Q3[i,t])
        file.write("%f\n" % Q4[i,t])
        file.write("%f\n" % Q5[i,t])
        # file.write("1.0\n")
        # file.write("1.0\n")
        # file.write("1.0\n")
        # file.write("1.0\n")
        # file.write("1.0\n")


file.close()

guess = np.loadtxt("initialguess.dat")
print(np.shape(guess), type(guess))
z = 0
t = 0
minen = scipy.optimize.minimize(en.calcenergy,guess,args=((sig,Lz,Lt,ks,kt,q0,z,t,s,alpha,beta,gamma,a,b,c)),options={ 'ftol': conv, 'maxiter': iters,'disp': True}, method='nelder-mead')
#minen = scipy.optimize.minimize(en.calcenergy,guess,args=((sig,Lz,Lt,ks,kt,q0,z,t,s,alpha,beta,gamma,a,b,c)),options={ 'ftol': conv, 'maxiter': iters,'disp': True}, method='L-BFGS-B')
print(minen)
print(np.shape(minen.x))
np.savetxt("energyarray.dat", minen.x)
#minimize(fun(emptyarray, tuple), guess, method)

Energy = np.zeros((Lt))
s = (b + math.sqrt(b**2 + 24*a*c))/(4.0*c)
qt = np.array([[2.0*s/3.0,0,0],[0,-s/3.0,0],[0,0,-s/3.0]])
Qt1,Qt2,Qt3,Qt4,Qt5= (2.0*s/3.0),0.0,0.0,-s/3.0,0.0
np.savetxt("differencearray.dat", guess-minen.x)
#
# for t in range(0,Lt):
#     ener = np.zeros((Lz))
#     for z in range(0,Lz):
#         ener[z] = en.calcenergy(minen.x,sig,Lz,Lt,ks,kt,q0,z,t,s,alpha,beta,gamma,a,b,c)
#         #print(en[z])
#     Energy[t] = np.sum(ener)
#     print(np.sum(ener))
# print("finish")
# # np.savetxt('energy.dat',Energy)
# plt.plot(Energy)
# plt.show()
