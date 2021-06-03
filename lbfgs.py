#Libraries
import math
import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal
import scipy
from scipy.optimize import minimize

import periodic
import calculategrad
import initialise
import update
import energy as en

zcomp = 1
Lz=101
Lt=400
gconv = 1e-6
fconv = 1e-9
iters = 5e6
funciters = 1000
lsiters = 1000
splitt = 4.0
splitz = 6.0
timesplit = int(Lt / splitt)
timesplitf = int((splitt-1)*(Lt / splitt))
lowz = int(Lz/splitz)
highz = int((splitz-1))*int(Lz/splitz)
dz = 1.0
coza,cozb = 0,0
#time energy
sig = 10.0
#Elastic constants
ks = 1e-12
kt = 1e-12

#Bulk constants
a=3.0
b=2.0
c=1.0
s = (b + math.sqrt(b**2 + 24*a*c))/(4.0*c)
print(s)

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

(alpha,beta,gamma,Q3,Q5) = initialise.initialise(splitt,splitz,lowz,highz,timesplit,timesplitf,alpha,beta,gamma,Q3,Q5,Lz,Lt,w0,w1,a,b,c)
if zcomp == 1:
    (alpha,beta,gamma,Q3,Q5) = initialise.redefine(splitt,splitz,lowz,highz,timesplit,timesplitf,alpha,beta,gamma,Q3,Q5,Lz,Lt,w0,w1,a,b,c)


Q1 = np.add(alpha,beta)
Q4 = np.subtract(alpha,beta)
#Q5 = -Q1-Q4
Q2 = gamma

file = open("initialguess.dat", 'w')
gradfile =  open("graden.dat", 'w')
for t in range(0,Lt):
    for i in range(0,Lz):
        file.write("%f\n" % Q1[i,t])
        file.write("%f\n" % Q2[i,t])
        file.write("%f\n" % Q3[i,t])
        file.write("%f\n" % Q4[i,t])
        file.write("%f\n" % Q5[i,t])
        gradfile.write("0.0\n")
        gradfile.write("0.0\n")
        gradfile.write("0.0\n")
        gradfile.write("0.0\n")
        gradfile.write("0.0\n")

file.close()
gradfile.close()
GradE = np.loadtxt("graden.dat")
guess = np.loadtxt("initialguess.dat")
original = guess

print(np.shape(guess), type(guess))
z = 0
t = 0
GradE = calculategrad.calcgrad(guess,original,GradE,sig,Lz,Lt,ks,kt,q0,z,t,s,alpha,beta,gamma,a,b,c,Q1,Q2,Q3,Q4,Q5)

Q3 = np.zeros([Lz,Lt])
Q5 = np.zeros([Lz,Lt])
Q1 = np.zeros([Lz,Lt])
Q2 = np.zeros([Lz,Lt])
Q4 = np.zeros([Lz,Lt])
#minen = scipy.optimize.minimize(en.calcenergy,guess,args=((original,sig,Lz,Lt,ks,kt,q0,z,t,s,alpha,beta,gamma,a,b,c)),options={ 'ftol': conv, 'maxiter': iters,'disp': True}, method='nelder-mead')
minen = scipy.optimize.minimize(en.calcenergy,guess,\
args=((original,GradE,sig,Lz,Lt,ks,kt,q0,z,t,s,alpha,beta,gamma,a,b,c,Q1,Q2,Q3,Q4,Q5)),\
#options={ 'ftol': fconv, 'gtol': gconv, 'maxls': lsiters, 'maxfun': funciters 'maxiter': iters,'disp': True},\
options={'disp': True, 'maxls': 1e7, 'gtol': 1e-05, 'maxiter': 35000, 'ftol': 2.220446049250313e-09, 'maxcor': 10, 'maxfun': 1e6},\
method='L-BFGS-B',jac=True)

np.savetxt("energyarray.dat", minen.x)
Energy = np.zeros((Lt))
s = (b + math.sqrt(b**2 + 24*a*c))/(4.0*c)
qt = np.array([[2.0*s/3.0,0,0],[0,-s/3.0,0],[0,0,-s/3.0]])
Qt1,Qt2,Qt3,Qt4,Qt5= (2.0*s/3.0),0.0,0.0,-s/3.0,0.0
np.savetxt("differencearray.dat", guess-minen.x)
