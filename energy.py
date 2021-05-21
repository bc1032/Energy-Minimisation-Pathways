#Libraries
import math
import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal
import scipy
from scipy.optimize import minimize

def calcenergy(guess,sig,Lz,Lt,ks,kt,q0,z,t,s,alpha,beta,gamma,a,b,c):
    #print("Energy Called")
    gradzQ1,gradzQ2,gradzQ3,gradzQ4,gradzQ5=np.zeros((Lz)),np.zeros((Lz)),np.zeros((Lz)),np.zeros((Lz)),np.zeros((Lz))
    ws = 1
    dz = 1
    dt = 1
    Q3 = np.zeros([Lz,Lt])
    Q5 = np.zeros([Lz,Lt])
    Q1 = np.zeros([Lz,Lt])
    Q2 = np.zeros([Lz,Lt])
    Q4 = np.zeros([Lz,Lt])
    original = guess
    for t in range(0,Lt):
        for z in range(0,Lz):
            Q1[z,t] = guess[5*t*Lz + 5*z]
            Q2[z,t] = guess[5*t*Lz + 5*z + 1]
            Q3[z,t] = guess[5*t*Lz + 5*z + 2]
            Q4[z,t] = guess[5*t*Lz + 5*z + 3]
            Q5[z,t] = guess[5*t*Lz + 5*z + 4]

    for t in range(0,Lt):

        gradxQ1,gradxQ2,gradxQ3,gradxQ4,gradxQ5=np.zeros((Lz)),np.zeros((Lz)),np.zeros((Lz)),np.zeros((Lz)),np.zeros((Lz))
        splay,twist,bend,surface,bulk = 0.0,0.0,0.0,0.0,0.0
        Qt1,Qt2,Qt3,Qt4,Qt5= (2.0*s/3.0),0.0,0.0,-s/3.0,0.0
        gradxQ1 = np.zeros((Lz))
        gradxQ2 = np.zeros((Lz))
        gradxQ3 = np.zeros((Lz))
        gradxQ4 = np.zeros((Lz))
        gradxQ5 = np.zeros((Lz))
        for z in range(0,Lz):
            if t == 0:
                Q1[z,t] = original[5*z]
                Q2[z,t] = original[5*z + 1]
                Q3[z,t] = original[5*z + 2]
                Q4[z,t] = original[5*z + 3]
                Q5[z,t] = original[5*z + 4]
            elif t == Lt-1:
                Q1[z,t] = original[5*t*Lz + 5*z]
                Q2[z,t] = original[5*t*Lz + 5*z + 1]
                Q3[z,t] = original[5*t*Lz + 5*z + 2]
                Q4[z,t] = original[5*t*Lz + 5*z + 3]
                Q5[z,t] = original[5*t*Lz + 5*z + 4]

            if z == 0:
                Q1[z,t] = original[5*t*Lz + 5*z]
                Q2[z,t] = original[5*t*Lz + 5*z + 1]
                Q3[z,t] = original[5*t*Lz + 5*z + 2]
                Q4[z,t] = original[5*t*Lz + 5*z + 3]
                Q5[z,t] = original[5*t*Lz + 5*z + 4]

            elif z == Lz-1:
                Q1[z,t] = original[5*t*Lz + 5*z]
                Q2[z,t] = original[5*t*Lz + 5*z + 1]
                Q3[z,t] = original[5*t*Lz + 5*z + 2]
                Q4[z,t] = original[5*t*Lz + 5*z + 3]
                Q5[z,t] = original[5*t*Lz + 5*z + 4]

            time = 0.0
            if t != 0 and t!= Lt-1:
                time = (sig/2.0)*((-Q3[z,t-1] + Q3[z,t+1])**2/(2.*dt**2) + (-Q5[z,t-1] + Q5[z,t+1])**2/(2.*dt**2)\
                 + (-alpha[z,t-1] + alpha[z,t+1])**2/dt**2 +\
                 ((-alpha[z,t-1] + alpha[z,t+1])/(2.*dt) - (-beta[z,t-1] + beta[z,t+1])/(2.*dt))**2\
                 + ((-alpha[z,t-1] + alpha[z,t+1])/(2.*dt) + (-beta[z,t-1] + beta[z,t+1])/(2.*dt))**2 +\
                 (-gamma[z,t-1] + gamma[z,t+1])**2/(2.*dt**2))

            if z == Lz-1:
                gradzQ1[z] =  (ws) * (3*Q1[z,t]-4*Q1[z-1,t]+Q1[z-2,t]) / (2.0*dz)
                gradzQ2[z] =  (ws) * (3*Q2[z,t]-4*Q2[z-1,t]+Q2[z-2,t]) / (2.0*dz)
                gradzQ3[z] =  (ws) * (3*Q3[z,t]-4*Q3[z-1,t]+Q3[z-2,t]) / (2.0*dz)
                gradzQ4[z] =  (ws) * (3*Q4[z,t]-4*Q4[z-1,t]+Q4[z-2,t]) / (2.0*dz)
                gradzQ5[z] =  (ws) * (3*Q5[z,t]-4*Q5[z-1,t]+Q5[z-2,t]) / (2.0*dz)

                surface = 0.5*ws*( (Q1[z,t]-Qt1)**2 + 2*(Q2[z,t]-Qt2)**2 \
                + 2*(Q3[z,t]-Qt3)**2+(Qt1 + Qt4 - Q1[z,t] - Q4[z,t])**2 \
                + (Q4[z,t] - Qt4)**2 +2*(Q5[z,t] - Qt5)**2 )

            elif z == 0:
                gradzQ1[z] =  (ws) * -(3*Q1[z,t]-4*Q1[z-1,t]+Q1[z-2,t]) / (2.0*dz)
                gradzQ2[z] =  (ws) * -(3*Q2[z,t]-4*Q2[z-1,t]+Q2[z-2,t]) / (2.0*dz)
                gradzQ3[z] =  (ws) * -(3*Q3[z,t]-4*Q3[z-1,t]+Q3[z-2,t]) / (2.0*dz)
                gradzQ4[z] =  (ws) * -(3*Q4[z,t]-4*Q4[z-1,t]+Q4[z-2,t]) / (2.0*dz)
                gradzQ5[z] =  (ws) * -(3*Q5[z,t]-4*Q5[z-1,t]+Q5[z-2,t]) / (2.0*dz)

                surface = 0.5*ws*( (Q1[z,t]-Qt1)**2 + 2*(Q2[z,t]-Qt2)**2 \
                + 2*(Q3[z,t]-Qt3)**2+(Qt1 + Qt4 - Q1[z,t] - Q4[z,t])**2 \
                + (Q4[z,t] - Qt4)**2 +2*(Q5[z,t] - Qt5)**2 )
            else:
                gradzQ1[z] = ((2.0/3.0)*(Q1[z+1,t]-Q1[z-1,t])+(1.0/12.0)*(Q1[z+1,t]-Q1[z-1,t])) / dz
                gradzQ2[z] = ((2.0/3.0)*(Q2[z+1,t]-Q2[z-1,t])+(1.0/12.0)*(Q2[z+1,t]-Q2[z-1,t])) / dz
                gradzQ3[z] = ((2.0/3.0)*(Q3[z+1,t]-Q3[z-1,t])+(1.0/12.0)*(Q3[z+1,t]-Q3[z-1,t])) / dz
                gradzQ4[z] = ((2.0/3.0)*(Q4[z+1,t]-Q4[z-1,t])+(1.0/12.0)*(Q4[z+1,t]-Q4[z-1,t])) / dz
                gradzQ5[z] = ((2.0/3.0)*(Q5[z+1,t]-Q5[z-1,t])+(1.0/12.0)*(Q5[z+1,t]-Q5[z-1,t])) / dz

            bulk += (-a*(Q1[z,t]*Q4[z,t] + Q1[z,t]**2 + Q2[z,t]**2 + Q3[z,t]**2 + Q4[z,t]**2 + Q5[z,t]**2) \
                    + b*(Q1[z,t]*(-Q2[z,t]**2 + Q4[z,t]**2 + Q5[z,t]**2)\
                    + (Q1[z,t]**2)*Q4[z,t] - 2*Q2[z,t]*Q3[z,t]*Q5[z,t] \
                    - (Q2[z,t]**2)*Q4[z,t] + (Q3[z,t]**2)*Q4[z,t]) \
                    + c*(Q1[z,t]*Q4[z,t] + Q1[z,t]**2 + Q2[z,t]**2 + Q3[z,t]**2 + Q4[z,t]**2 + Q5[z,t]**2)**2 )

            twist += 0.5*kt*( (-2*q0*Q1[z,t]+gradzQ2[z])**2 \
                    + 2*(2*q0*Q3[z,t] - gradzQ5[z])*(2*q0*Q3[z,t]+gradxQ2[z]) \
                    + 2*(2*q0*Q2[z,t]-gradzQ4[z])*(2*q0*Q2[z,t]+gradzQ1[z] - gradxQ3[z]) \
                    + 2*(2*q0*Q5[z,t]+gradxQ4[z])*(2*q0*Q5[z,t]+gradzQ3[z]+gradxQ1[z]-gradxQ4[z]) \
                    + (2*q0*Q4[z,t]+gradzQ2[z]-gradxQ5[z])**2 + (-2*q0*(Q1[z,t]+Q4[z,t]) + gradxQ5[z])**2)

            splay += 0.5*ks*((gradzQ3[z]+gradxQ1[z])**2 \
                    + ( gradzQ5[z] + gradxQ2[z])**2 + (gradzQ1[z] + gradzQ4[z] - gradxQ3[z])**2 )

    energy = bulk + splay + twist + surface + time
            #print(bulk,splay,twist,surface)

        #print(energy)
    return(energy)
