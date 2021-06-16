#Libraries
import math
import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal
import scipy
from scipy.optimize import minimize
import calculategrad

def calcenergy(guess,original,GradE,sig,Lz,Lt,ks,kt,q0,z,t,s,alpha,beta,gamma,a,b,c,Q1,Q2,Q3,Q4,Q5,ws,timeen,splay,twist,bend,surface,bulk):
    #print("Energy Called")
    gradzQ1,gradzQ2,gradzQ3,gradzQ4,gradzQ5=np.zeros((Lz)),np.zeros((Lz)),np.zeros((Lz)),np.zeros((Lz)),np.zeros((Lz))
    gradxQ1 = np.zeros((Lz))
    gradxQ2 = np.zeros((Lz))
    gradxQ3 = np.zeros((Lz))
    gradxQ4 = np.zeros((Lz))
    gradxQ5 = np.zeros((Lz))
    dz = 1
    dt = 1
    Qt1,Qt2,Qt3,Qt4,Qt5= (2.0*s/3.0),0.0,0.0,-s/3.0,0.0

    for t in range(1,Lt-1):
        for z in range(0,Lz):
            Q1[z,t] = guess[5*t*Lz + 5*z]
            Q2[z,t] = guess[5*t*Lz + 5*z + 1]
            Q3[z,t] = guess[5*t*Lz + 5*z + 2]
            Q4[z,t] = guess[5*t*Lz + 5*z + 3]
            Q5[z,t] = guess[5*t*Lz + 5*z + 4]

    for t in range(1,Lt-1):

        for z in range(0,Lz):
            if t == 0:
                Q1[z,t] = original[5*z]
                Q2[z,t] = original[5*z + 1]
                Q3[z,t] = original[5*z + 2]
                Q4[z,t] = original[5*z + 3]
                Q5[z,t] = original[5*z + 4]
            if t == Lt-1:
                Q1[z,t] = original[5*t*Lz + 5*z]
                Q2[z,t] = original[5*t*Lz + 5*z + 1]
                Q3[z,t] = original[5*t*Lz + 5*z + 2]
                Q4[z,t] = original[5*t*Lz + 5*z + 3]
                Q5[z,t] = original[5*t*Lz + 5*z + 4]

            if t != 0 and t != Lt-1:
                timeen += (sig*((-Q1[z,t-1] + Q1[z,t+1])**2/(4.*dz**2) + (-Q2[z,t-1] + Q2[z,t+1])**2/(2.*dz**2)\
                        + (-Q3[z,t-1] + Q3[z,t+1])**2/(2.*dz**2) + \
                        (-Q4[z,t-1] + Q4[z,t+1])**2/(4.*dz**2) + \
                        (-(-Q1[z,t-1] + Q1[z,t+1])/(2.*dz) - (-Q4[z,t-1] + Q4[z,t+1])/(2.*dz))**2 +\
                        (-Q5[z,t-1] + Q5[z,t+1])**2/(2.*dz**2)))/2.

            if z == 0:
                Q1[z,t] = original[5*t*Lz + 5*z]
                Q2[z,t] = original[5*t*Lz + 5*z + 1]
                Q3[z,t] = original[5*t*Lz + 5*z + 2]
                Q4[z,t] = original[5*t*Lz + 5*z + 3]
                Q5[z,t] = original[5*t*Lz + 5*z + 4]

                gradzQ1[z] =  (ws) * (3*Q1[z,t]-4*Q1[z+1,t]+Q1[z+2,t]) / (2.0*dz)
                gradzQ2[z] =  (ws) * (3*Q2[z,t]-4*Q2[z+1,t]+Q2[z+2,t]) / (2.0*dz)
                gradzQ3[z] =  (ws) * (3*Q3[z,t]-4*Q3[z+1,t]+Q3[z+2,t]) / (2.0*dz)
                gradzQ4[z] =  (ws) * (3*Q4[z,t]-4*Q4[z+1,t]+Q4[z+2,t]) / (2.0*dz)
                gradzQ5[z] =  (ws) * (3*Q5[z,t]-4*Q5[z+1,t]+Q5[z+2,t]) / (2.0*dz)

                surface += 0.5*ws*((Q1[z,t]-Qt1)**2 + 2*(Q2[z,t]-Qt2)**2 \
                + 2*(Q3[z,t]-Qt3)**2+(Qt1 + Qt4 - Q1[z,t] - Q4[z,t])**2 \
                + (Q4[z,t] - Qt4)**2 +2*(Q5[z,t] - Qt5)**2 )

                bulk += -(a*(Q1[z,t]**2 + Q2[z,t]**2 + Q3[z,t]**2 + Q1[z,t]*Q4[z,t] + Q4[z,t]**2 + Q5[z,t]**2)) +\
                        c*(Q1[z,t]**2 + Q2[z,t]**2 + Q3[z,t]**2 + Q1[z,t]*Q4[z,t] + Q4[z,t]**2 + Q5[z,t]**2)**2 +\
                        b*(Q1[z,t]**2*Q4[z,t] + (-Q2[z,t]**2 + Q3[z,t]**2)*Q4[z,t] - 2*Q2[z,t]*Q3[z,t]*Q5[z,t] +\
                        Q1[z,t]*(-Q2[z,t]**2 + Q4[z,t]**2 + Q5[z,t]**2))
            if z == Lz-1:
                Q1[z,t] = original[5*t*Lz + 5*z]
                Q2[z,t] = original[5*t*Lz + 5*z + 1]
                Q3[z,t] = original[5*t*Lz + 5*z + 2]
                Q4[z,t] = original[5*t*Lz + 5*z + 3]
                Q5[z,t] = original[5*t*Lz + 5*z + 4]
                gradzQ1[z] =  (ws) * -(3*Q1[z,t]-4*Q1[z-1,t]+Q1[z-2,t]) / (2.0*dz)
                gradzQ2[z] =  (ws) * -(3*Q2[z,t]-4*Q2[z-1,t]+Q2[z-2,t]) / (2.0*dz)
                gradzQ3[z] =  (ws) * -(3*Q3[z,t]-4*Q3[z-1,t]+Q3[z-2,t]) / (2.0*dz)
                gradzQ4[z] =  (ws) * -(3*Q4[z,t]-4*Q4[z-1,t]+Q4[z-2,t]) / (2.0*dz)
                gradzQ5[z] =  (ws) * -(3*Q5[z,t]-4*Q5[z-1,t]+Q5[z-2,t]) / (2.0*dz)

                surface += 0.5*ws*( (Q1[z,t]-Qt1)**2 + 2*(Q2[z,t]-Qt2)**2 \
                + 2*(Q3[z,t]-Qt3)**2+(Qt1 + Qt4 - Q1[z,t] - Q4[z,t])**2 \
                + (Q4[z,t] - Qt4)**2 +2*(Q5[z,t] - Qt5)**2 )

                bulk += -(a*(Q1[z,t]**2 + Q2[z,t]**2 + Q3[z,t]**2 + Q1[z,t]*Q4[z,t] + Q4[z,t]**2 + Q5[z,t]**2)) +\
                        c*(Q1[z,t]**2 + Q2[z,t]**2 + Q3[z,t]**2 + Q1[z,t]*Q4[z,t] + Q4[z,t]**2 + Q5[z,t]**2)**2 +\
                        b*(Q1[z,t]**2*Q4[z,t] + (-Q2[z,t]**2 + Q3[z,t]**2)*Q4[z,t] - 2*Q2[z,t]*Q3[z,t]*Q5[z,t] +\
                        Q1[z,t]*(-Q2[z,t]**2 + Q4[z,t]**2 + Q5[z,t]**2))

            else:

                gradzQ1[z] = ((2.0/3.0)*(Q1[z+1,t]-Q1[z-1,t])+(1.0/12.0)*(Q1[z+1,t]-Q1[z-1,t])) / dz
                gradzQ2[z] = ((2.0/3.0)*(Q2[z+1,t]-Q2[z-1,t])+(1.0/12.0)*(Q2[z+1,t]-Q2[z-1,t])) / dz
                gradzQ3[z] = ((2.0/3.0)*(Q3[z+1,t]-Q3[z-1,t])+(1.0/12.0)*(Q3[z+1,t]-Q3[z-1,t])) / dz
                gradzQ4[z] = ((2.0/3.0)*(Q4[z+1,t]-Q4[z-1,t])+(1.0/12.0)*(Q4[z+1,t]-Q4[z-1,t])) / dz
                gradzQ5[z] = ((2.0/3.0)*(Q5[z+1,t]-Q5[z-1,t])+(1.0/12.0)*(Q5[z+1,t]-Q5[z-1,t])) / dz


                bulk += -(a*(Q1[z,t]**2 + Q2[z,t]**2 + Q3[z,t]**2 + Q1[z,t]*Q4[z,t] + Q4[z,t]**2 + Q5[z,t]**2)) +\
                        c*(Q1[z,t]**2 + Q2[z,t]**2 + Q3[z,t]**2 + Q1[z,t]*Q4[z,t] + Q4[z,t]**2 + Q5[z,t]**2)**2 +\
                        b*(Q1[z,t]**2*Q4[z,t] + (-Q2[z,t]**2 + Q3[z,t]**2)*Q4[z,t] - 2*Q2[z,t]*Q3[z,t]*Q5[z,t] +\
                        Q1[z,t]*(-Q2[z,t]**2 + Q4[z,t]**2 + Q5[z,t]**2))


                twist += (kt*((-Q1[z-1,t] + Q1[z+1,t])**2/(4.*dt**2) - \
                    (2*q0*Q1[z,t]*(-Q2[z-1,t] + Q2[z+1,t]))/dt + \
                    (-Q2[z-1,t] + Q2[z+1,t])**2/(2.*dt**2) + \
                    (-Q3[z-1,t] + Q3[z+1,t])**2/(4.*dt**2) + \
                    (2*q0*(-Q2[z-1,t] + Q2[z+1,t])*Q4[z,t])/dt + \
                    (-Q4[z-1,t] + Q4[z+1,t])**2/(4.*dt**2) + \
                    4*q0*Q2[z,t]*((-Q1[z-1,t] + Q1[z+1,t])/(2.*dt) - \
                    (-Q4[z-1,t] + Q4[z+1,t])/(2.*dt)) + \
                    (2*q0*(-Q3[z-1,t] + Q3[z+1,t])*Q5[z,t])/dt - \
                    (2*q0*Q3[z,t]*(-Q5[z-1,t] + Q5[z+1,t]))/dt + \
                    (-Q5[z-1,t] + Q5[z+1,t])**2/(4.*dt**2)))/2.

                splay += (ks*((-Q3[z-1,t] + Q3[z+1,t])**2/(4.*dt**2) + \
                    (-(-Q1[z-1,t] + Q1[z+1,t])/(2.*dt) - (-Q4[z-1,t] + Q4[z+1,t])/(2.*dt))**2\
                    + (-Q5[z-1,t] + Q5[z+1,t])**2/(4.*dt**2)))/2.

    energy = (bulk + splay + twist + timeen + surface)# / (Lz*Lt)
    calculategrad.calcgrad(guess,original,GradE,sig,Lz,Lt,ks,kt,q0,z,t,s,alpha,beta,gamma,a,b,c,Q1,Q2,Q3,Q4,Q5,ws)
    #print(bulk,splay,twist,surface)

        #print(energy)
    return(energy,GradE)
    #return(energy)
