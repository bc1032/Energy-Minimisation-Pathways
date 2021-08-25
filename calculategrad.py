import math
import numpy as np

def calcgrad(guess,original,GradE,sig,Lz,Lt,ks,kt,q0,z,t,s,alpha,beta,gamma,a,b,c,Q1,Q2,Q3,Q4,Q5,ws,timeen,splay,twist,bend,surface,bulk):
    dz = 1
    dt = 1
    scale = 1.0#/((Lz*Lt))
    Qt1,Qt2,Qt3,Qt4,Qt5 = (2.0*s/3.0),0.0,0.0,-s/3.0,0.0

    #Calculate Bulk Gradient Energy Term and Elastic Energy Terms
    for t in range(2, Lt-2):
        for z in range(2,Lz-2):

            GradE[5*t*Lz + 5*z] = scale*-(dt**2*(ks + kt)*Q1[z-2,t] + 2*dz**2*sig*Q1[z,t-2] - 2*(dt**2*(ks + kt) + 2*dz**2*sig)*Q1[z,t] +\
                2*dz**2*sig*Q1[z,t+2] + dt**2*(ks + kt)*Q1[z+2,t] - 8*dt**2*dz*kt*q0*Q2[z-1,t] +\
                8*dt**2*dz*kt*q0*Q2[z+1,t] + dt**2*ks*Q4[z-2,t] + dz**2*sig*Q4[z,t-2] - 2*dt**2*ks*Q4[z,t] -\
                2*dz**2*sig*Q4[z,t] + dz**2*sig*Q4[z,t+2] + dt**2*ks*Q4[z+2,t])/(4.*dt*dz)

            GradE[5*t*Lz + 5*z + 1] = scale*-2*dt*kt*q0*Q1[z-1,t] + 2*dt*kt*q0*Q1[z+1,t] -\
                (dt**2*kt*Q2[z-2,t] + dt**2*kt*(-2*Q2[z,t] + Q2[z+2,t]) +\
                dz*(dz*sig*(Q2[z,t-2] - 2*Q2[z,t] + Q2[z,t+2]) + 4*dt**2*kt*q0*(-Q4[z-1,t] + Q4[z+1,t])))/(2.*dt*dz)

            GradE[5*t*Lz + 5*z + 2] = scale*(-(dt**2*(ks + kt)*Q3[z-2,t]) + dt**2*(ks + kt)*(2*Q3[z,t] - Q3[z+2,t]) -\
                2*dz*(dz*sig*(Q3[z,t-2] - 2*Q3[z,t] + Q3[z,t+2]) + 4*dt**2*kt*q0*(-Q5[z-1,t] + Q5[z+1,t])))/(4.*dt*dz)

            GradE[5*t*Lz + 5*z + 3] = scale*-(dt**2*ks*Q1[z-2,t] + dz**2*sig*Q1[z,t-2] - 2*(dt**2*ks + dz**2*sig)*Q1[z,t] + dz**2*sig*Q1[z,t+2] +\
                dt**2*ks*Q1[z+2,t] + 8*dt**2*dz*kt*q0*Q2[z-1,t] - 8*dt**2*dz*kt*q0*Q2[z+1,t] + dt**2*ks*Q4[z-2,t] +\
                dt**2*kt*Q4[z-2,t] + 2*dz**2*sig*Q4[z,t-2] - 2*dt**2*ks*Q4[z,t] - 2*dt**2*kt*Q4[z,t] - 4*dz**2*sig*Q4[z,t] +\
                2*dz**2*sig*Q4[z,t+2] + dt**2*(ks + kt)*Q4[z+2,t])/(4.*dt*dz)

            GradE[5*t*Lz + 5*z + 4] = scale*-(8*dt**2*dz*kt*q0*Q3[z-1,t] - 8*dt**2*dz*kt*q0*Q3[z+1,t] + dt**2*ks*Q5[z-2,t] + dt**2*kt*Q5[z-2,t] +\
                2*dz**2*sig*Q5[z,t-2] - 2*dt**2*ks*Q5[z,t] - 2*dt**2*kt*Q5[z,t] - 4*dz**2*sig*Q5[z,t] + 2*dz**2*sig*Q5[z,t+2] +\
                dt**2*(ks + kt)*Q5[z+2,t])/(4.*dt*dz)

            GradE[5*t*Lz + 5*z] += scale*-(b*Q2[z,t]**2) + (2*Q1[z,t] + Q4[z,t])*\
                (-a + b*Q4[z,t] + 2*c*(Q1[z,t]**2 + Q2[z,t]**2 + Q3[z,t]**2 + Q1[z,t]*Q4[z,t] + Q4[z,t]**2)) +\
                (b + 4*c*Q1[z,t] + 2*c*Q4[z,t])*Q5[z,t]**2

            GradE[5*t*Lz + 5*z + 1] += scale*4*c*Q2[z,t]**3 - 2*b*Q3[z,t]*Q5[z,t] - 2*Q2[z,t]*\
                (a - 2*c*Q1[z,t]**2 + b*Q4[z,t] + Q1[z,t]*(b - 2*c*Q4[z,t]) - 2*c*(Q3[z,t]**2 + Q4[z,t]**2 + Q5[z,t]**2))

            GradE[5*t*Lz + 5*z + 2] +=  scale*4*c*Q3[z,t]**3 - 2*b*Q2[z,t]*Q5[z,t] + Q3[z,t]*\
                (-2*a + 2*b*Q4[z,t] + 4*c*(Q1[z,t]**2 + Q2[z,t]**2 + Q1[z,t]*Q4[z,t] + Q4[z,t]**2 + Q5[z,t]**2))

            GradE[5*t*Lz + 5*z + 3] +=  scale*-(a*(Q1[z,t] + 2*Q4[z,t])) + b*(Q1[z,t]**2 - Q2[z,t]**2 + Q3[z,t]**2 + 2*Q1[z,t]*Q4[z,t]) +\
                2*c*(Q1[z,t] + 2*Q4[z,t])*(Q1[z,t]**2 + Q2[z,t]**2 + Q3[z,t]**2 + Q1[z,t]*Q4[z,t] + Q4[z,t]**2 + Q5[z,t]**2)

            GradE[5*t*Lz + 5*z + 4] += scale*-2*b*Q2[z,t]*Q3[z,t] + 2*(-a + 2*c*Q1[z,t]**2 + Q1[z,t]*(b + 2*c*Q4[z,t]) +\
                2*c*(Q2[z,t]**2 + Q3[z,t]**2 + Q4[z,t]**2))*Q5[z,t] + 4*c*Q5[z,t]**3

    #GradE = GradE
    return(GradE)
