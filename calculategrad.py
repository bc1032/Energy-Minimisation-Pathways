import math
import numpy as np

def calcgrad(guess,original,GradE,sig,Lz,Lt,ks,kt,q0,z,t,s,alpha,beta,gamma,a,b,c,Q1,Q2,Q3,Q4,Q5,ws):
    dz = 1
    dt = 1
    scale = 1.0#/(Lz*Lt)
    Qt1,Qt2,Qt3,Qt4,Qt5 = (2.0*s/3.0),0.0,0.0,-s/3.0,0.0

    #Calculate Bulk Gradient Energy Term and Elastic Energy Terms
    for t in range(1, Lt-1):


        for z in range(0,Lz):

            if z == 0:

                GradE[5*t*Lz + 5*z] =  scale*(dt*(5*(ks + kt)*Q1[0,t-1] + 5*(ks + kt)*Q1[0,t] + 5*ks*Q1[0,t+1] + 5*kt*Q1[0,t+1] - 6*ks*Q1[1,t-1] -\
                    6*kt*Q1[1,t-1] - 6*ks*Q1[1,t] - 6*kt*Q1[1,t] - 6*ks*Q1[1,t+1] - 6*kt*Q1[1,t+1] + ks*Q1[2,t-1] +\
                    kt*Q1[2,t-1] + ks*Q1[2,t] + kt*Q1[2,t] + ks*Q1[2,t+1] + kt*Q1[2,t+1] - 10*dz*kt*q0*Q2[1,t-1] -\
                    10*dz*kt*q0*Q2[1,t] - 10*dz*kt*q0*Q2[1,t+1] + 2*dz*kt*q0*Q2[2,t-1] + 2*dz*kt*q0*Q2[2,t] +\
                    2*dz*kt*q0*Q2[2,t+1] + 5*ks*Q4[0,t-1] + 5*ks*Q4[0,t] + 5*ks*Q4[0,t+1] - 6*ks*Q4[1,t-1] -\
                    6*ks*Q4[1,t] - 6*ks*Q4[1,t+1] + ks*Q4[2,t-1] + ks*Q4[2,t] + ks*Q4[2,t+1]))/(2.*dz)\

                GradE[5*t*Lz + 5*z + 1] =  scale*(dt*kt*(5*Q2[0,t-1] + 5*Q2[0,t] + 5*Q2[0,t+1] - 6*Q2[1,t-1] - 6*Q2[1,t] - 6*Q2[1,t+1] +\
                    Q2[2,t-1] + Q2[2,t] + Q2[2,t+1] +\
                    dz*q0*(5*Q1[1,t-1] + 5*Q1[1,t] + 5*Q1[1,t+1] - Q1[2,t-1] - Q1[2,t] - Q1[2,t+1] -\
                    5*Q4[1,t-1] - 5*Q4[1,t] - 5*Q4[1,t+1] + Q4[2,t-1] + Q4[2,t] + Q4[2,t+1])))/dz

                GradE[5*t*Lz + 5*z + 2] =  scale*(dt*(5*(ks + kt)*Q3[0,t-1] + 5*(ks + kt)*Q3[0,t] + 5*ks*Q3[0,t+1] + 5*kt*Q3[0,t+1] - 6*ks*Q3[1,t-1] -\
                    6*kt*Q3[1,t-1] - 6*ks*Q3[1,t] - 6*kt*Q3[1,t] - 6*ks*Q3[1,t+1] - 6*kt*Q3[1,t+1] + ks*Q3[2,t-1] +\
                    kt*Q3[2,t-1] + ks*Q3[2,t] + kt*Q3[2,t] + ks*Q3[2,t+1] + kt*Q3[2,t+1] - 10*dz*kt*q0*Q5[1,t-1] -\
                    10*dz*kt*q0*Q5[1,t] - 10*dz*kt*q0*Q5[1,t+1] + 2*dz*kt*q0*Q5[2,t-1] + 2*dz*kt*q0*Q5[2,t] +\
                    2*dz*kt*q0*Q5[2,t+1]))/(2.*dz)



                GradE[5*t*Lz + 5*z + 3] = scale*(dt*(5*ks*Q1[0,t-1] + 5*ks*Q1[0,t] + 5*ks*Q1[0,t+1] - 6*ks*Q1[1,t-1] - 6*ks*Q1[1,t] -\
                    6*ks*Q1[1,t+1] + ks*Q1[2,t-1] + ks*Q1[2,t] + ks*Q1[2,t+1] +\
                    2*dz*kt*q0*(5*Q2[1,t-1] + 5*Q2[1,t] + 5*Q2[1,t+1] - Q2[2,t-1] - Q2[2,t] - Q2[2,t+1]) +\
                    5*ks*Q4[0,t-1] + 5*kt*Q4[0,t-1] + 5*ks*Q4[0,t] + 5*kt*Q4[0,t] + 5*ks*Q4[0,t+1] + 5*kt*Q4[0,t+1] -\
                    6*ks*Q4[1,t-1] - 6*kt*Q4[1,t-1] - 6*ks*Q4[1,t] - 6*kt*Q4[1,t] - 6*ks*Q4[1,t+1] - 6*kt*Q4[1,t+1] +\
                    ks*Q4[2,t-1] + kt*Q4[2,t-1] + ks*Q4[2,t] + kt*Q4[2,t] + (ks + kt)*Q4[2,t+1]))/(2.*dz)


                GradE[5*t*Lz + 5*z + 4] = scale*(dt*(ks*(5*Q5[0,t-1] + 5*Q5[0,t] + 5*Q5[0,t+1] - 6*Q5[1,t-1] - 6*Q5[1,t] - 6*Q5[1,t+1] +\
                    Q5[2,t-1] + Q5[2,t] + Q5[2,t+1]) +\
                    kt*(2*dz*q0*(5*Q3[1,t-1] + 5*Q3[1,t] + 5*Q3[1,t+1] - Q3[2,t-1] - Q3[2,t] - Q3[2,t+1]) +\
                    5*Q5[0,t-1] + 5*Q5[0,t] + 5*Q5[0,t+1] - 6*Q5[1,t-1] - 6*Q5[1,t] - 6*Q5[1,t+1] +\
                    Q5[2,t-1] + Q5[2,t] + Q5[2,t+1])))/(2.*dz)


                GradE[5*t*Lz + 5*z] += scale*(ws/2.0)*(2*(Q1[z,t]-Qt1)-2*(Qt1+Qt4-Q1[z,t]-Q4[z,t]))

                GradE[5*t*Lz + 5*z + 1] +=  scale*2.0*ws*(Q2[z,t]-Qt2)

                GradE[5*t*Lz + 5*z + 2] +=  scale*2.0*ws*(Q3[z,t]-Qt3)

                GradE[5*t*Lz + 5*z + 3] += scale*(ws/2.0)*(2*(Q4[z,t]-Qt4)-2*(Qt1+Qt4-Q1[z,t]-Q4[z,t]))

                GradE[5*t*Lz + 5*z + 4] += scale*2.0*ws*(Q5[z,t]-Qt5)


            elif z == 1:

                GradE[5*t*Lz + 5*z] =  scale*-(dt*(12*(ks + kt)*Q1[0,t-1] + 12*(ks + kt)*Q1[0,t] + 12*ks*Q1[0,t+1] + 12*kt*Q1[0,t+1] -\
                    17*ks*Q1[1,t-1] - 17*kt*Q1[1,t-1] - 17*ks*Q1[1,t] - 17*kt*Q1[1,t] - 17*ks*Q1[1,t+1] -\
                    17*kt*Q1[1,t+1] + 4*ks*Q1[2,t-1] + 4*kt*Q1[2,t-1] + 4*ks*Q1[2,t] + 4*kt*Q1[2,t] + 4*ks*Q1[2,t+1] +\
                    4*kt*Q1[2,t+1] + ks*Q1[3,t-1] + kt*Q1[3,t-1] + ks*Q1[3,t] + kt*Q1[3,t] + ks*Q1[3,t+1] +\
                    kt*Q1[3,t+1] - 20*dz*kt*q0*Q2[0,t-1] - 20*dz*kt*q0*Q2[0,t] - 20*dz*kt*q0*Q2[0,t+1] +\
                    8*dz*kt*q0*Q2[2,t-1] + 8*dz*kt*q0*Q2[2,t] + 8*dz*kt*q0*Q2[2,t+1] + 12*ks*Q4[0,t-1] + 12*ks*Q4[0,t] +\
                    12*ks*Q4[0,t+1] - 17*ks*Q4[1,t-1] - 17*ks*Q4[1,t] - 17*ks*Q4[1,t+1] + 4*ks*Q4[2,t-1] +\
                    4*ks*Q4[2,t] + 4*ks*Q4[2,t+1] + ks*Q4[3,t-1] + ks*Q4[3,t] + ks*Q4[3,t+1]))/(4.*dz)

                GradE[5*t*Lz + 5*z + 1] =  scale*-(dt*kt*(12*Q2[0,t-1] + 12*Q2[0,t] + 12*Q2[0,t+1] - 17*Q2[1,t-1] - 17*Q2[1,t] - 17*Q2[1,t+1] +\
                    4*Q2[2,t-1] + 4*Q2[2,t] + 4*Q2[2,t+1] + Q2[3,t-1] + Q2[3,t] + Q2[3,t+1] +\
                    2*dz*q0*(5*Q1[0,t-1] + 5*Q1[0,t] + 5*Q1[0,t+1] - 2*Q1[2,t-1] - 2*Q1[2,t] - 2*Q1[2,t+1] -\
                    5*Q4[0,t-1] - 5*Q4[0,t] - 5*Q4[0,t+1] + 2*(Q4[2,t-1] + Q4[2,t] + Q4[2,t+1]))))/(2.*dz)

                GradE[5*t*Lz + 5*z + 2] =  scale*-(dt*(12*(ks + kt)*Q3[0,t-1] + 12*(ks + kt)*Q3[0,t] + 12*ks*Q3[0,t+1] + 12*kt*Q3[0,t+1] -\
                    17*ks*Q3[1,t-1] - 17*kt*Q3[1,t-1] - 17*ks*Q3[1,t] - 17*kt*Q3[1,t] - 17*ks*Q3[1,t+1] -\
                    17*kt*Q3[1,t+1] + 4*ks*Q3[2,t-1] + 4*kt*Q3[2,t-1] + 4*ks*Q3[2,t] + 4*kt*Q3[2,t] + 4*ks*Q3[2,t+1] +\
                    4*kt*Q3[2,t+1] + ks*Q3[3,t-1] + kt*Q3[3,t-1] + ks*Q3[3,t] + kt*Q3[3,t] + ks*Q3[3,t+1] +\
                    kt*Q3[3,t+1] - 20*dz*kt*q0*Q5[0,t-1] - 20*dz*kt*q0*Q5[0,t] - 20*dz*kt*q0*Q5[0,t+1] +\
                    8*dz*kt*q0*Q5[2,t-1] + 8*dz*kt*q0*Q5[2,t] + 8*dz*kt*q0*Q5[2,t+1]))/(4.*dz)


                GradE[5*t*Lz + 5*z + 3] =  scale*-(dt*(12*ks*Q1[0,t-1] + 12*ks*Q1[0,t] + 12*ks*Q1[0,t+1] - 17*ks*Q1[1,t-1] - 17*ks*Q1[1,t] -\
                    17*ks*Q1[1,t+1] + 4*ks*Q1[2,t-1] + 4*ks*Q1[2,t] + 4*ks*Q1[2,t+1] + ks*Q1[3,t-1] + ks*Q1[3,t] +\
                    ks*Q1[3,t+1] + 4*dz*kt*q0*(5*Q2[0,t-1] + 5*Q2[0,t] + 5*Q2[0,t+1] -\
                    2*(Q2[2,t-1] + Q2[2,t] + Q2[2,t+1])) + 12*ks*Q4[0,t-1] + 12*kt*Q4[0,t-1] + 12*ks*Q4[0,t] +\
                    12*kt*Q4[0,t] + 12*ks*Q4[0,t+1] + 12*kt*Q4[0,t+1] - 17*ks*Q4[1,t-1] - 17*kt*Q4[1,t-1] -\
                    17*ks*Q4[1,t] - 17*kt*Q4[1,t] - 17*ks*Q4[1,t+1] - 17*kt*Q4[1,t+1] + 4*ks*Q4[2,t-1] +\
                    4*kt*Q4[2,t-1] + 4*ks*Q4[2,t] + 4*kt*Q4[2,t] + 4*ks*Q4[2,t+1] + 4*kt*Q4[2,t+1] + ks*Q4[3,t-1] +\
                    kt*Q4[3,t-1] + ks*Q4[3,t] + kt*Q4[3,t] + (ks + kt)*Q4[3,t+1]))/(4.*dz)

                GradE[5*t*Lz + 5*z + 4] =  scale*(dt*(ks*(12*Q5[0,t-1] + 12*Q5[0,t] + 12*Q5[0,t+1] - 17*Q5[1,t-1] - 17*Q5[1,t] - 17*Q5[1,t+1] +\
                    4*Q5[2,t-1] + 4*Q5[2,t] + 4*Q5[2,t+1] + Q5[3,t-1] + Q5[3,t] + Q5[3,t+1]) +\
                    kt*(4*dz*q0*(5*Q3[0,t-1] + 5*Q3[0,t] + 5*Q3[0,t+1] - 2*(Q3[2,t-1] + Q3[2,t] + Q3[2,t+1])) +\
                    12*Q5[0,t-1] + 12*Q5[0,t] + 12*Q5[0,t+1] - 17*Q5[1,t-1] - 17*Q5[1,t] - 17*Q5[1,t+1] +\
                    4*Q5[2,t-1] + 4*Q5[2,t] + 4*Q5[2,t+1] + Q5[3,t-1] + Q5[3,t] + Q5[3,t+1])))/(4.*dz)

            elif z == 2:

                GradE[5*t*Lz + 5*z] =  scale*-(dt*(-2*(ks + kt)*Q1[0,t-1] - 2*(ks + kt)*Q1[0,t] - 2*ks*Q1[0,t+1] - 2*kt*Q1[0,t+1] + 4*ks*Q1[1,t-1] +\
                    4*kt*Q1[1,t-1] + 4*ks*Q1[1,t] + 4*kt*Q1[1,t] + 4*ks*Q1[1,t+1] + 4*kt*Q1[1,t+1] - 3*ks*Q1[2,t-1] -\
                    3*kt*Q1[2,t-1] - 3*ks*Q1[2,t] - 3*kt*Q1[2,t] - 3*ks*Q1[2,t+1] - 3*kt*Q1[2,t+1] + ks*Q1[4,t-1] +\
                    kt*Q1[4,t-1] + ks*Q1[4,t] + kt*Q1[4,t] + ks*Q1[4,t+1] + kt*Q1[4,t+1] + 4*dz*kt*q0*Q2[0,t-1] +\
                    4*dz*kt*q0*Q2[0,t] + 4*dz*kt*q0*Q2[0,t+1] - 8*dz*kt*q0*Q2[1,t-1] - 8*dz*kt*q0*Q2[1,t] -\
                    8*dz*kt*q0*Q2[1,t+1] + 8*dz*kt*q0*Q2[3,t-1] + 8*dz*kt*q0*Q2[3,t] + 8*dz*kt*q0*Q2[3,t+1] -\
                    2*ks*Q4[0,t-1] - 2*ks*Q4[0,t] - 2*ks*Q4[0,t+1] + 4*ks*Q4[1,t-1] + 4*ks*Q4[1,t] + 4*ks*Q4[1,t+1] -\
                    3*ks*Q4[2,t-1] - 3*ks*Q4[2,t] - 3*ks*Q4[2,t+1] + ks*Q4[4,t-1] + ks*Q4[4,t] + ks*Q4[4,t+1]))/(4.*dz)

                GradE[5*t*Lz + 5*z + 1] =  scale*(dt*kt*(2*Q2[0,t-1] + 2*Q2[0,t] + 2*Q2[0,t+1] - 4*Q2[1,t-1] - 4*Q2[1,t] - 4*Q2[1,t+1] +\
                    3*Q2[2,t-1] + 3*Q2[2,t] + 3*Q2[2,t+1] - Q2[4,t-1] - Q2[4,t] - Q2[4,t+1] +\
                    2*dz*q0*(Q1[0,t-1] + Q1[0,t] + Q1[0,t+1] - 2*Q1[1,t-1] - 2*Q1[1,t] - 2*Q1[1,t+1] +\
                    2*Q1[3,t-1] + 2*Q1[3,t] + 2*Q1[3,t+1] - Q4[0,t-1] - Q4[0,t] - Q4[0,t+1] +\
                    2*(Q4[1,t-1] + Q4[1,t] + Q4[1,t+1] - Q4[3,t-1] - Q4[3,t] - Q4[3,t+1]))))/(2.*dz)

                GradE[5*t*Lz + 5*z + 2] =  scale*-(dt*(-2*(ks + kt)*Q3[0,t-1] - 2*(ks + kt)*Q3[0,t] - 2*ks*Q3[0,t+1] - 2*kt*Q3[0,t+1] + 4*ks*Q3[1,t-1] +\
                    4*kt*Q3[1,t-1] + 4*ks*Q3[1,t] + 4*kt*Q3[1,t] + 4*ks*Q3[1,t+1] + 4*kt*Q3[1,t+1] - 3*ks*Q3[2,t-1] -\
                    3*kt*Q3[2,t-1] - 3*ks*Q3[2,t] - 3*kt*Q3[2,t] - 3*ks*Q3[2,t+1] - 3*kt*Q3[2,t+1] + ks*Q3[4,t-1] +\
                    kt*Q3[4,t-1] + ks*Q3[4,t] + kt*Q3[4,t] + ks*Q3[4,t+1] + kt*Q3[4,t+1] + 4*dz*kt*q0*Q5[0,t-1] +\
                    4*dz*kt*q0*Q5[0,t] + 4*dz*kt*q0*Q5[0,t+1] - 8*dz*kt*q0*Q5[1,t-1] - 8*dz*kt*q0*Q5[1,t] -\
                    8*dz*kt*q0*Q5[1,t+1] + 8*dz*kt*q0*Q5[3,t-1] + 8*dz*kt*q0*Q5[3,t] + 8*dz*kt*q0*Q5[3,t+1]))/(4.*dz)

                GradE[5*t*Lz + 5*z + 3] =  scale*(dt*(2*ks*Q1[0,t-1] + 2*ks*Q1[0,t] + 2*ks*Q1[0,t+1] - 4*ks*Q1[1,t-1] - 4*ks*Q1[1,t] - 4*ks*Q1[1,t+1] +\
                    3*ks*Q1[2,t-1] + 3*ks*Q1[2,t] + 3*ks*Q1[2,t+1] - ks*Q1[4,t-1] - ks*Q1[4,t] - ks*Q1[4,t+1] +\
                    4*dz*kt*q0*(Q2[0,t-1] + Q2[0,t] + Q2[0,t+1] -\
                    2*(Q2[1,t-1] + Q2[1,t] + Q2[1,t+1] - Q2[3,t-1] - Q2[3,t] - Q2[3,t+1])) + 2*ks*Q4[0,t-1] +\
                    2*kt*Q4[0,t-1] + 2*ks*Q4[0,t] + 2*kt*Q4[0,t] + 2*ks*Q4[0,t+1] + 2*kt*Q4[0,t+1] - 4*ks*Q4[1,t-1] -\
                    4*kt*Q4[1,t-1] - 4*ks*Q4[1,t] - 4*kt*Q4[1,t] - 4*ks*Q4[1,t+1] - 4*kt*Q4[1,t+1] + 3*ks*Q4[2,t-1] +\
                    3*kt*Q4[2,t-1] + 3*ks*Q4[2,t] + 3*kt*Q4[2,t] + 3*ks*Q4[2,t+1] + 3*kt*Q4[2,t+1] - ks*Q4[4,t-1] -\
                    kt*Q4[4,t-1] - ks*Q4[4,t] - kt*Q4[4,t] - (ks + kt)*Q4[4,t+1]))/(4.*dz)

                GradE[5*t*Lz + 5*z + 4] =  scale*(dt*(ks*(2*Q5[0,t-1] + 2*Q5[0,t] + 2*Q5[0,t+1] - 4*Q5[1,t-1] - 4*Q5[1,t] - 4*Q5[1,t+1] +\
                    3*Q5[2,t-1] + 3*Q5[2,t] + 3*Q5[2,t+1] - Q5[4,t-1] - Q5[4,t] - Q5[4,t+1]) +\
                    kt*(4*dz*q0*(Q3[0,t-1] + Q3[0,t] + Q3[0,t+1] -\
                    2*(Q3[1,t-1] + Q3[1,t] + Q3[1,t+1] - Q3[3,t-1] - Q3[3,t] - Q3[3,t+1])) + 2*Q5[0,t-1] +\
                    2*Q5[0,t] + 2*Q5[0,t+1] - 4*Q5[1,t-1] - 4*Q5[1,t] - 4*Q5[1,t+1] + 3*Q5[2,t-1] + 3*Q5[2,t] +\
                    3*Q5[2,t+1] - Q5[4,t-1] - Q5[4,t] - Q5[4,t+1])))/(4.*dz)

            elif z == Lz-1:

                GradE[5*t*Lz + 5*z] =  scale*(dt*((ks + kt)*Q1[Lz-3,t-1] + (ks + kt)*Q1[Lz-3,t] + ks*Q1[Lz-3,t+1] + kt*Q1[Lz-3,t+1] -\
                    6*ks*Q1[Lz-2,t-1] - 6*kt*Q1[Lz-2,t-1] - 6*ks*Q1[Lz-2,t] - 6*kt*Q1[Lz-2,t] -\
                    6*ks*Q1[Lz-2,t+1] - 6*kt*Q1[Lz-2,t+1] + 5*ks*Q1[Lz-1,t-1] + 5*kt*Q1[Lz-1,t-1] +\
                    5*ks*Q1[Lz-1,t] + 5*kt*Q1[Lz-1,t] + 5*ks*Q1[Lz-1,t+1] + 5*kt*Q1[Lz-1,t+1] -\
                    2*dz*kt*q0*Q2[Lz-3,t-1] - 2*dz*kt*q0*Q2[Lz-3,t] - 2*dz*kt*q0*Q2[Lz-3,t+1] +\
                    10*dz*kt*q0*Q2[Lz-2,t-1] + 10*dz*kt*q0*Q2[Lz-2,t] + 10*dz*kt*q0*Q2[Lz-2,t+1] +\
                    ks*Q4[Lz-3,t-1] + ks*Q4[Lz-3,t] + ks*Q4[Lz-3,t+1] - 6*ks*Q4[Lz-2,t-1] -\
                    6*ks*Q4[Lz-2,t] - 6*ks*Q4[Lz-2,t+1] + 5*ks*Q4[Lz-1,t-1] + 5*ks*Q4[Lz-1,t] +\
                    5*ks*Q4[Lz-1,t+1]))/(2.*dz)


                GradE[5*t*Lz + 5*z + 1] =  scale*(dt*kt*(Q2[Lz-3,t-1] + Q2[Lz-3,t] + Q2[Lz-3,t+1] - 6*Q2[Lz-2,t-1] - 6*Q2[Lz-2,t] -\
                    6*Q2[Lz-2,t+1] + 5*(Q2[Lz-1,t-1] + Q2[Lz-1,t] + Q2[Lz-1,t+1]) +\
                    dz*q0*(Q1[Lz-3,t-1] + Q1[Lz-3,t] + Q1[Lz-3,t+1] - 5*Q1[Lz-2,t-1] - 5*Q1[Lz-2,t] -\
                    5*Q1[Lz-2,t+1] - Q4[Lz-3,t-1] - Q4[Lz-3,t] - Q4[Lz-3,t+1] +\
                    5*(Q4[Lz-2,t-1] + Q4[Lz-2,t] + Q4[Lz-2,t+1]))))/dz

                GradE[5*t*Lz + 5*z + 2] =  scale*(dt*((ks + kt)*Q3[Lz-3,t-1] + (ks + kt)*Q3[Lz-3,t] + ks*Q3[Lz-3,t+1] + kt*Q3[Lz-3,t+1] -\
                    6*ks*Q3[Lz-2,t-1] - 6*kt*Q3[Lz-2,t-1] - 6*ks*Q3[Lz-2,t] - 6*kt*Q3[Lz-2,t] -\
                    6*ks*Q3[Lz-2,t+1] - 6*kt*Q3[Lz-2,t+1] + 5*ks*Q3[Lz-1,t-1] + 5*kt*Q3[Lz-1,t-1] +\
                    5*ks*Q3[Lz-1,t] + 5*kt*Q3[Lz-1,t] + 5*ks*Q3[Lz-1,t+1] + 5*kt*Q3[Lz-1,t+1] -\
                    2*dz*kt*q0*Q5[Lz-3,t-1] - 2*dz*kt*q0*Q5[Lz-3,t] - 2*dz*kt*q0*Q5[Lz-3,t+1] +\
                    10*dz*kt*q0*Q5[Lz-2,t-1] + 10*dz*kt*q0*Q5[Lz-2,t] + 10*dz*kt*q0*Q5[Lz-2,t+1]))/(2.*dz)

                GradE[5*t*Lz + 5*z + 3] =  scale*(dt*(ks*Q1[Lz-3,t-1] + ks*Q1[Lz-3,t] + ks*Q1[Lz-3,t+1] - 6*ks*Q1[Lz-2,t-1] -\
                    6*ks*Q1[Lz-2,t] - 6*ks*Q1[Lz-2,t+1] + 5*ks*Q1[Lz-1,t-1] + 5*ks*Q1[Lz-1,t] +\
                    5*ks*Q1[Lz-1,t+1] + 2*dz*kt*q0*(Q2[Lz-3,t-1] + Q2[Lz-3,t] + Q2[Lz-3,t+1] -\
                    5*(Q2[Lz-2,t-1] + Q2[Lz-2,t] + Q2[Lz-2,t+1])) + ks*Q4[Lz-3,t-1] +\
                    kt*Q4[Lz-3,t-1] + ks*Q4[Lz-3,t] + kt*Q4[Lz-3,t] + ks*Q4[Lz-3,t+1] + kt*Q4[Lz-3,t+1] -\
                    6*ks*Q4[Lz-2,t-1] - 6*kt*Q4[Lz-2,t-1] - 6*ks*Q4[Lz-2,t] - 6*kt*Q4[Lz-2,t] -\
                    6*ks*Q4[Lz-2,t+1] - 6*kt*Q4[Lz-2,t+1] + 5*ks*Q4[Lz-1,t-1] + 5*kt*Q4[Lz-1,t-1] +\
                    5*ks*Q4[Lz-1,t] + 5*kt*Q4[Lz-1,t] + 5*(ks + kt)*Q4[Lz-1,t+1]))/(2.*dz)

                GradE[5*t*Lz + 5*z + 4] =  scale*(dt*(ks*(Q5[Lz-3,t-1] + Q5[Lz-3,t] + Q5[Lz-3,t+1] - 6*Q5[Lz-2,t-1] - 6*Q5[Lz-2,t] -\
                    6*Q5[Lz-2,t+1] + 5*(Q5[Lz-1,t-1] + Q5[Lz-1,t] + Q5[Lz-1,t+1])) +\
                    kt*(2*dz*q0*(Q3[Lz-3,t-1] + Q3[Lz-3,t] + Q3[Lz-3,t+1] -\
                    5*(Q3[Lz-2,t-1] + Q3[Lz-2,t] + Q3[Lz-2,t+1])) + Q5[Lz-3,t-1] + Q5[Lz-3,t] +\
                    Q5[Lz-3,t+1] - 6*Q5[Lz-2,t-1] - 6*Q5[Lz-2,t] - 6*Q5[Lz-2,t+1] +\
                    5*(Q5[Lz-1,t-1] + Q5[Lz-1,t] + Q5[Lz-1,t+1]))))/(2.*dz)

                GradE[5*t*Lz + 5*z] +=  scale*(ws/2.0)*(2*(Q1[z,t]-Qt1)-2*(Qt1+Qt4-Q1[z,t]-Q4[z,t]))

                GradE[5*t*Lz + 5*z + 1] +=  scale*2.0*ws*(Q2[z,t]-Qt2)

                GradE[5*t*Lz + 5*z + 2] +=  scale*2.0*ws*(Q3[z,t]-Qt3)

                GradE[5*t*Lz + 5*z + 3] += scale*(ws/2.0)*(2*(Q4[z,t]-Qt4)-2*(Qt1+Qt4-Q1[z,t]-Q4[z,t]))

                GradE[5*t*Lz + 5*z + 4] += scale*2.0*ws*(Q5[z,t]-Qt5)

            elif z == Lz-2:

                GradE[5*t*Lz + 5*z] =  scale*(dt*(-((ks + kt)*Q1[Lz-4,t-1]) - (ks + kt)*Q1[Lz-4,t] - ks*Q1[Lz-4,t+1] - kt*Q1[Lz-4,t+1] +\
                    ks*Q1[Lz-2,t-1] + kt*Q1[Lz-2,t-1] + ks*Q1[Lz-2,t] + kt*Q1[Lz-2,t] + ks*Q1[Lz-2,t+1] +\
                    kt*Q1[Lz-2,t+1] + 8*dz*kt*q0*Q2[Lz-3,t-1] + 8*dz*kt*q0*Q2[Lz-3,t] + 8*dz*kt*q0*Q2[Lz-3,t+1] -\
                    4*dz*kt*q0*Q2[Lz-1,t-1] - 4*dz*kt*q0*Q2[Lz-1,t] - 4*dz*kt*q0*Q2[Lz-1,t+1] - ks*Q4[Lz-4,t-1] -\
                    ks*Q4[Lz-4,t] - ks*Q4[Lz-4,t+1] + ks*Q4[Lz-2,t-1] + ks*Q4[Lz-2,t] + ks*Q4[Lz-2,t+1]))/(4.*dz)

                GradE[5*t*Lz + 5*z + 1] =  scale*(dt*kt*(-Q2[Lz-4,t-1] - Q2[Lz-4,t] - Q2[Lz-4,t+1] + Q2[Lz-2,t-1] + Q2[Lz-2,t] +\
                    Q2[Lz-2,t+1] + 2*dz*q0*(-2*Q1[Lz-3,t-1] - 2*Q1[Lz-3,t] - 2*Q1[Lz-3,t+1] +\
                    Q1[Lz-1,t-1] + Q1[Lz-1,t] + Q1[Lz-1,t+1] + 2*Q4[Lz-3,t-1] + 2*Q4[Lz-3,t] +\
                    2*Q4[Lz-3,t+1] - Q4[Lz-1,t-1] - Q4[Lz-1,t] - Q4[Lz-1,t+1])))/(2.*dz)

                GradE[5*t*Lz + 5*z + 2] =  scale*(dt*(-((ks + kt)*Q3[Lz-4,t-1]) - (ks + kt)*Q3[Lz-4,t] - ks*Q3[Lz-4,t+1] - kt*Q3[Lz-4,t+1] +\
                    ks*Q3[Lz-2,t-1] + kt*Q3[Lz-2,t-1] + ks*Q3[Lz-2,t] + kt*Q3[Lz-2,t] + ks*Q3[Lz-2,t+1] +\
                    kt*Q3[Lz-2,t+1] + 8*dz*kt*q0*Q5[Lz-3,t-1] + 8*dz*kt*q0*Q5[Lz-3,t] + 8*dz*kt*q0*Q5[Lz-3,t+1] -\
                    4*dz*kt*q0*Q5[Lz-1,t-1] - 4*dz*kt*q0*Q5[Lz-1,t] - 4*dz*kt*q0*Q5[Lz-1,t+1]))/(4.*dz)

                GradE[5*t*Lz + 5*z + 3] =  scale*(dt*(ks*(-Q1[Lz-4,t-1] - Q1[Lz-4,t] - Q1[Lz-4,t+1] + Q1[Lz-2,t-1] + Q1[Lz-2,t] +\
                    Q1[Lz-2,t+1]) + 4*dz*kt*q0*(-2*Q2[Lz-3,t-1] - 2*Q2[Lz-3,t] - 2*Q2[Lz-3,t+1] +\
                    Q2[Lz-1,t-1] + Q2[Lz-1,t] + Q2[Lz-1,t+1]) - (ks + kt)*Q4[Lz-4,t-1] -\
                    (ks + kt)*Q4[Lz-4,t] - (ks + kt)*Q4[Lz-4,t+1] + (ks + kt)*Q4[Lz-2,t-1] + (ks + kt)*Q4[Lz-2,t] +\
                    (ks + kt)*Q4[Lz-2,t+1]))/(4.*dz)

                GradE[5*t*Lz + 5*z + 4] =  scale*(dt*(4*dz*kt*q0*(-2*Q3[Lz-3,t-1] - 2*Q3[Lz-3,t] - 2*Q3[Lz-3,t+1] + Q3[Lz-1,t-1] +\
                    Q3[Lz-1,t] + Q3[Lz-1,t+1]) - (ks + kt)*Q5[Lz-4,t-1] - (ks + kt)*Q5[Lz-4,t] -\
                    (ks + kt)*Q5[Lz-4,t+1] + (ks + kt)*Q5[Lz-2,t-1] + (ks + kt)*Q5[Lz-2,t] + (ks + kt)*Q5[Lz-2,t+1]))/(4.*dz)

            elif z == Lz-3:

                GradE[5*t*Lz + 5*z] =  scale*-(dt*((ks + kt)*Q1[Lz-5,t-1] + (ks + kt)*Q1[Lz-5,t] + ks*Q1[Lz-5,t+1] + kt*Q1[Lz-5,t+1] -\
                    2*ks*Q1[Lz-3,t-1] - 2*kt*Q1[Lz-3,t-1] - 2*ks*Q1[Lz-3,t] - 2*kt*Q1[Lz-3,t] -\
                    2*ks*Q1[Lz-3,t+1] - 2*kt*Q1[Lz-3,t+1] + ks*Q1[Lz-1,t-1] + kt*Q1[Lz-1,t-1] +\
                    ks*Q1[Lz-1,t] + kt*Q1[Lz-1,t] + ks*Q1[Lz-1,t+1] + kt*Q1[Lz-1,t+1] -\
                    8*dz*kt*q0*Q2[Lz-4,t-1] - 8*dz*kt*q0*Q2[Lz-4,t] - 8*dz*kt*q0*Q2[Lz-4,t+1] +\
                    8*dz*kt*q0*Q2[Lz-2,t-1] + 8*dz*kt*q0*Q2[Lz-2,t] + 8*dz*kt*q0*Q2[Lz-2,t+1] +\
                    ks*Q4[Lz-5,t-1] + ks*Q4[Lz-5,t] + ks*Q4[Lz-5,t+1] - 2*ks*Q4[Lz-3,t-1] -\
                    2*ks*Q4[Lz-3,t] - 2*ks*Q4[Lz-3,t+1] + ks*Q4[Lz-1,t-1] + ks*Q4[Lz-1,t] + ks*Q4[Lz-1,t+1]))/(4.*dz)

                GradE[5*t*Lz + 5*z + 1] =  scale*-(dt*kt*(Q2[Lz-5,t-1] + Q2[Lz-5,t] + Q2[Lz-5,t+1] - 2*Q2[Lz-3,t-1] - 2*Q2[Lz-3,t] -\
                    2*Q2[Lz-3,t+1] + Q2[Lz-1,t-1] + Q2[Lz-1,t] + Q2[Lz-1,t+1] +\
                    4*dz*q0*(Q1[Lz-4,t-1] + Q1[Lz-4,t] + Q1[Lz-4,t+1] - Q1[Lz-2,t-1] - Q1[Lz-2,t] -\
                    Q1[Lz-2,t+1] - Q4[Lz-4,t-1] - Q4[Lz-4,t] - Q4[Lz-4,t+1] + Q4[Lz-2,t-1] +\
                    Q4[Lz-2,t] + Q4[Lz-2,t+1])))/(2.*dz)

                GradE[5*t*Lz + 5*z + 2] =  scale*-(dt*((ks + kt)*Q3[Lz-5,t-1] + (ks + kt)*Q3[Lz-5,t] + ks*Q3[Lz-5,t+1] + kt*Q3[Lz-5,t+1] -\
                    2*ks*Q3[Lz-3,t-1] - 2*kt*Q3[Lz-3,t-1] - 2*ks*Q3[Lz-3,t] - 2*kt*Q3[Lz-3,t] -\
                    2*ks*Q3[Lz-3,t+1] - 2*kt*Q3[Lz-3,t+1] + ks*Q3[Lz-1,t-1] + kt*Q3[Lz-1,t-1] +\
                    ks*Q3[Lz-1,t] + kt*Q3[Lz-1,t] + ks*Q3[Lz-1,t+1] + kt*Q3[Lz-1,t+1] -\
                    8*dz*kt*q0*Q5[Lz-4,t-1] - 8*dz*kt*q0*Q5[Lz-4,t] - 8*dz*kt*q0*Q5[Lz-4,t+1] +\
                    8*dz*kt*q0*Q5[Lz-2,t-1] + 8*dz*kt*q0*Q5[Lz-2,t] + 8*dz*kt*q0*Q5[Lz-2,t+1]))/(4.*dz)

                GradE[5*t*Lz + 5*z + 3] =  scale*-(dt*(ks*Q1[Lz-5,t-1] + ks*Q1[Lz-5,t] + ks*Q1[Lz-5,t+1] - 2*ks*Q1[Lz-3,t-1] - \
                    2*ks*Q1[Lz-3,t] - 2*ks*Q1[Lz-3,t+1] + ks*Q1[Lz-1,t-1] + ks*Q1[Lz-1,t] +\
                    ks*Q1[Lz-1,t+1] + 8*dz*kt*q0*(Q2[Lz-4,t-1] + Q2[Lz-4,t] + Q2[Lz-4,t+1] -\
                    Q2[Lz-2,t-1] - Q2[Lz-2,t] - Q2[Lz-2,t+1]) + ks*Q4[Lz-5,t-1] + kt*Q4[Lz-5,t-1] +\
                    ks*Q4[Lz-5,t] + kt*Q4[Lz-5,t] + ks*Q4[Lz-5,t+1] + kt*Q4[Lz-5,t+1] - 2*ks*Q4[Lz-3,t-1] -\
                    2*kt*Q4[Lz-3,t-1] - 2*ks*Q4[Lz-3,t] - 2*kt*Q4[Lz-3,t] - 2*ks*Q4[Lz-3,t+1] -\
                    2*kt*Q4[Lz-3,t+1] + ks*Q4[Lz-1,t-1] + kt*Q4[Lz-1,t-1] + ks*Q4[Lz-1,t] + kt*Q4[Lz-1,t] +\
                    (ks + kt)*Q4[Lz-1,t+1]))/(4.*dz)

                GradE[5*t*Lz + 5*z + 4] =  scale*(dt*(8*dz*kt*q0*(-Q3[Lz-4,t-1] - Q3[Lz-4,t] - Q3[Lz-4,t+1] + Q3[Lz-2,t-1] + Q3[Lz-2,t] + \
                    Q3[Lz-2,t+1]) - (ks + kt)*Q5[Lz-5,t-1] - (ks + kt)*Q5[Lz-5,t] - (ks + kt)*Q5[Lz-5,t+1] + \
                    2*(ks + kt)*Q5[Lz-3,t-1] + 2*(ks + kt)*Q5[Lz-3,t] + 2*(ks + kt)*Q5[Lz-3,t+1] - \
                    (ks + kt)*Q5[Lz-1,t-1] - (ks + kt)*Q5[Lz-1,t] - (ks + kt)*Q5[Lz-1,t+1]))/(4.*dz)

            else:

                GradE[5*t*Lz + 5*z] += scale*dt*kt*q0*(Q2[z-2,t-1] + Q2[z-2,t]\
                      + Q2[z-2,t+1] + Q2[z-1,t-1] + Q2[z-1,t] + Q2[z-1,t+1]\
                      - Q2[z+1,t-1] - Q2[z+1,t] - Q2[z+1,t+1] - Q2[z+2,t-1] - Q2[z+2,t] - Q2[z+2,t+1])

                GradE[5*t*Lz + 5*z + 1] += scale*dt*kt*q0*(-Q1[z-2,t-1] - Q1[z-2,t] - Q1[z-2,t+1] - Q1[z-1,t-1] - Q1[z-1,t] - \
                        Q1[z-1,t+1] + Q1[z+1,t-1] + Q1[z+1,t] + Q1[z+1,t+1] + Q1[z+2,t-1] + Q1[z+2,t] + \
                        Q1[z+2,t+1] + Q4[z-2,t-1] + Q4[z-2,t] + Q4[z-2,t+1] + Q4[z-1,t-1] + Q4[z-1,t] +\
                        Q4[z-1,t+1] - Q4[z+1,t-1] - Q4[z+1,t] - Q4[z+1,t+1] - Q4[z+2,t-1] - Q4[z+2,t] - \
                        Q4[z+2,t+1])

                GradE[5*t*Lz + 5*z + 2] += scale*dt*kt*q0*(Q5[z-2,t-1] + Q5[z-2,t] + Q5[z-2,t+1] + Q5[z-1,t-1] + Q5[z-1,t] + \
                        Q5[z-1,t+1] - Q5[z+1,t-1] - Q5[z+1,t] - Q5[z+1,t+1] - Q5[z+2,t-1] - Q5[z+2,t] - \
                        Q5[z+2,t+1])

                GradE[5*t*Lz + 5*z + 3] += scale*dt*kt*q0*(-Q2[z-2,t-1] - Q2[z-2,t] - Q2[z-2,t+1] - Q2[z-1,t-1] - Q2[z-1,t] - \
                        Q2[z-1,t+1] + Q2[z+1,t-1] + Q2[z+1,t] + Q2[z+1,t+1] + Q2[z+2,t-1] + Q2[z+2,t] + \
                        Q2[z+2,t+1])

                GradE[5*t*Lz + 5*z + 4] += scale*dt*kt*q0*(-Q3[z-2,t-1] - Q3[z-2,t] - Q3[z-2,t+1] - Q3[z-1,t-1] - Q3[z-1,t] - \
                        Q3[z-1,t+1] + Q3[z+1,t-1] + Q3[z+1,t] + Q3[z+1,t+1] + Q3[z+2,t-1] + Q3[z+2,t] + \
                        Q3[z+2,t+1])

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

    #GradE = GradE*1.0/(1.0*Lz*Lt)
    return(GradE)
