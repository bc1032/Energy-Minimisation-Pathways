#Libraries
import math
import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal
import scipy
from scipy.optimize import minimize

def initialise(splitt,alpha, beta, gamma,Q3,Q5, Length, Time, winding0, winding1, a, b, c):
    timesplit = int(Time / splitt)
    timesplitf = int((splitt-1)*(Time / splitt))
    lowz = int(Length/5)
    highz = 4*int(Length/5)

    strength = (b + math.sqrt(b**2 + 24*a*c))/(4.0*c)
    for i in range(0,Length):
        theta = math.pi/2.0
        j = 0
        n = 0.0
        N = timesplitf-timesplit
        const = 0.0
        for j in range(0,Time):
            alpha[i,j] = -strength*(1.0/12.0)*(1.0 + 3.0*math.cos(2*theta))
            beta[i,j] = 0.5*strength*(math.sin(theta)**2)*(\
                ((Time-1-j)/(Time-1))*math.cos(winding0*i*2*math.pi/(Length-1))\
            +   ((j)/(Time-1))*       math.cos(winding1*i*2*math.pi/(Length-1))\
            )
            gamma[i,j] = 0.5*strength*(math.sin(theta)**2)*(\
                ((Time-1-j)/(Time-1))*math.sin(winding0*i*2*math.pi/(Length-1))\
            +   ((j)/(Time-1))*       math.sin(winding1*i*2*math.pi/(Length-1))\
            )
            Q3[i,j] = strength*(math.sin(theta))*math.cos(theta)*(\
                ((Time-1-j)/(Time-1))*math.cos(winding0*i*math.pi/(Length-1))\
            +   ((j)/(Time-1))*       math.cos(winding1*i*math.pi/(Length-1)))
            Q5[i,j] = strength*(math.sin(theta))*math.cos(theta)*(\
                ((Time-1-j)/(Time-1))*math.sin(winding0*i*math.pi/(Length-1))\
            +   ((j)/(Time-1))*       math.sin(winding1*i*math.pi/(Length-1))\
            )


        Nt = timesplitf-timesplit
        nt = 0
        for j in range(timesplit,timesplitf):
        #if i >= lowz and i <= highz:
            N = highz-lowz
            n = lowz
            theta = math.pi/2.0 + math.sin(nt*math.pi/((Nt-1)))*math.pi/2.0

            for i in range(lowz,highz):
                alpha[i,j] = -strength*(1.0/12.0)*(1.0 + 3.0*math.cos(2*theta))
                beta[i,j] = 0.5*strength*(math.sin(theta)**2)*(\
                    ((Time-1-j)/(Time-1))*math.cos(winding0*i*2*math.pi/(Length-1))\
                +   ((j)/(Time-1))*       math.cos(winding1*i*2*math.pi/(Length-1))\
                )
                gamma[i,j] = 0.5*strength*(math.sin(theta)**2)*(\
                    ((Time-1-j)/(Time-1))*math.sin(winding0*i*2*math.pi/(Length-1))\
                +   ((j)/(Time-1))*       math.sin(winding1*i*2*math.pi/(Length-1))\
                )
                Q3[i,j] = strength*(math.sin(theta))*math.cos(theta)*(\
                    ((Time-1-j)/(Time-1))*math.cos(winding0*i*math.pi/(Length-1))\
                +   ((j)/(Time-1))*       math.cos(winding1*i*math.pi/(Length-1)))
                Q5[i,j] = strength*(math.sin(theta))*math.cos(theta)*(\
                    ((Time-1-j)/(Time-1))*math.sin(winding0*i*math.pi/(Length-1))\
                +   ((j)/(Time-1))*       math.sin(winding1*i*math.pi/(Length-1))\
                )

                n += 1.0
            nt += 1.0

    return(alpha, beta, gamma,Q3,Q5, Length, Time, winding0, winding1, a, b, c)
