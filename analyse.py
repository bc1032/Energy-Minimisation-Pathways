import numpy as np
import pandas as pd
import pylab as pl
import sys
import os
import matplotlib.pyplot as plt
import math as m
from mpl_toolkits.mplot3d import Axes3D


lz,lt = 101,400
frames = 1
data = np.loadtxt("energyarray.dat")
#data = np.loadtxt("initialguess.dat")
if not os.path.exists('data'):
    os.makedirs('data')

#for i in range(0, lt):
#    if not os.path.exists('../spherical/%d' % i):
#        os.makedirs('../spherical/%d' % i)

r,theta,phi=[],[],[]
def cart2sph(x,y,z):
    XsqPlusYsq = x**2 + y**2
    r2 = m.sqrt(XsqPlusYsq + z**2)               # r
    elev = m.atan2(z,m.sqrt(XsqPlusYsq))     # theta
    az = m.atan2(y,x)                          # phi
    r.append(r2)
    theta.append(elev)
    phi.append(az)
    phifile.write("%f\n" % (az))
    rfile.write("%f\n" % r2 )
    thetafile.write("%f\n" % elev )
    sphfile.write("%f   %f  %f\n" % (r2,elev,az))
    #print(r2,elev,az)
    return (r2, elev, az)

cent = int(lt*lz/2.0)


#data = np.loadtxt("string")
print("File Loaded")
mod = []
nx,ny,nz=[],[],[]
phifile = open("data/phi.dat", 'w')
rfile = open("data/r.dat", 'w')
thetafile = open("data/theta.dat", 'w')
sphfile = open("data/sphericalcoords.dat", 'w')
cartfile = open("data/cartesiancoords.dat", 'w')
evalfile = open("data/eval.dat", 'w')
threedfile1 = open("data/3dfile1.dat", 'w')
threedfile2 = open("data/3dfile2.dat", 'w')
threedfile3 = open("data/3dfile3.dat", 'w')
sumqfile = open("data/sumQ.dat",'w')
eval1file = open("data/eval1.dat",'w')
eval2file = open("data/eval2.dat",'w')
eval3file = open("data/eval3.dat",'w')
nxfile = open("data/nx.dat",'w')
nyfile = open("data/ny.dat",'w')
nzfile = open("data/nz.dat",'w')
Q1file = open("data/Q1.dat",'w')
Q2file = open("data/Q2.dat",'w')
Q3file = open("data/Q3.dat",'w')
Q4file = open("data/Q4.dat",'w')
Q5file = open("data/Q5.dat",'w')
modqfile = open("data/modQ.dat", 'w')
Q1m,Q2m,Q3m,Q4m,Q5m = [],[],[],[],[]
eval1,eval2,eval3 = [],[],[]
for n in range(0,frames):
    print(n)
    for i in range(0,lt):
        for k in range(0,(lz)):

            Q1 = data[5*i*lz + 5*k + 5*n*lt*lz]
            Q2 = data[5*i*lz + 5*k + 5*n*lt*lz + 1]
            Q3 = data[5*i*lz + 5*k + 5*n*lt*lz + 2]
            Q4 = data[5*i*lz + 5*k + 5*n*lt*lz + 3]
            Q5 = data[5*i*lz + 5*k + 5*n*lt*lz + 4]
            M = np.array([[Q1,Q2,Q3],[Q2,Q4,Q5],[Q3,Q5,-Q1-Q4]])

            mod.append((Q1**2 + Q2**2 + Q3**2 + Q4**2 + Q5**2 + Q1*Q4))
            modqfile.write("%f\n" % (Q1**2 + Q2**2 + Q3**2 + Q4**2 + Q5**2 + Q1*Q4))
            sumqfile.write("%f\n" % (Q1+Q2+Q3+Q4+Q5))

            eval, evec = np.linalg.eig(M)
            idx = eval.argsort()[::-1]
            eval = eval[idx]
            evec = evec[:,idx]
            x = abs(evec[0,0])
            y = abs(evec[1,0])
            z = abs(evec[2,0])

            nx.append(x)
            ny.append(y)
            nz.append(z)
            eval1file.write("%f\n" % eval[0])
            eval2file.write("%f\n" % eval[1])
            eval3file.write("%f\n" % eval[2])
            Q1m.append(Q1)
            Q2m.append(Q2)
            Q3m.append(Q3)
            Q4m.append(Q4)
            Q5m.append(Q5)
            Q1file.write("%f\n" % Q1)
            Q2file.write("%f\n" % Q2)
            Q3file.write("%f\n" % Q3)
            Q4file.write("%f\n" % Q4)
            Q5file.write("%f\n" % Q5)
            eval1.append(eval[0])
            eval2.append(eval[1])
            eval3.append(eval[2])
            evalfile.write("%f  %f  %f\n" % (eval[0],eval[1],eval[2]))
            nxfile.write("%f\n" % x)
            nyfile.write("%f\n" % y)
            nzfile.write("%f\n" % z)
            cartfile.write("%f  %f  %f\n" % (x,y,z))
            cart2sph(x,y,z)
            threedfile1.write("%f\n" % (Q1-Q4) )
            threedfile2.write("%f\n" % (Q2 ))
            threedfile3.write("%f\n" % (Q1+Q4) )
thetafile.close()
phifile.close()
rfile.close()
sphfile.close()
cartfile.close()
Q1file.close()
Q2file.close()
Q3file.close()
Q4file.close()
Q5file.close()
modqfile.close()
sumqfile.close()
evalfile.close()
eval1file.close()
eval2file.close()
eval3file.close()
threedfile1.close()
threedfile2.close()
threedfile3.close()
nxfile.close()
nyfile.close()
nzfile.close()
