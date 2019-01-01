# This is module for for solving a coupled first odrder differential equation
""" The general form of the coupled differential equation is
            dy/dx = f1(x,y,z)
            dz/dx = f2(x,y,z)
Since the ODE's are coupled we need to solve the K's parallely """

from math import *
from matplotlib.pylab import *

def ODE_runge_4th_CD(f1,f2,y0,z0,x_min,x_max):
    x = linspace(x_min,x_max,(x_max-x_min)*100+1)
    h = x[1]-x[0]
    y = zeros(len(x))
    z = zeros(len(x))
    y[0] = y0
    z[0] = z0
    for i in range(len(x)-1):
        ky1 = f1(x[i],y[i],z[i])*h
        kz1 = f2(x[i],y[i],z[i])*h
        ky2 = f1(x[i] + 0.5*h,y[i] + 0.5*ky1,z[i] + 0.5*kz1)*h
        kz2 = f2(x[i] + 0.5*h,y[i] + 0.5*ky1,z[i] + 0.5*kz1)*h
        ky3 = f1(x[i] + 0.5*h,y[i] + 0.5*ky2,z[i] + 0.5*kz2)*h
        kz3 = f2(x[i] + 0.5*h,y[i] + 0.5*ky2,z[i] + 0.5*kz2)*h
        ky4 = f1(x[i] + h,y[i] + ky3,z[i] + kz3)*h
        kz4 = f2(x[i] + h,y[i] + ky3,z[i] + kz3)*h
        y[i+1] = y[i] + (ky1+ky2+ky2+ky3+ky3+ky4)/6
        z[i+1] = z[i] + (kz1+kz2+kz2+kz3+kz3+kz4)/6
    return z,y,x
