import numpy as np
import xarray as xr
from matplotlib import pyplot as plt
import pandas as pd
from matplotlib import ticker
import os, sys
from scipy.optimize import minimize
import WetBulb

# Generate temperature range (in degrees C)
N_steps = 500

min_temp = 25
max_temp = 50

temp = np.linspace(min_temp,max_temp,N_steps)

# Generate specific humidity range (in kg/kg)
min_q = 0
max_q = 0.02

q = np.linspace(min_q, max_q, N_steps)

# Calculate wet bulb temperature grid
pres = 1000*100 # set general surface pressure to 1000 hPa

# Make gridded arrays
temp_grid, q_grid = np.meshgrid(temp, q) 
pres_grid = np.ones_like(temp_grid)*pres
    
# For users who wish to derive stickiness based on a different humid heat metric which combines temperature and humidity,
# calculate that metric here on the temperature and humidity grid calculated above.
Twb,Teq,epott = WetBulb.WetBulb(temp_grid, pres_grid, q_grid, HumidityMode = 0)

def stickiness(v):
    
    # Calculate stick using input vector
    C_00 = v[0]
    C_10 = v[1]
    C_01 = v[2]
    C_11 = v[3]
    C_20 = v[4]
    C_02 = v[5]
    C_21 = v[6]
    C_12 = v[7]
    C_22 = v[8]
    C_30 = v[9]
    C_03 = v[10]
    C_31 = v[11]
    C_13 = v[12]
    C_32 = v[13]
    C_23 = v[14]
    C_33 = v[15]

    stick = (C_00 + C_10*q_grid + C_01*temp_grid + C_11*q_grid*temp_grid + 
             C_20*q_grid**2 + C_02*temp_grid**2 + C_21*(q_grid**2)*temp_grid + 
             C_12*(temp_grid**2)*q_grid + C_22*(q_grid**2)*(temp_grid**2) + 
             C_30*q_grid**3 + C_03*temp_grid**3 + C_31*(q_grid**3)*temp_grid + 
             C_13*(temp_grid**3)*q_grid + C_32*(q_grid**3)*(temp_grid**2) + 
             C_23*(temp_grid**3)*(q_grid**2) + C_33*(q_grid**3)*(temp_grid**3))
    
    lam1 = 0.8
    lam2 = 1 - lam1
    
    dstick_dT, dstick_dq = np.gradient(stick)
    
    # Users employing a different humid heat metric can change 'Twb' to their metric of choice
    dTW_dT, dTW_dq = np.gradient(Twb)
    
    int1 = (dstick_dq / dTW_dq + dstick_dT/dTW_dT)**2
    int2 = (dstick_dq/dTW_dq - 1)**2
    
    dq = (max_q - min_q)/N_steps
    dT = (max_temp - min_temp)/N_steps
    
    eps = (lam1*np.sum(int1)*dq*dT + lam2*np.sum(int2)*dq*dT)**0.5
    
    print('error = ' + eps)
    
    return eps

guess = [3.92578019e+00, -7.81268471e+02,  2.36368404e-01, -5.09311865e-01,
         7.02523267e+03,  2.23047656e-03, -2.44331804e+01, -7.54293880e-02,
         4.54184507e-01, -4.15274027e+02, -2.98769414e-05, -2.16736557e+01, 
         1.44044908e-03, -1.17047276e+01, -2.15204251e-02,  3.09008029e-01]

result = minimize(stickiness,guess, method='nelder-mead', options = {'xatol': 1e-08,'disp': True, 'maxfev': 100000})
print('The final array that minimizes mean squared error is:'+ str(result.x))