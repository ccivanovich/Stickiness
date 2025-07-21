# Calculates stickiness based off of the derivation using wet bulb
# temperature as described in Ivanovich et al. 2024. Requires 
# inputs of dry bulb temperature (C) and specific humidity (kg/kg).
# Final units of stickiness are in degrees C.
#
# Reference: Ivanovich, C., Sobel, A., Horton, R., Raymond, C., 2024.
#          Stickiness: A New Variable to Characterize the Temperature
#          and Humidity Contributions toward Humid Heat, Journal of the
#          Atmospheric Sciences.

def stickiness(T,q):
    v = [1.199557787089589, -7.75269363e+02, 3.02416275e-01, 3.08627625e+00,
              7.74095672e+03, 1.78241619e-03, -2.38012318e+02, -4.98714576e-01,
              2.40170949e+01, 7.18600120e+03, -2.75682633e-05, 4.29814844e+02,
              7.01883383e-03, -2.83672170e+02, -3.67109506e-01, 5.09430371e+00]
    
    stickiness = -1*(v[0] + v[1]*q + v[2]*T + v[3]*q*T + 
             v[4]*q**2 + v[5]*T**2 + v[6]*(q**2)*T + 
             v[7]*(T**2)*q + v[8]*(q**2)*(T**2) + 
             v[9]*q**3 + v[10]*T**3 + v[11]*(q**3)*T + 
             v[12]*(T**3)*q + v[13]*(q**3)*(T**2) + 
             v[14]*(T**3)*(q**2) + v[15]*(q**3)*(T**3))
    
    return stickiness

