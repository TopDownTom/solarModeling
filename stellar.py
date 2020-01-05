# Import the necessary modules and files
# Imported files (astroConstants and rk4Function) must be in same director as this file, else an absolute path must be given
import time, math
from array import *
import numpy as np
import matplotlib.pyplot as plt
from astroConstants import *
from rk4Function import *

# Define each derivative
def dP(y,r):
    dP = -G*rho*y[1]/r**2
    return dP
def dM(y,r):
    dM = 4*pi*rho*r**2
    return dM
def dL(y,r):
    dL = 4*pi*rho*r**2*epsilon
    return dL
def dTrad(y,r):
    dTrad =  (- 3*kappaBar*rho*y[2] ) / ( 4*a*c*((y[3])**3) * 4*pi*(r**2) ) 
    return dTrad
def dTconv(y,r):
    dTconv =  ( (1/gamma) - 1 ) * ( mu*mH/k * G*y[1]/r**2 ) 
    return dTconv

# initial conditions
r = int(1e5) # beginning radius at 10^5 cm
y = np.array([2.34e17,0.001,0.001,1.527e7]) # Press,Mass,Lum,Temp at core

# Calculate an initial density
rho = mu*mH/k * y[0]/y[3]

# Depending which energy density you use must be changed in loop below as well. 
# Unsure at this time which is the better expression for energy density, however results were closest with what is currently un-commented.

#epsilonPP = .241e4 * rho* X**2 * (y[3]/1e6)**(-2/3) * np.exp( -33.8*(y[3]/1e6)**(-1/3) )
#epsilonCNO = 8.67e24 * rho * X*Z * (y[3]/1e6)**(-2/3) * np.exp( -152.28*(y[3]/1e6)**(-1/3) )
epsilonPP = e0pp * rho * X**2 * (y[3]/1e6)**4
epsilonCNO = e0cno * rho * X * Z * (y[3]/1e6)**(19.9)
epsilon = epsilonPP + epsilonCNO
#epsilon = 0.045*rho*.5**2/np.sqrt(mup)*(y[3]**(-2/3))*np.exp(-(3.395e3)/(y[3]**(1/3)))


# initialize lists
pressureN = [y[0]]
massN = [y[1]]
luminosityN = [y[2]]
dLumN = [dL(y,r)]
temperatureN = [y[3]]
densityN = [rho]
epsilonN = [epsilon]
rN = [r]

# Max num steps and step size
N = int(1e11)
step = int(1e7)

# solve each function
for i in range(r,N):

# Solve for energy production

#    epsilonPP = .241e4 * rho* X**2 * (y[3]/1e6)**(-2/3) * np.exp( -33.8*(y[3]/1e6)**(-1/3) )
#    epsilonCNO = 8.67e24 * rho * X*Z * (y[3]/1e6)**(-2/3) * np.exp( -152.28*(y[3]/1e6)**(-1/3) )
# As power laws centered at 1.5e7 K
    epsilonPP = e0pp * rho * X**2 * (y[3]/1e6)**4
    epsilonCNO = e0cno * rho * X * Z * (y[3]/1e6)**(19.9)
#    epsilon = epsilonPP + epsilonCNO
    epsilon = .045*rho*.5**2/np.sqrt(mup)*(y[3]**(-2/3))*np.exp(-(3.395e3)/(y[3]**(1/3)))

# Update density at each step
    rho = mu*mH/k * y[0]/y[3]

# Opacity
    if ( 3e3 <= y[3] <= 6e3 ) and ( 1e-7 <= rho*1e3 <= 1e-2 ): 
        kappaH = 7.9e-35*(Z/.02)*rho**(1/2)*y[3]**9
    else:
        kappaH = 0
    kappaES = 0.2*(1+X)
    kappaBF = 4.34e22*(gbf/t)*Z*(1+X)*rho/(y[3]**(3.5))
    kappaFF = 3.68e19*(gff)*(1-Z)*(1+X)*rho/(y[3]**(3.5))
    kappaBar = np.average(kappaES + kappaH + kappaBF + kappaFF)

# The actual integration occurs for each of the stellar structure equations here
    press = rk4(dP,y,r,step)
    y[0] = press[0]

    mass = rk4(dM,y,r,step)
    y[1] = mass[1]

    lum = rk4(dL,y,r,step)
    y[2] = lum[2]

# Criteria for whether the stars transfers energy radiatively or convectively
# This is an area for potential investigation/tweaking
    convCriteria = np.log(y[0]) / np.log(y[3])

    if convCriteria > (gamma/(gamma-1)):
        temp = rk4(dTrad,y,r,step)
        y[3] = temp[3]
    else:
        temp = rk4(dTconv,y,r,step)
        y[3] = temp[3]

# Stop either at the observed radius of the Sun or if the pressure (y[0]) goes to zero
    if r > int(6.95e10):
        break
    if y[0] < 0:
        break

# Increment the values of each derivative
    y = y
    r = r + step

# Store values of each
    pressureN.append(y[0])
    massN.append(y[1])
    luminosityN.append(y[2])
    dLumN.append(dL(y,r))
    temperatureN.append(temp[3])
    densityN.append(rho)
    epsilonN.append(epsilon)
    rN.append(r)

#### Post-loop data processing and plotting

# Print the final values for the variables of interest wherever the code stops
print()
print()
print("Final Values:")
print("-----------------")
print("Radius: {:.2e} cm".format(r))
#print("Pressure: {:.2e} Ba".format(y[0]))
print("Rho: {:.2e} kg m^-3".format(rho*1e3))
print("Pressure: {:.2e} N m^-2".format(y[-1]/10))
print("Mass: {:.2e} kg".format(y[1]/1e3))
print("Luminosity: {:.2e} W".format(y[2]/1e7))
print("Max dL: {:.2e} W m^-1".format(max(dLumN)*1e-5))
print("Temperature: {:.2e} K".format(y[3]))
print("-----------------")
print()
print()

# Use this section to print out some reference values at the radius the code stops at
fracStar = r/R
print("Expected Vals at this Radius:")
print("-----------------")
print("Mass: {:.2e} kg".format(fracStar*mStar/1e3))
print("Lum: {:.2e} W".format(fracStar*lStar))
print("-----------------")
print()
print()

## This section converts several variables into MKS units before plotting
maxLum = [ i / max(luminosityN) for i in luminosityN ] # scale the max lum if we want the luminosity plotted from 0->1
massKg = [ i / 1e3 for i in massN ] # g -> kg
densityKg = [ i * 1e3 for i in densityN ] # g cm^-3 -> kg m^-3
scaledLum = [ i * 1e-7 for i in luminosityN ] # erg s^-1 -> W
epsilonW = [ i * 1e-4 for i in epsilonN ] # erg s^-1 g^-1 -> W kg^-1
solarRadii = [ i / int(6.95e10) for i in rN ] # x-axis in fractions of the solar radius
stellarRadii = [ i / rN[-1] for i in rN ] # x-axis in fractions of the achieved stellar radius


#### Plotting

# Font size definitions
SMALL_SIZE = 26
MEDIUM_SIZE = 28
BIGGER_SIZE = 30

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)   # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# Uncomment which variables you want to plot

#plt.title("Mass")
#plt.plot(solarRadii,massKg)
#plt.text(.2,4.8e29,"Expected: 5.6e29 kg", fontsize = 22)
#plt.text(.2,4e29,"Actual: 8.1e29 kg", fontsize = 22)
#plt.ylabel("Mass [kg]")
#plt.xlabel("Radius [r/Rstar]")
#plt.show()

#plt.title("Density")
#plt.plot(solarRadii,densityN)
#plt.ylabel("Density [kg/m^3]")
#plt.xlabel("Radius [r/Rstar]") 
#plt.show()

#plt.title("Luminosity")
#plt.plot(solarRadii,scaledLum)
#plt.text(.18,2.5e26,"Expected: 1.1e26 W", fontsize=22)
#plt.text(.18,2.0e26,"Actual: 3.5e26 W", fontsize=22)
#plt.ylabel("Luminosity [W]")
#plt.xlabel("Radius [r/Rstar]")
#plt.show()

#plt.title(r'$\Delta$'"Luminosity")
#plt.plot(solarRadii,scaleddLum,label="scaled dL")
#plt.text(.18,5e18,"Expected: ~4.2e18 W/m", fontsize=22)
#plt.text(.18,4.2e18,"Actual: 6.2e18 W/m",fontsize=22)
#plt.ylabel("dL/dR [Watts/m]")
#plt.xlabel("Radius [r/Rstar]")
#plt.show()

#plt.title("Pressure")
#plt.ylabel("Pressure [N/m^2]")
#plt.xlabel("Radius [r/Rstar]")
#plt.plot(solarRadii,pressureN)
#plt.show()

#plt.title("Energy Density   " r'$\epsilon$')
#plt.ylabel("Energy Density [W/kg]")
#plt.xlabel("Radius [r/Rstar]")
#plt.plot(solarRadii,epsilonW)
#plt.show()

#plt.title("Temperature")
#plt.plot(solarRadii,temperatureN)
#plt.text(.13,1.4e7,"Expected surface temp: ~5800 K", fontsize=22)
#plt.text(.13,1.2e7, "Actual surface temp: 3290 K", fontsize=22)
#plt.ylabel("Temperature [K]")
#plt.xlabel("Radius [r/Rstar]")
#plt.show()
