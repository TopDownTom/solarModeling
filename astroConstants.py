import math
import numpy as np

exp = math.exp
e = math.exp
sqrt = np.sqrt
pi = math.pi
G = 6.67e-8 # cm^3 g^-1 s^-2
c = 2.99792458e10 # cm s^-1
R = 6.96e10 # cm
mStar = 2e33 # g
lStar = 3.9e26 # W
k = 1.380658e-16 # erg K^-1
rhoc = 150 # g cm^-3
mp = 1.6726231e-24 # g
mn = 1.6749286e-24 # g
mH = 1.6733e-24 # g
alpha = R/5 # number

# Opacity - Caroll & Ostlie ch. 9
X = .70 # % hydrogen abundance
Y = .28 # % helium abundance
Z = .02 # % metals abundance
Xa = .7
Xb = .7
Aa = 1
Ab = 1
gamma = 5/3 # ratio of specific heats for ionized hydrogen
gbf = 1 # Gaunt factor Caroll & Ostlie p. 250 
gff = 1 # Gaunt factor
t = 100 # guillotine factor Caroll and & Ostlie p. 250
muI = 0.62 # Caroll & Ostlie p. 293
muN = 1.30 # Caroll & Ostlie p. 293
mu = 0.84 
mup = mp/2

# Energy Production - Caroll & Ostlie ch. 10
a = 7.565767e-15 # Caroll & Ostlie p. 234   [ erg cm^-3 K^-4 ]
S0 = 4e-43 # Dan Maoz p.54
Q = 26.73e6 # eV p-p chain release Dan Maoz p.54
Eg = .500e6 # eV for 2 protons Dan Maoz p.54
e0pp = 1.08e-5 # erg cm^3 s^-1 g^-2 # Caroll & Ostlie p. 311
e0cno = 8.24e-26 # erg cm^3 s^-1 g^-2 # Caroll & Ostlie p. 312
muI = .62 # Caroll & Ostlie 293

