#! /usr/bin/env python
# Visualize different Inductances that can be obtained from
# the data the onelab inductor model provides.
# Alexander Krimm, <alex@wirew0rm.de>, 2015

import numpy as np
import matplotlib.pyplot as pl
import scipy.interpolate as interp
LI = np.loadtxt(open("nonlinear.csv", "rb"), delimiter=",")
# i = LI[:-1, 0]
# l = LI[:-1, 1]/1000
# f = LI[:-1, 2]
i = LI[:, 0]
l = LI[:, 1]/1000
f = LI[:, 2]
interpolation = interp.UnivariateSpline(i, l, s=1e-6)
ldi = interpolation(i, 1)
interpolationFlux = interp.UnivariateSpline(i, f, s=1e-6)
fdi = interpolationFlux(i, 1)
lchord = f/i
ld2 = l + ldi*i
# Brauer Kennlinie
k1 = 20
k2 = 1
k3 = 0.003
lbr = 1.78/(k1+k2*np.exp(k3*i**2))+0.025
interpolationbr = interp.UnivariateSpline(i, lbr, s=1e-6)
ldbr = lbr - i*1.78*2*k2*k3*i*np.exp(k3*i**2)/(k1+k2*np.exp(k3*i**2))**2
ldbr2 = lbr + interpolationbr(i, 1)*i
fbr = lbr*i
pl.subplot(2, 1, 1)
pl.plot(i, l, label='L (from onelab)')
pl.plot(i, interpolation(i), linestyle='dashed', label='Lspline')
pl.plot(i, ld2, label='Ld = L + dLdi*i')
pl.plot(i, fdi, label='Ld = dPhi/di')
pl.plot(i, lchord, label='Lchord = Phi/i')
pl.plot(i, lbr, label='L_Brauer')
pl.plot(i, ldbr, label='Ld_Brauer')
pl.plot(i, ldbr2, linestyle='dashed', label='Ld_Brauer_interpolationderive')
pl.ylabel('Inductance [H]')
pl.legend()
pl.subplot(2, 1, 2)
pl.plot(i, f, label='Flux')
pl.plot(i, fbr, label='Flux = L_Br*i')
pl.plot(i, interpolationFlux(i), linestyle='dashed', label='Flux Spline')
pl.xlabel('Current [A]')
pl.ylabel('Flux [Wb]')
pl.legend()
pl.show()
