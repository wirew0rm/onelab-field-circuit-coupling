#! /usr/bin/env python
# coding=utf-8

# Emulates a gmsh subclient returning precalculated data from a
# csv File int the format i, l
# The data is interpolated using a spline function. Axis symmetry
# is assumed so the data should contain 0 up to a reasonable high value
# for i.
# Data can be for example be obtained by running the interactive.m
# solver after a parameter sweep over i.
#
# Alexander Krimm <alex@stud.tu-darmstadt.de>
# TU-Darmstadt, 2015

import onelab
import numpy as np
import scipy.interpolate as interp

c = onelab.client(__file__)

I = c.getNumber("LICurve/I")

LI = np.loadtxt(open("nonlinear.csv", "rb"), delimiter=",")
interpolation = interp.UnivariateSpline(LI[:, 0], LI[:, 1]/1000, s=1e-6)
interpolationFlux = interp.UnivariateSpline(LI[:, 0], LI[:, 2], s=1e-6)
c.setNumber("LICurve/L", value=np.asscalar(interpolation(abs(I))))
c.setNumber("LICurve/Ldi", value=np.asscalar(interpolation(abs(I), 1)))
c.setNumber("LICurve/Flux", value=np.sign(I) *
            np.asscalar(interpolationFlux(abs(I))))
c.setNumber("LICurve/Ld", value=np.asscalar(interpolationFlux(abs(I), 1)))
