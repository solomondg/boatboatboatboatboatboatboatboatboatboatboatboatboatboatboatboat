from math import cos, sin, atan2

import numpy as np
import scipy as sp
from scipy.linalg import norm
import matplotlib.pyplot as plt
import pint
from pint import UnitRegistry
from sklearn.preprocessing import normalize
from collections import namedtuple

Vec2d = namedtuple('Vec2d', ['x', 'y'])

ureg = UnitRegistry()
rad = ureg.radian
deg = ureg.degree
m_s = ureg.meter / ureg.second
inch = ureg.inch

windDir = np.asarray([1.3, .32])
windDir = np.asarray([0, 1])  # Should result in a 100% vertical sail?
windStrength = 10 * m_s
windVec = (windDir / norm(windDir)) * windStrength

windVec = Vec2d(x=windVec[0], y=windVec[1])

sailC_x = 0.05
sailC_y = 1.25

mainSailWidth = 1 * ureg.foot
mainSailHeight = 5 * ureg.foot
mainSailThickness = .25 * inch
ancillarySailWidth = 6 * inch
ancillarySailHeight = 3 * ureg.foot
ancillarySailThickness = .25 * inch
ancillarySailPivotOffset = 3 * ureg.foot

mainSailArea = mainSailWidth * mainSailHeight
ancillarySailArea = ancillarySailWidth * ancillarySailHeight

sailAngle = 0 * deg
ancillarySailActuateAngle = 0 * deg

angVel = 0 * ureg.rad / ureg.second

sailMOI = 0.005

t0 = 0
tEnd = 1000
timeStep = 0.01
times = [x * timeStep for x in range(t0, int(tEnd / timeStep))]


def drag(C, u, A, Rho=1.225 * ureg.kg / (ureg.meter ** 3)):
    # C is drag coefficient
    # Rho is fluid density (kg/m^3)
    # u is flow rate (m/s)
    # A is area (m^2)
    return 1 / 2 * Rho * (u ** 2) * C * A


# kg/m^3 * m^2/s^2 * m^2
# (kg*m^4)/(m^3*s^2)
# kg*m/s^2


def presentedArea(width, height, thick, angle):
    # angle of 0, 180 results in w*h
    # angle of 90, -90 results in t*h
    angle += 90 * deg
    # cos(angle) = proj/width
    # proj = cos(angle)*width
    return abs(cos(angle.to(rad)) * width * height) + abs(sin(angle.to(rad)) * thick * height)


def effectiveCod(angle, c_x, c_y):
    angle += 90 * deg
    return abs(cos(angle.to(rad)) * c_x) + abs(sin(angle.to(rad)) * c_y)


times = [0]
for t in times:
    sailAngle = 30 * deg
    effectiveWindAngle = sailAngle - atan2(windVec.y.magnitude, windVec.x.magnitude) * rad

    sailArea = presentedArea(mainSailWidth, mainSailHeight, mainSailThickness, effectiveWindAngle)
    ancillaryArea = presentedArea(ancillarySailWidth, ancillarySailHeight, ancillarySailThickness,
                                  effectiveWindAngle + ancillarySailActuateAngle)

    print("Theoretical main:", mainSailWidth * mainSailHeight)
    print("Theoretical ancillary:", ancillarySailWidth * ancillarySailHeight)
    print("Actual main:", sailArea)
    print("Actual ancillary:", ancillaryArea)
    sailCod = effectiveCod(sailAngle, sailC_x, sailC_y)
    ancillaryCod = effectiveCod(sailAngle + ancillarySailActuateAngle, sailC_x, sailC_y)
    print("Sail effective CoD:", sailCod)
    print("Ancillary effective CoD:", ancillaryCod)

    # Calculate forces on boat
    mainSailAbsForce = drag(sailCod, windStrength, sailArea.to(ureg.meter ** 2)).to('newton')
    ancillarySailAbsForce = drag(ancillaryCod, windStrength, ancillaryArea.to(ureg.meter ** 2)).to('newton')
    print("Force on main sail:", mainSailAbsForce)
    print("Force on ancillary sail:", ancillarySailAbsForce)
