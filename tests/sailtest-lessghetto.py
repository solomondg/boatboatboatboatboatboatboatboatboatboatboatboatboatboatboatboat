import math
from types import SimpleNamespace

import numpy as np
import scipy as sp
from scipy.linalg import norm
import matplotlib.pyplot as plt
import pint
from pint import UnitRegistry
from sklearn.preprocessing import normalize
from collections import namedtuple

from mathutils import *

unit = UnitRegistry()

wind = SimpleNamespace(
    angle=Rotation2d(0.0, 1.0),
    speed=3 * unit.meters / unit.seconds
)

boat = SimpleNamespace(
    heading=Rotation2d(1.0, 0.0),
    velocity=0.0 * unit.meter / unit.second,
    angularVelocity=0.0 * unit.radian / unit.second,
    MOI=0.01 * unit.kg * (unit.meter ** 2)
)

mainSail = SimpleNamespace(
    width=1 * unit.foot,
    height=5 * unit.foot,
    CoD=1.25,
    relativeHeading=Rotation2d(0.0, 1.0),
    angularVelocity=0.0 * unit.radian / unit.second,
    MOI=0.18 * unit.kg * (unit.meter ** 2),
    CoF=1E-2
)
mainSail.area = mainSail.width * mainSail.height

ancillarySail = SimpleNamespace(
    width=6 * unit.inch,
    height=3 * unit.foot,
    pivotOffset=3 * unit.foot,
    CoD=1.25,
    relativeHeading=Rotation2d(1.0, 0.0)
)
ancillarySail.area = ancillarySail.width * ancillarySail.height

rudder = SimpleNamespace(
    width=4 * unit.inch,
    height=1 * unit.foot,
    CoD=1.0,
    relativeHeading=Rotation2d(1.0, 0.0),
)
rudder.area = rudder.width * rudder.height


#
#
# /----------|
# <  ==== == |
# \----------|
#   ^ main
#         ^ ancillary
#
#  ↑   ↑   ↑ wind

# Ancillary should stay [1 0], sail should move to [0 -1] or [0 1]

def calcSailPercentExposed(windAngle: Rotation2d, sailAngle: Rotation2d) -> float:
    return windAngle.rotation.dot(sailAngle.rotation)


def drag(C, u, A, Rho):
    # C is drag coefficient
    # Rho is fluid density (kg/m^3)
    # u is flow rate (m/s)
    # A is presented area (m^2)
    return 1 / 2 * Rho * (u ** 2) * C * A


def airDrag(C, u, A):
    return drag(C, u, A, 1.225 * unit.kg / (unit.meter ** 3))


def waterDrag(C, u, A):
    return drag(C, u, A, 997 * unit.kg / (unit.meter ** 3))


sailAngleSetpoint = Rotation2d.fromDegrees(-90)

csv = ""
dt = 0.1 * unit.second
simTime = 30 * unit.second
for i in np.linspace(0, simTime, int(simTime * 1 / dt)):  # 60 seconds at dt=0.01s
    mainSail.relativeHeading = \
        mainSail.relativeHeading.rotateBy(Rotation2d.fromRadians(mainSail.angularVelocity * dt))

    mainSailPresentedArea = mainSail.area * \
                            calcSailPercentExposed(wind.angle,
                                                   mainSail
                                                   .relativeHeading
                                                   .rotateBy(boat.heading)
                                                   .normal)
    ancillarySailPresentedArea = ancillarySail.area * calcSailPercentExposed(wind.angle,
                                                                             ancillarySail
                                                                             .relativeHeading
                                                                             .rotateBy(mainSail.relativeHeading)
                                                                             .rotateBy(boat.heading)
                                                                             .normal)

    mainSailAoa = math.fabs(wind.angle.degrees - mainSail.relativeHeading.rotateBy(boat.heading).degrees)
    print(mainSailAoa)

    break







    mainForce = airDrag(mainSail.CoD, wind.speed, mainSailPresentedArea).to(unit.lbf)
    ancillaryForce = airDrag(ancillarySail.CoD, wind.speed, ancillarySailPresentedArea).to(unit.lbf)

    mainSailTorque = (ancillaryForce * ancillarySail.pivotOffset).to(unit.meter * unit.newton)
    mainSail.angularVelocity = mainSail.angularVelocity + (mainSailTorque / mainSail.MOI).to(
        unit.radian / unit.second ** 2) * dt - mainSail.angularVelocity * mainSail.CoF
    mainSailAngularAccel = (mainSailTorque / mainSail.MOI).to(unit.radian / unit.second ** 2)

    # print("Main presented area: ", mainSailPresentedArea)
    # print("Ancillary presented area: ", ancillarySailPresentedArea)
    # print("Main force: ", mainForce)
    # print("Ancillary force: ", ancillaryForce)

    # print("Sail angular velocity:", mainSail.angularVelocity)
    # print("Sail angular acceleration:", mainSailAngularAccel)

    # print("Sail heading:", str(mainSail.relativeHeading.degrees * unit.degrees))

    # ancillarySail.relativeHeading = Rotation2d.fromDegrees(
    #    min(45, max((Rotation2d(0, 1).degrees - mainSail.relativeHeading.degrees) * -.1, -45)))

    csv += str(i) + "," + str(mainSail.relativeHeading.degrees) + "," + str(mainForce.m) + "," + str(
        math.degrees(mainSail.angularVelocity.m)) + "," + str(ancillarySail.relativeHeading.degrees) + \
           "," + str(Rotation2d(0, -1).degrees - mainSail.relativeHeading.degrees) + "\n"

    error = sailAngleSetpoint.degrees - mainSail.relativeHeading.degrees
    ancillarySail.relativeHeading = Rotation2d.fromDegrees(
        min(45, max(error * .9 + (mainSail.angularVelocity * unit.dt).m * 5.3, -45))
    )

    break

with open('plot.csv', 'w') as f:
    f.write(csv)
