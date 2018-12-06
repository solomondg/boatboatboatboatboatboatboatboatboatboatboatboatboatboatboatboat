import math
from types import SimpleNamespace
from typing import Tuple, List

import numpy as np
import scipy as sp
from scipy.linalg import norm
import matplotlib.pyplot as plt
import pint
from pint import UnitRegistry
from sklearn.preprocessing import normalize
from collections import namedtuple

from aoa_lookup import getCl, getCd
from units import unit
from mathutils import Rotation2d, Vector2d, airDrag, calcLift_Air


# CCW coordinate system
#         -y
#     /----------|
# +x <   ==== == |
#     \----------|
#         ^ main sail
#             ^ ancillary sail

# Mainsail relativeHeading is relative to boar
# [1 0] is pointing forwards
# [0 1] is pointing right
#

# Ancillary sail relative heading is relative to main sail
# (main,ancillary)
# [1 0], [1 0] is both main and ancillary pointing forwards
# [0.707 0.707] [1 0] is both main and ancillary in line but pointing 45 degrees to port (aka ancillary hanging over starboard side)
# [0.707 0.707] [0.707 -0.707] is main at 45 degrees to port but ancillary parallel with boat
#


# TODO: Calc drag, make it act opposite to air flow direction
class Sail:
    # NACA 0012 airfoil
    mainSail = SimpleNamespace(
        width=1 * unit.foot,
        height=5 * unit.foot,
        thicc=4 * unit.inch,
        sideCoD=1.25,
        frontCoD=0.6,
        relativeHeading=Rotation2d(1.0, 0.0),
        angularVelocity=0.0 * unit.radian / unit.second,
        MOI=0.50 * unit.kg * (unit.meter ** 2),
        ViscFriction=1e-0 * unit.newton * unit.meter * unit.second
    )

    # NACA 0012 airfoil
    ancillarySail = SimpleNamespace(
        width=6 * unit.inch,
        height=3 * unit.foot,
        thicc=4 * unit.inch,
        pivotOffset=3 * unit.foot,
        sideCoD=1.25,
        frontCoD=0.6,
        relativeHeading=Rotation2d(1.0, 0.0),
        maxAngle=Rotation2d.fromDegrees(45),
        minAngle=Rotation2d.fromDegrees(-45)
    )

    def __init__(self):
        self.mainSail.sideArea = self.mainSail.width * self.mainSail.height
        self.ancillarySail.sideArea = self.ancillarySail.width * self.ancillarySail.height
        self.mainSail.frontArea = self.mainSail.thicc * self.mainSail.height
        self.ancillarySail.frontArea = self.ancillarySail.thicc * self.ancillarySail.height

    @property
    def trailingEdgeAngle(self):
        return self.ancillarySail.relativeHeading

    @trailingEdgeAngle.setter
    def trailingEdgeAngle(self, newAngle: Rotation2d):
        assert self.ancillarySail.minAngle < newAngle.degrees < self.ancillarySail.maxAngle
        self.ancillarySail.relativeHeading = newAngle

    def _calcSailPercentExposed(self, boatRelativeWindAngle: Rotation2d, sailAngle: Rotation2d) -> float:
        """
        :param boatRelativeWindAngle: Angle of wind relative to boat
        :param sailAngle: Angle of sail relative to boat
        :return: Percent of sail that's exposed to wind (as the angle of attack becomes smaller, less is exposed)
        """
        return boatRelativeWindAngle.rotation.dot(sailAngle.rotation)

    def update(self, dt, boatRelativeWindAngle: Rotation2d, boatRelativeWindSpeed) -> Vector2d:
        """
        :param dt: Update delta time, in seconds
        :param boatRelativeWindAngle: Wind angle relative to boat ([-1 0] is a headwind)
        :param boatRelativeWindSpeed: Wind strength, in m/s
        :return: Boat oriented impulse vector direction + magnitude, in newtons
        """

        # Calculate angle of attack of sails to get the coefficient of lift and drag
        mainSailAoA = abs(self.mainSail.relativeHeading
                          # .rotateBy(boatHeading)
                          .degrees - boatRelativeWindAngle.degrees) * unit.degree
        ancillarySailAoA = abs(self.ancillarySail.relativeHeading
                               .rotateBy(self.mainSail.relativeHeading)
                               # .rotateBy(boatHeading)
                               .degrees - boatRelativeWindAngle.degrees) * unit.degree
        mainCl, mainCd = getCl(mainSailAoA), getCd(mainSailAoA)
        ancillaryCl, ancillaryCd = getCl(ancillarySailAoA), getCd(ancillarySailAoA)

        # Calculate effective area of the side (pushy bit) of the sail that the wind's hitting. Not necessary now, might be later (for drag).
        # mainSailSidePresentedArea = self.mainSail.sideArea * \
        #                            self._calcSailPercentExposed(boatRelativeWindAngle,
        #                                                         self.mainSail
        #                                                         .relativeHeading
        #                                                         # .rotateBy(boatHeading)
        #                                                         .normal)
        # ancillarySailSidePresentedArea = self.ancillarySail.sideArea * \
        #                                 self._calcSailPercentExposed(boatRelativeWindAngle,
        #                                                              self.ancillarySail
        #                                                              .relativeHeading
        #                                                              .rotateBy(self.mainSail.relativeHeading)
        #                                                              # .rotateBy(boatHeading)
        #                                                              .normal)

        mainSailLift = calcLift_Air(airVelocity=boatRelativeWindSpeed,
                                    wingArea=self.mainSail.sideArea,
                                    cl=mainCl)
        ancillarySailLift = calcLift_Air(airVelocity=boatRelativeWindSpeed,
                                         wingArea=self.ancillarySail.sideArea,
                                         cl=ancillaryCl)

        # print(mainSailLift.to(unit.newton))
        # print(ancillarySailLift.to(unit.newton))

        # AoA is side invarient, now we need to figure out which way our force is actually applied
        mainSailForce = mainSailLift.to(unit.newton)
        forceOnBoatDirection = None
        if mainSailAoA == 0:
            forceOnBoatDirection = self.mainSail.relativeHeading
            mainSailForce = 0 * unit.newton
        elif self.mainSail.relativeHeading.degrees - boatRelativeWindAngle.degrees > 0:
            forceOnBoatDirection = self.mainSail.relativeHeading.normal
        else:
            forceOnBoatDirection = self.mainSail.relativeHeading.normal.rotateBy(Rotation2d.fromDegrees(180.0))

        # Same deal, but just a torque multiplier this time, since the sail assembly can only pivot and not translate xy
        ancillarySailForce = ancillarySailLift.to(unit.newton)
        torqueSign = None
        if ancillarySailAoA == 0:
            torqueSign = 0
        elif self.ancillarySail.relativeHeading \
                .rotateBy(self.mainSail.relativeHeading).degrees - boatRelativeWindAngle.degrees > 0:
            torqueSign = -1
        else:
            torqueSign = 1

        # Now we have to calculate the torque that the ancillary sail puts out on the sail assembly
        torque = (self.ancillarySail.pivotOffset * ancillarySailForce * torqueSign).to(unit.newton * unit.meter)
        angularAcceleration = (torque - self.mainSail.angularVelocity * self.mainSail.ViscFriction) / self.mainSail.MOI
        # t' = -tb/J
        self.mainSail.angularVelocity += angularAcceleration * dt

        self.mainSail.relativeHeading = self.mainSail.relativeHeading \
            .rotateBy(
            Rotation2d.fromDegrees(
                (self.mainSail.angularVelocity * dt).to(unit.degrees).m
            )
        )

        ret = Vector2d(0.0, 0.0)
        ret.x = forceOnBoatDirection.cos * mainSailForce
        ret.y = forceOnBoatDirection.sin * mainSailForce
        return ret


if __name__ == '__main__':
    sail = Sail()
    sail.mainSail.relativeHeading = Rotation2d.fromDegrees(0)
    sail.ancillarySail.relativeHeading = Rotation2d.fromDegrees(30)
    X, Y = [], []
    dt = 0.01

    for i in range(1000):
        try:
            sail.update(dt * unit.second,
                        boatRelativeWindAngle=Rotation2d(1.0, 0.0),
                        boatRelativeWindSpeed=5 * unit.meter / unit.second,
                        )
            X.append(i * dt)
            Y.append(sail.mainSail.relativeHeading.degrees)
        except:
            break

    import matplotlib.pyplot as plt

    plt.plot(X, Y)
    plt.show()
