import math
from types import SimpleNamespace
from typing import Tuple

import numpy as np
import scipy as sp
from scipy.linalg import norm
import matplotlib.pyplot as plt
import pint
from pint import UnitRegistry
from sklearn.preprocessing import normalize
from collections import namedtuple

from main import unit
from mathutils import Rotation2d, Vector2d, airDrag


# TODO: Make actual lift calculations instead of doing 100% angle of attack
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
        MOI=0.18 * unit.kg * (unit.meter ** 2),
        CoF=1E-2
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

    def setTrailingEdgeAngle(self, newAngle: Rotation2d):
        assert self.ancillarySail.minAngle < newAngle.degrees < self.ancillarySail.maxAngle
        self.ancillarySail.relativeHeading = newAngle

    def update(self, dt, windAngle: Rotation2d, windStrength, boatHeading: Rotation2d) -> Vector2d:
        """
        :param dt: Update delta time, in seconds
        :param windAngle: Wind angle relative to boat ( [0 1] is a headwind)
        :param windStrength: Wind strength, in m/s
        :return: Boat oriented impulse vector direction + magnitude, in newtons
        """

        # Calculate effective area of the side (pushy bit) of the sail that the wind's hitting
        mainSailSidePresentedArea = self.mainSail.sideArea * \
                                    _calcSailPercentExposed(windAngle,
                                                            self.mainSail
                                                            .relativeHeading
                                                            .rotateBy(boatHeading)
                                                            .normal)
        ancillarySailSidePresentedArea = self.ancillarySail.sideArea * \
                                         _calcSailPercentExposed(windAngle,
                                                                 self.ancillarySail
                                                                 .relativeHeading
                                                                 .rotateBy(self.mainSail.relativeHeading)
                                                                 .rotateBy(boatHeading)
                                                                 .normal)
        # Calculate effective area of the front (useless bit) of the sail that the wind's hitting
        mainSailFrontPresentedArea = self.mainSail.frontArea * \
                                     (1 - _calcSailPercentExposed(windAngle,
                                                                  self.mainSail
                                                                  .relativeHeading
                                                                  .rotateBy(boatHeading)
                                                                  .normal)
                                      )
        ancillarySailFrontPresentedArea = self.ancillarySail.frontArea * \
                                          (1 - _calcSailPercentExposed(windAngle,
                                                                       self.ancillarySail
                                                                       .relativeHeading
                                                                       .rotateBy(self.mainSail.relativeHeading)
                                                                       .rotateBy(boatHeading)
                                                                       .normal)
                                           )

        mainSailSideForce = airDrag(C=self.mainSail.sideCoD, u=windStrength, A=mainSailSidePresentedArea)
        ancillarySailSideForce = airDrag(C=self.ancillarySail.sideCoD, u=windStrength, A=ancillarySailSidePresentedArea)
        mainSailFrontForce = airDrag(C=self.mainSail.sideCoD, u=windStrength, A=mainSailFrontPresentedArea)
        ancillarySailFrontForce = airDrag(C=self.ancillarySail.sideCoD, u=windStrength, A=ancillarySailFrontPresentedArea)





def _calcSailPercentExposed(windAngle: Rotation2d, sailAngle: Rotation2d) -> float:
    return windAngle.rotation.dot(sailAngle.rotation)
