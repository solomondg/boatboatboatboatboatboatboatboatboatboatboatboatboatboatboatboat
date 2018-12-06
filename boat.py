import math
from types import SimpleNamespace
import numpy as np
import scipy as sp

from rudder import Rudder
from sail import Sail
from units import unit
from mathutils import Rotation2d, Vector2d, Translation2d, RigidTransform2d
from wind import Wind


class Boat:
    orientation = RigidTransform2d(
        Translation2d(0.0 * unit.meter, 0.0 * unit.meter),
        Rotation2d(1.0, 0.0)
    )

    # Boat relative ([1.0, 0.0] is always 1 m/s forwards)
    velocity = Vector2d(0.0 * unit.meter / unit.second, 0.0 * unit.meter / unit.second)
    angularVelocity = 0 * unit.radian / unit.second

    MOI = 10 * unit.kg * unit.meter ** 2

    fwdCrossSectionArea = 1 * unit.foot * 2 * unit.foot
    sideCrossSectionArea = 8 * unit.foot * 2 * unit.foot

    sail = Sail()
    rudder = Rudder()

    forwardCd = 0.15
    sideCd = 1.2
    viscFriction = 1e1 * unit.newton * unit.meter * unit.second

    def __init__(self):
        pass

    @property
    def spankerAngle(self):
        return self.sail.trailingEdgeAngle

    @spankerAngle.setter
    def spankerAngle(self, newAngle: Rotation2d):
        self.sail.trailingEdgeAngle = newAngle

    @property
    def rudderAngle(self):
        return self.rudder.angle

    @rudderAngle.setter
    def rudderAngle(self, newAngle: Rotation2d):
        self.rudder.angle = newAngle

    def update(self, *, dt):
        forceFromSail = self.sail.update(
            dt=dt,
            boatRelativeWindAngle=Wind.angle.rotateBy(self.orientation.rot.rotation)
        )
