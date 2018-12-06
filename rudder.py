from aoa_lookup import getCd, getCl
from mathutils import Rotation2d, calcLift_Water, Vector2d
from units import unit


class Rudder:
    width = 6 * unit.inch
    height = 2 * unit.foot
    _angle = Rotation2d(1.0, 0.0)  # [1 0] is inline with boat hull

    angleRangeDegrees = [-45, 45]

    def __init__(self):
        self.area = self.width * self.height

    @property
    def angle(self):
        return self._angle

    @angle.setter
    def angle(self, new: Rotation2d):
        self._angle = Rotation2d.fromDegrees(
            max(self.angleRangeDegrees[0], min(self.angleRangeDegrees[-1], new.degrees))
        )

    def calcForce(self, boatRelativeWaterAngle: Rotation2d, waterVelocity) -> Vector2d:
        """
        :param boatRelativeWaterAngle: Angle of water current, relative to boat
        :param waterVelocity: Magnitude of water current, relative to boat
        :return: Boat oriented impulse vector direction + magnitude, in newtons
        """
        aoa = abs(self._angle.degrees - boatRelativeWaterAngle.degrees) * unit.degree
        cL, cD = getCl(aoa), getCd(aoa)

        lift = calcLift_Water(waterVelocity=waterVelocity, foilArea=self.area, cl=cL).to(unit.newton)

        forceDir = self._angle.normal
        torqueSign = None
        if aoa == 0:
            torqueSign = 0
        elif (self._angle.degrees - boatRelativeWaterAngle.degrees) > 0:
            torqueSign = 1
        else:
            torqueSign = -1

        lift *= torqueSign
        if (lift.m < 0):
            lift *= -1
            forceDir = forceDir.rotateBy(Rotation2d.fromDegrees(180.0))

        ret = Vector2d(0.0, 0.0)
        ret.x = forceDir.cos * lift
        ret.y = forceDir.sin * lift
        return ret


if __name__ == '__main__':
    rudder = Rudder()
    rudder.angle = Rotation2d(1.0, 0.0).rotateBy(Rotation2d.fromDegrees(10))

    print(rudder.calcForce(
        boatRelativeWaterAngle=Rotation2d(1.0, 0.0),
        waterVelocity=(1 * unit.m / unit.s)
    ))
