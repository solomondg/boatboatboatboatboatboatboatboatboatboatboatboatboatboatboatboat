import numpy as np
import math

from units import unit

kEpsilon = 1e-8


class Point2d:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __str__(self):
        return "x: " + str(self.x) + " y: " + str(self.y)

    def __repr__(self):
        return self.__str__()

    def as_list(self):
        return [self.x, self.y]

    def as_vector(self):
        return Vector2d(self.x, self.y)

    def __sub__(self, other):
        return Point2d(self.x - other.x, self.y - other.y)


class Vector2d:
    @staticmethod
    def from_ndarray(nd: np.ndarray):
        assert nd.shape == (2,)
        return Vector2d(nd[0], nd[1])

    def __init__(self, x, y):
        self.x = x
        self.y = y

    def as_ndarray(self):
        return np.asarray([self.x, self.y])

    def normalize(self):
        return Vector2d.from_ndarray(self.as_ndarray() / (np.linalg.norm(self.as_ndarray())))

    def __str__(self):
        return "[" + str(self.x) + " " + str(self.y) + "]"

    def as_point(self):
        return Point2d(self.x, self.y)

    def dot(self, other) -> float:
        return float(np.dot(self.as_ndarray(), other.as_ndarray()))


class Rotation2d:
    @staticmethod
    def fromRadians(rad: float):
        return Rotation2d(math.cos(rad), math.sin(rad))

    @staticmethod
    def fromDegrees(deg: float):
        return Rotation2d.fromRadians(math.radians(deg))

    rotation = Vector2d(1.0, 0.0)

    def __init__(self, cos=1.0, sin=0.0, *, vec=None):
        if vec is not None:
            self.rotation = vec
        else:
            self.rotation = Vector2d(cos, sin)

    @property
    def tan(self):
        return self.rotation.y / self.rotation.x if self.rotation.x > kEpsilon \
            else math.inf if self.rotation.y >= 0.0 \
            else -math.inf

    @property
    def radians(self):
        return math.atan2(self.rotation.y, self.rotation.x)

    @property
    def degrees(self):
        return math.degrees(self.radians)

    @property
    def theta(self):
        return self.radians

    @property
    def sin(self):
        return self.rotation.y

    @property
    def cos(self):
        return self.rotation.x

    @property
    def normal(self):
        return Rotation2d(self.sin, -self.cos)

    def rotateBy(self, other):
        if type(other) == type(self):
            return Rotation2d(vec=rotateVec2d(Vector2d(self.cos, self.sin), Vector2d(other.cos, other.sin)))
        else:
            return Rotation2d(vec=rotateVec2d(Vector2d(self.cos, self.sin), other))


class Translation2d:
    pos = Vector2d(0, 0)

    @property
    def x(self):
        return self.pos.x

    @property
    def y(self):
        return self.pos.y

    def __init__(self, x=0, y=0, *, vec=None):
        if vec is not None:
            self.pos = vec
        else:
            self.pos = Vector2d(x, y)

    @property
    def norm(self):
        return np.linalg.norm(self.pos.as_ndarray())

    @property
    def normalized(self):
        return Translation2d(vec=self.pos.normalize())

    def translateBy(self, other):
        return Translation2d(self.x + other.x, self.y + other.y)

    def __add__(self, other):
        return self.translateBy(other)

    def __neg__(self):
        return Translation2d(-self.x, -self.y)

    def __sub__(self, other):
        return self + -other

    def __mul__(self, other):
        return Translation2d(self.x * other, self.y * other)

    def __truediv__(self, other):
        return Translation2d(self.x / other, self.y / other)

    def __str__(self):
        return self.pos.__str__()


class RigidTransform2d:
    pos: Translation2d
    rot: Rotation2d

    def __init__(self, pos: Translation2d = Translation2d(), rot: Rotation2d = Rotation2d(), *, copy=None):
        if copy is not None:
            self.pos = copy.pos
            self.rot = copy.rot
        else:
            self.pos = pos
            self.rot = rot

    @property
    def x(self):
        return self.pos.x

    @property
    def y(self):
        return self.pos.y

    @property
    def cos(self):
        return self.rot.cos

    @property
    def sin(self):
        return self.rot.sin

    @property
    def theta(self):
        return self.rot.theta


def rotateVec2d(src: Vector2d, other: Vector2d) -> Vector2d:
    norm = other.normalize()
    rotmtx = np.array(
        [
            [norm.x, -norm.y],
            [norm.y, norm.x]
        ]
    )
    #    return rotmtx * np.asarray([
    #        [src.x],
    #        [src.y]
    #    ])

    return Vector2d.from_ndarray(np.matmul(rotmtx, src.as_ndarray()).reshape(1, 2)[0])


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

if __name__ == '__main__':
    print(
        rotateVec2d(
            src=Vector2d(1.0, 0.0),
            other=Vector2d(math.sqrt(2) / 2, math.sqrt(2) / 2)
        )
    )
    print(Translation2d(1.0, 5.0) / 5.0)
