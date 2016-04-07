import math

class Quaternion(object):

    def __init__(self, a, x, y, z):
        self.a = a
        self.x = x
        self.y = y
        self.z = z

    @classmethod
    def from_matrix(cls, matrix):

    @classmethod
    def from_euler(cls, values, axes):

    @classmethod
    def from_quaternion(cls, quaternion):

    @classmethod
    def from_translation(cls, translation):

    @classmethod
    def from_axis_angle(cls, axis, angle):
        sin_result = math.sin(angle/2)
        return cls(math.cos(angle/2), axis[0]*sin_result, axis[1]*sin_result, axis[2]*sin_result)

    def __add__(self, other):
        return Quaternion(self.a + other.a, self.x + other.x, self.y + other.y, self.z + other.z)

    def __mul__(self, other):
        if isinstance(other, Quaternion):
            return Quaternion((self.a * other.a) - (self.x * other.x) - (self.y * other.y) - (self.z * other.z),
                              (self.a * other.x) + (self.x * other.a) + (self.y * other.z) - (self.z * other.y),
                              (self.a * other.y) - (self.y * other.z) + (self.y * other.a) + (self.z * other.y),
                              (self.a * other.z) + (self.x * other.y) - (self.y * other.x) + (self.z * other.a))
        elif isinstance(other,(int, float)):
            return Quaternion(self.a * other, self.x * other, self.y * other, self.z * other)
        else:
            return NotImplemented

    def conjugate(self):
        return Quaternion(self.a, -self.x, -self.y, -self.z)
    # Test: conjugate of conjugate equals original
    # Test: (pq)* = q*p*

    def norm(self):
        return math.sqrt(math.pow(self.a, 2) + math.pow(self.x, 2) + math.pow(self.y, 2) + math.pow(self.z, 2))
    # Test: norm(x) = sqrt(x x*) = sqrt(x* x) > 0 & Real
    # Test: norm(a q) = abs(real) norm(q) given a is scalar
    # Test: norm(p q) = norm(p) norm(q)
    # Test: norm(p + a p1 + q + aq1) - (p+q)) = a norm(p1 + q1) given a is scalar

    def __sub__(self, other):
        return Quaternion(self.a - other.a, self.x - other.x, self.y - other.y, self.z - other.z)

    def distance(self, other):
        return (self - other).norm()

    def unit(self):
        q/norm(q)

    def reciprocal(self):
        if q non zero:
            q*/(norm(x)^2)
    # Test: recip(q) q = 1 & q recip(q) = 1

    def compose(self, other):

    def dot(self, other):
        return (self.x * other.x) + (self.y * other.y) + (self.z * other.z)
