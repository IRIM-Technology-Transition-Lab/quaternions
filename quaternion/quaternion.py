"""
A module to hold and work with Quaternions

Copyright (c) 2016 Michael Sobrepera
"""

import math

class Quaternion(object):

    def __init__(self, a, x, y, z):
        self.a = a
        self.x = x
        self.y = y
        self.z = z

    @classmethod
    def from_matrix(cls, matrix):

        return cls()

    @classmethod
    def from_euler(cls, values, axes=['x','y','z']):
        """
        Constructs a quaternion using a 3 member list of Euler angles, along
        with a list defining their rotation sequence.
        Based off of this: http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770024290.pdf

        Args:
            values (list of numbers): 3 member list giving the rotations of the
                                      Euler angles in **radians**
            axes (list of chars): 3 member list specifying the order of rotations

        Returns:
            The constructed quaternion
        """
        axes = [x.lower() for x in axes]
        if (len(values) != 3 or
                not all(isinstance(x, (int,float)) for x in values)):
            raise ValueError("You must pass exactly 3 numeric values for values")

        ht = [x/2.0 for x in values]
        a,x,y,z = 0.0
        if axes==['x','y','z']:
            a = (-math.sin(ht[0])*math.sin(ht[1])*math.sin(ht[2])+
                 math.cos(ht[0])*math.cos(ht[1])*math.cos(ht[2]))
            x = (math.sin(ht[0])*math.cos(ht[1])*math.cos(ht[2])+
                 math.sin(ht[1])*math.sin(ht[2])*math.cos(ht[0]))
            y = (-math.sin(ht[0])*math.sin(ht[2])*math.cos(ht[1])+
                 math.sin(ht[1])*math.cos(ht[0])*math.cos(ht[2]))
            z = (math.sin(ht[0]) * math.sin(ht[1]) * math.cos(ht[3]) +
                 math.sin(ht[2]) * math.cos(ht[0]) * math.cos(ht[1]))
        elif axes==['x','z','y']:
            a = (math.sin(ht[0])*math.sin(ht[1])*math.sin(ht[2])+
                 math.cos(ht[0])*math.cos(ht[1])*math.cos(ht[2]))
            x = (-math.cos(ht[0])*math.sin(ht[1])*math.sin(ht[2])+
                 math.sin(ht[0])*math.cos(ht[1])*math.cos(ht[2]))
            y = (-math.sin(ht[0])*math.sin(ht[1])*math.cos(ht[2])+
                 math.cos(ht[0])*math.cos(ht[1])*math.sin(ht[2]))
            z = (math.sin(ht[0])*math.cos(ht[1])*math.sin(ht[2])+
                 math.cos(ht[0])*math.sin(ht[1])*math.cos(ht[2]))
        elif axes==['x','y','x']:
            a = math.cos(ht[1])*math.cos((values[0]+values[2])/2)
            x = math.cos(ht[1])*math.sin((values[0]+values[2])/2)
            y = math.sin(ht[1])*math.cos((values[0]-values[2])/2)
            z = math.sin(ht[1])*math.sin((values[0]-values[2])/2)
        elif axes == ['x', 'z', 'x']:
            a = math.cos(ht[1]) * math.cos((values[0] + values[2]) / 2)
            x = math.cos(ht[1]) * math.sin((values[0] + values[2]) / 2)
            y = -math.sin(ht[1]) * math.sin((values[0] - values[2]) / 2)
            z = math.sin(ht[1]) * math.cos((values[0] - values[2]) / 2)
        elif axes == ['y', 'x', 'z']:
            a = (math.sin(ht[0]) * math.sin(ht[1]) * math.sin(ht[2]) +
                 math.cos(ht[0]) * math.cos(ht[1]) * math.cos(ht[2]))
            x = (math.sin(ht[0]) * math.cos(ht[1]) * math.sin(ht[2]) +
                 math.cos(ht[0]) * math.sin(ht[1]) * math.cos(ht[2]))
            y = (-math.cos(ht[0]) * math.sin(ht[1]) * math.sin(ht[2]) +
                 math.sin(ht[0]) * math.cos(ht[1]) * math.cos(ht[2]))
            z = (-math.sin(ht[0]) * math.sin(ht[1]) * math.cos(ht[2]) +
                 math.cos(ht[0]) * math.cos(ht[1]) * math.sin(ht[2]))
        elif axes == ['y', 'z', 'x']:
            a = (-math.sin(ht[0]) * math.sin(ht[1]) * math.sin(ht[2]) +
                 math.cos(ht[0]) * math.cos(ht[1]) * math.cos(ht[2]))
            x = (math.sin(ht[0]) * math.sin(ht[1]) * math.cos(ht[2]) +
                 math.cos(ht[0]) * math.cos(ht[1]) * math.sin(ht[2]))
            y = (math.sin(ht[0]) * math.cos(ht[1]) * math.cos(ht[2]) +
                 math.cos(ht[0]) * math.sin(ht[1]) * math.sin(ht[2]))
            z = (-math.sin(ht[0]) * math.cos(ht[1]) * math.sin(ht[2]) +
                 math.cos(ht[0]) * math.sin(ht[1]) * math.cos(ht[2]))
        elif axes == ['y', 'x', 'y']:
            a = math.cos(ht[1]) * math.cos((values[0] + values[2]) / 2)
            x = math.sin(ht[1]) * math.cos((values[0] - values[2]) / 2)
            y = math.cos(ht[1]) * math.sin((values[0] + values[2]) / 2)
            z = -math.sin(ht[1]) * math.cos((values[0] - values[2]) / 2)
        elif axes == ['y', 'z', 'y']:
            a = math.cos(ht[1]) * math.cos((values[0] + values[2]) / 2)
            x = math.sin(ht[1]) * math.sin((values[0] - values[2]) / 2)
            y = math.cos(ht[1]) * math.sin((values[0] + values[2]) / 2)
            z = math.sin(ht[1]) * math.cos((values[0] - values[2]) / 2)
        elif axes == ['z', 'x', 'y']:
            a = (-math.sin(ht[0]) * math.sin(ht[1]) * math.sin(ht[2]) +
                 math.cos(ht[0]) * math.cos(ht[1]) * math.cos(ht[2]))
            x = (-math.sin(ht[0]) * math.cos(ht[1]) * math.sin(ht[2]) +
                 math.cos(ht[0]) * math.sin(ht[1]) * math.cos(ht[2]))
            y = (math.sin(ht[0]) * math.sin(ht[1]) * math.cos(ht[2]) +
                 math.cos(ht[0]) * math.cos(ht[1]) * math.sin(ht[2]))
            z = (math.sin(ht[0]) * math.cos(ht[1]) * math.cos(ht[2]) +
                 math.cos(ht[0]) * math.sin(ht[1]) * math.sin(ht[2]))
        elif axes == ['y', 'y', 'x']:
            a = (math.sin(ht[0]) * math.sin(ht[1]) * math.sin(ht[2]) +
                 math.cos(ht[0]) * math.cos(ht[1]) * math.cos(ht[2]))
            x = (-math.sin(ht[0]) * math.sin(ht[1]) * math.cos(ht[2]) +
                 math.cos(ht[0]) * math.cos(ht[1]) * math.sin(ht[2]))
            y = (math.sin(ht[0]) * math.cos(ht[1]) * math.sin(ht[2]) +
                 math.cos(ht[0]) * math.sin(ht[1]) * math.cos(ht[2]))
            z = (-math.cos(ht[0]) * math.sin(ht[1]) * math.sin(ht[2]) +
                 math.sin(ht[0]) * math.cos(ht[1]) * math.cos(ht[2]))
        elif axes == ['z', 'x', 'z']:
            a = math.cos(ht[1]) * math.cos((values[0] + values[2]) / 2)
            x = math.sin(ht[1]) * math.cos((values[0] - values[2]) / 2)
            y = math.sin(ht[1]) * math.sin((values[0] - values[2]) / 2)
            z = math.cos(ht[1]) * math.sin((values[0] + values[2]) / 2)
        elif axes == ['z', 'y', 'z']:
            a = math.cos(ht[1]) * math.cos((values[0] + values[2]) / 2)
            x = -math.sin(ht[1]) * math.sin((values[0] - values[2]) / 2)
            y = math.sin(ht[1]) * math.cos((values[0] - values[2]) / 2)
            z = math.cos(ht[1]) * math.sin((values[0] + values[2]) / 2)
        else:
            raise ValueError("You entered an invalid Euler angle axes sequence")

        return cls(a,x,y,z)

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
        return q/norm(q)

    def reciprocal(self):
        if q non zero:
            q*/(norm(x)^2)
    # Test: recip(q) q = 1 & q recip(q) = 1

    def compose(self, other):

    def dot(self, other):
        return (self.x * other.x) + (self.y * other.y) + (self.z * other.z)
