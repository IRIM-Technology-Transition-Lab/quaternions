"""
A module to hold and work with Quaternions

Copyright (c) 2016 Michael Sobrepera
"""

import math

class Quaternion(object):

    def __init__(self, w, x, y, z):
        self.w = w
        self.x = x
        self.y = y
        self.z = z

    @classmethod
    def from_matrix(cls, matrix):
        """Generates a Quaternion from a rotation matrix

        Args:
            matrix: A rotation matrix, a list of 3, 3 member lists of numbers

        Returns:
            The constructed quaternion
        """
        if not isinstance(matrix, list) and all(isinstance(sub, list) and all(
            isinstance(val, (int, float)) for val in sub) for sub in matrix):
            raise ValueError("input must be a list of 3 3 member lists of numbers")
        w = math.sqrt(1+matrix[0][0]+matrix[1][1]+matrix[2][2])/2
        x = (matrix[2][1]-matrix[1][2])/(4*w)
        y = (matrix[0][2]-matrix[2][0])/(4*w)
        z = (matrix[1][0]-matrix[0][1])/(4*w)

        return cls(w,x,y,z)

    @classmethod
    def from_euler(cls, values, axes=['x','y','z']):
        """
        Constructs w quaternion using w 3 member list of Euler angles, along
        with w list defining their rotation sequence.
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
        if axes==['x','y','z']:
            w = (-math.sin(ht[0])*math.sin(ht[1])*math.sin(ht[2])+
                 math.cos(ht[0])*math.cos(ht[1])*math.cos(ht[2]))
            x = (math.sin(ht[0])*math.cos(ht[1])*math.cos(ht[2])+
                 math.sin(ht[1])*math.sin(ht[2])*math.cos(ht[0]))
            y = (-math.sin(ht[0])*math.sin(ht[2])*math.cos(ht[1])+
                 math.sin(ht[1])*math.cos(ht[0])*math.cos(ht[2]))
            z = (math.sin(ht[0]) * math.sin(ht[1]) * math.cos(ht[3]) +
                 math.sin(ht[2]) * math.cos(ht[0]) * math.cos(ht[1]))
        elif axes==['x','z','y']:
            w = (math.sin(ht[0])*math.sin(ht[1])*math.sin(ht[2])+
                 math.cos(ht[0])*math.cos(ht[1])*math.cos(ht[2]))
            x = (-math.cos(ht[0])*math.sin(ht[1])*math.sin(ht[2])+
                 math.sin(ht[0])*math.cos(ht[1])*math.cos(ht[2]))
            y = (-math.sin(ht[0])*math.sin(ht[1])*math.cos(ht[2])+
                 math.cos(ht[0])*math.cos(ht[1])*math.sin(ht[2]))
            z = (math.sin(ht[0])*math.cos(ht[1])*math.sin(ht[2])+
                 math.cos(ht[0])*math.sin(ht[1])*math.cos(ht[2]))
        elif axes==['x','y','x']:
            w = math.cos(ht[1])*math.cos((values[0]+values[2])/2)
            x = math.cos(ht[1])*math.sin((values[0]+values[2])/2)
            y = math.sin(ht[1])*math.cos((values[0]-values[2])/2)
            z = math.sin(ht[1])*math.sin((values[0]-values[2])/2)
        elif axes == ['x', 'z', 'x']:
            w = math.cos(ht[1]) * math.cos((values[0] + values[2]) / 2)
            x = math.cos(ht[1]) * math.sin((values[0] + values[2]) / 2)
            y = -math.sin(ht[1]) * math.sin((values[0] - values[2]) / 2)
            z = math.sin(ht[1]) * math.cos((values[0] - values[2]) / 2)
        elif axes == ['y', 'x', 'z']:
            w = (math.sin(ht[0]) * math.sin(ht[1]) * math.sin(ht[2]) +
                 math.cos(ht[0]) * math.cos(ht[1]) * math.cos(ht[2]))
            x = (math.sin(ht[0]) * math.cos(ht[1]) * math.sin(ht[2]) +
                 math.cos(ht[0]) * math.sin(ht[1]) * math.cos(ht[2]))
            y = (-math.cos(ht[0]) * math.sin(ht[1]) * math.sin(ht[2]) +
                 math.sin(ht[0]) * math.cos(ht[1]) * math.cos(ht[2]))
            z = (-math.sin(ht[0]) * math.sin(ht[1]) * math.cos(ht[2]) +
                 math.cos(ht[0]) * math.cos(ht[1]) * math.sin(ht[2]))
        elif axes == ['y', 'z', 'x']:
            w = (-math.sin(ht[0]) * math.sin(ht[1]) * math.sin(ht[2]) +
                 math.cos(ht[0]) * math.cos(ht[1]) * math.cos(ht[2]))
            x = (math.sin(ht[0]) * math.sin(ht[1]) * math.cos(ht[2]) +
                 math.cos(ht[0]) * math.cos(ht[1]) * math.sin(ht[2]))
            y = (math.sin(ht[0]) * math.cos(ht[1]) * math.cos(ht[2]) +
                 math.cos(ht[0]) * math.sin(ht[1]) * math.sin(ht[2]))
            z = (-math.sin(ht[0]) * math.cos(ht[1]) * math.sin(ht[2]) +
                 math.cos(ht[0]) * math.sin(ht[1]) * math.cos(ht[2]))
        elif axes == ['y', 'x', 'y']:
            w = math.cos(ht[1]) * math.cos((values[0] + values[2]) / 2)
            x = math.sin(ht[1]) * math.cos((values[0] - values[2]) / 2)
            y = math.cos(ht[1]) * math.sin((values[0] + values[2]) / 2)
            z = -math.sin(ht[1]) * math.cos((values[0] - values[2]) / 2)
        elif axes == ['y', 'z', 'y']:
            w = math.cos(ht[1]) * math.cos((values[0] + values[2]) / 2)
            x = math.sin(ht[1]) * math.sin((values[0] - values[2]) / 2)
            y = math.cos(ht[1]) * math.sin((values[0] + values[2]) / 2)
            z = math.sin(ht[1]) * math.cos((values[0] - values[2]) / 2)
        elif axes == ['z', 'x', 'y']:
            w = (-math.sin(ht[0]) * math.sin(ht[1]) * math.sin(ht[2]) +
                 math.cos(ht[0]) * math.cos(ht[1]) * math.cos(ht[2]))
            x = (-math.sin(ht[0]) * math.cos(ht[1]) * math.sin(ht[2]) +
                 math.cos(ht[0]) * math.sin(ht[1]) * math.cos(ht[2]))
            y = (math.sin(ht[0]) * math.sin(ht[1]) * math.cos(ht[2]) +
                 math.cos(ht[0]) * math.cos(ht[1]) * math.sin(ht[2]))
            z = (math.sin(ht[0]) * math.cos(ht[1]) * math.cos(ht[2]) +
                 math.cos(ht[0]) * math.sin(ht[1]) * math.sin(ht[2]))
        elif axes == ['y', 'y', 'x']:
            w = (math.sin(ht[0]) * math.sin(ht[1]) * math.sin(ht[2]) +
                 math.cos(ht[0]) * math.cos(ht[1]) * math.cos(ht[2]))
            x = (-math.sin(ht[0]) * math.sin(ht[1]) * math.cos(ht[2]) +
                 math.cos(ht[0]) * math.cos(ht[1]) * math.sin(ht[2]))
            y = (math.sin(ht[0]) * math.cos(ht[1]) * math.sin(ht[2]) +
                 math.cos(ht[0]) * math.sin(ht[1]) * math.cos(ht[2]))
            z = (-math.cos(ht[0]) * math.sin(ht[1]) * math.sin(ht[2]) +
                 math.sin(ht[0]) * math.cos(ht[1]) * math.cos(ht[2]))
        elif axes == ['z', 'x', 'z']:
            w = math.cos(ht[1]) * math.cos((values[0] + values[2]) / 2)
            x = math.sin(ht[1]) * math.cos((values[0] - values[2]) / 2)
            y = math.sin(ht[1]) * math.sin((values[0] - values[2]) / 2)
            z = math.cos(ht[1]) * math.sin((values[0] + values[2]) / 2)
        elif axes == ['z', 'y', 'z']:
            w = math.cos(ht[1]) * math.cos((values[0] + values[2]) / 2)
            x = -math.sin(ht[1]) * math.sin((values[0] - values[2]) / 2)
            y = math.sin(ht[1]) * math.cos((values[0] - values[2]) / 2)
            z = math.cos(ht[1]) * math.sin((values[0] + values[2]) / 2)
        else:
            raise ValueError("You entered an invalid Euler angle axes sequence")

        return cls(w,x,y,z)

    @classmethod
    def from_quaternion(cls, quaternion):
        """
        Constructs a quaternion from another quaternion.

        Args:
            quaternion (Quaternion): The quaternion to copy

        Returns:
            The newly constructed quaternion
        """
        if not isinstance(quaternion, Quaternion):
            return ValueError("input must be a quaternion")
        return cls(quaternion.w, quaternion.x, quaternion.y, quaternion.z)

    @classmethod
    def from_translation(cls, translation):

    @classmethod
    def from_axis_angle(cls, axis, angle):
        sin_result = math.sin(angle/2)
        return cls(math.cos(angle/2), axis[0]*sin_result, axis[1]*sin_result, axis[2]*sin_result)

    def __add__(self, other):
        if not isinstance(other, Quaternion):
            raise ValueError("Must pass in another Quaternion")
        return Quaternion(self.w + other.w, self.x + other.x, self.y + other.y, self.z + other.z)

    def __mul__(self, other):
        if isinstance(other, Quaternion):
            return Quaternion((self.w * other.w) - (self.x * other.x) - (self.y * other.y) - (self.z * other.z),
                              (self.w * other.x) + (self.x * other.w) + (self.y * other.z) - (self.z * other.y),
                              (self.w * other.y) - (self.y * other.z) + (self.y * other.w) + (self.z * other.y),
                              (self.w * other.z) + (self.x * other.y) - (self.y * other.x) + (self.z * other.w))
        elif isinstance(other,(int, float)):
            return Quaternion(self.w * other, self.x * other, self.y * other, self.z * other)
        else:
            return NotImplemented

    def conjugate(self):
        return Quaternion(self.w, -self.x, -self.y, -self.z)
    # Test: conjugate of conjugate equals original
    # Test: (pq)* = q*p*

    def norm(self):
        return math.sqrt(math.pow(self.w, 2) + math.pow(self.x, 2) + math.pow(self.y, 2) + math.pow(self.z, 2))
    # Test: norm(x) = sqrt(x x*) = sqrt(x* x) > 0 & Real
    # Test: norm(w q) = abs(real) norm(q) given w is scalar
    # Test: norm(p q) = norm(p) norm(q)
    # Test: norm(p + w p1 + q + aq1) - (p+q)) = w norm(p1 + q1) given w is scalar

    def __sub__(self, other):
        if not isinstance(other, Quaternion):
            raise ValueError("Must pass in another Quaternion")
        return Quaternion(self.w - other.w, self.x - other.x, self.y - other.y, self.z - other.z)

    def distance(self, other):
        if not isinstance(other, Quaternion):
            raise ValueError("Must pass in another Quaternion")
        return (self - other).norm()

    def unit(self):
        return self/norm(self)

    def reciprocal(self):
        if q non zero:
            q*/(norm(x)^2)
    # Test: recip(q) q = 1 & q recip(q) = 1

    def compose(self, other):

    def dot(self, other):
        return (self.x * other.x) + (self.y * other.y) + (self.z * other.z)
