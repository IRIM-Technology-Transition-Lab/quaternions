"""
A module to hold and work with Quaternions
Repo at: https://github.com/mjsobrep/quaternions
"""

from __future__ import division
import math


class Quaternion(object):
    """ A class to hold and work with quaternions

    Attributes:
        w:  The real component of the quaternion
        x:  The i component of the quaternion
        y:  The y component of the quaternion
        z:  The z component of the quaternion
        norm_error: The maximum deviation from magnitude 1 for which a
                    quaternion will be considered normalized
    """

    def __init__(self, w, x, y, z, norm_error=0.00000001):
        """ Constructs a quaternion from individual components

        Args:
            w (int,float): The real component of the quaternion
            x (int,float): The i component of the quaternion
            y (int,float): The j component of the quaternion
            z (int,float): The k component of the quaternion
            norm_error (float): The maximum deviation from magnitude 1 for
                                which a quaternion will be considered
                                normalized
        """
        self.w = w
        self.x = x
        self.y = y
        self.z = z
        self.norm_error = norm_error

    @classmethod
    def from_matrix(cls, matrix):
        """Generates a Quaternion from a rotation matrix

        Args:
            matrix: A rotation matrix, a list of 3, 3 member lists of numbers

        Returns:
            The constructed quaternion
        """
        if not isinstance(matrix, list) and all(isinstance(sub, list) and all(
                isinstance(val, (int, float)) for val in sub)
                for sub in matrix):
            raise ValueError(
                "input must be a list of 3 3 member lists of numbers")
        w = math.sqrt(1+matrix[0][0]+matrix[1][1]+matrix[2][2])/2
        x = (matrix[2][1]-matrix[1][2])/(4*w)
        y = (matrix[0][2]-matrix[2][0])/(4*w)
        z = (matrix[1][0]-matrix[0][1])/(4*w)

        return cls(w, x, y, z)

    @classmethod
    def from_euler(cls, values, axes=['x', 'y', 'z']):
        """
        Constructs w quaternion using w 3 member list of Euler angles, along
        with w list defining their rotation sequence.
        Based off of this:
        http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770024290.pdf

        Args:
            values (list of numbers): 3 member list giving the rotations of the
                                      Euler angles in **radians**
            axes (list of chars): 3 member list specifying the order of
                                  rotations

        Returns:
            The constructed quaternion
        """
        axes = [x.lower() for x in axes]
        if (len(values) != 3 or
                not all(isinstance(x, (int, float)) for x in values)):
            raise ValueError(
                "You must pass exactly 3 numeric values for values")

        ht = [x/2.0 for x in values]
        if axes == ['x', 'y', 'z']:
            w = (-math.sin(ht[0])*math.sin(ht[1])*math.sin(ht[2]) +
                 math.cos(ht[0])*math.cos(ht[1])*math.cos(ht[2]))
            x = (math.sin(ht[0])*math.cos(ht[1])*math.cos(ht[2]) +
                 math.sin(ht[1])*math.sin(ht[2])*math.cos(ht[0]))
            y = (-math.sin(ht[0])*math.sin(ht[2])*math.cos(ht[1]) +
                 math.sin(ht[1])*math.cos(ht[0])*math.cos(ht[2]))
            z = (math.sin(ht[0]) * math.sin(ht[1]) * math.cos(ht[2]) +
                 math.sin(ht[2]) * math.cos(ht[0]) * math.cos(ht[1]))
        elif axes == ['x', 'z', 'y']:
            w = (math.sin(ht[0])*math.sin(ht[1])*math.sin(ht[2]) +
                 math.cos(ht[0])*math.cos(ht[1])*math.cos(ht[2]))
            x = (-math.cos(ht[0])*math.sin(ht[1])*math.sin(ht[2]) +
                 math.sin(ht[0])*math.cos(ht[1])*math.cos(ht[2]))
            y = (-math.sin(ht[0])*math.sin(ht[1])*math.cos(ht[2]) +
                 math.cos(ht[0])*math.cos(ht[1])*math.sin(ht[2]))
            z = (math.sin(ht[0])*math.cos(ht[1])*math.sin(ht[2]) +
                 math.cos(ht[0])*math.sin(ht[1])*math.cos(ht[2]))
        elif axes == ['x', 'y', 'x']:
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
            raise ValueError(
                "You entered an invalid Euler angle axes sequence")

        return cls(w, x, y, z)

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
        """
        Generates a quaternion from a translation. Meant to be used along with
        the dual operator to generate dual quaternions

        Args:
            translation: A list of three numbers representing the (x,y,z)
            translation

        Returns:
            The quaternion created from the translation
        """
        if (len(translation) != 3 or
                not all(isinstance(x, (int, float)) for x in translation)):
            raise ValueError(
                "You must pass exactly 3 numeric values for the translation")
        return cls(0, translation[0]/2, translation[1]/2, translation[2]/2)

    @classmethod
    def from_axis_angle(cls, axis, angle):
        """
        Constructs a quaternion from axis-angle measurements.

        Args:
            axis: A list of three numbers representing the axis of rotation.
            angle: A single number representing the rotation about the axis in
                   **radians**

        Returns:
            The constructed quaternion
        """
        if (len(axis) != 3 or
                not all(isinstance(x, (int, float)) for x in axis)):
            raise ValueError(
                "You must pass exactly 3 numeric values for the axis")
        if not isinstance(angle, (int, float)):
            raise ValueError("The angle value must be a single numeric value")
        sin_result = math.sin(angle/2)
        return cls(math.cos(angle/2), axis[0]*sin_result, axis[1]*sin_result,
                   axis[2]*sin_result)

    @classmethod
    def from_rotation_vector(cls, vect):
        """
        Constructs a quaternion from a rotation vector.

        Args:
            vect: A rotation vector with angle as the magnitude and vector for
                  the vector.

        Returns:
            The constructed quaternion
        """
        angle = math.sqrt(sum(element**2 for element in vect))
        if angle == 0:
            axis = [1, 0, 0]
        else:
            axis = [element/angle for element in vect]
        return cls.from_axis_angle(axis, angle)

    def __add__(self, other):
        """
        Add together two quaternions.

        Args:
            other (Quaternion): The quaternion to add.

        Returns:
            The result of the addition.
        """
        if not isinstance(other, Quaternion):
            raise ValueError("Must pass in another Quaternion")
        return Quaternion(self.w + other.w, self.x + other.x, self.y + other.y,
                          self.z + other.z)

    def __mul__(self, other):
        """
        Multiply two quaternions.

        Args:
            other (Quaternion): The quaternion to multiply by. The caller will
                                be on the left side and the callee argument
                                will be on the right side

        Returns:
            The result of the multiplication.
        """
        if isinstance(other, Quaternion):
            return Quaternion((self.w * other.w) - (self.x * other.x) -
                              (self.y * other.y) - (self.z * other.z),
                              (self.w * other.x) + (self.x * other.w) +
                              (self.y * other.z) - (self.z * other.y),
                              (self.w * other.y) - (self.x * other.z) +
                              (self.y * other.w) + (self.z * other.x),
                              (self.w * other.z) + (self.x * other.y) -
                              (self.y * other.x) + (self.z * other.w))
        elif isinstance(other, (int, float)):
            return Quaternion(self.w * other, self.x * other, self.y * other,
                              self.z * other)
        else:
            return NotImplemented

    def conjugate(self):
        """Return the conjugate of the quaternion.

        Returns:
            The conjugate of the quaternion
        """
        return Quaternion(self.w, -self.x, -self.y, -self.z)

    def norm(self):
        """Return the norm of the quaternion

        Returns:
            The norm of the quaternion
        """
        return math.sqrt(math.pow(self.w, 2) + math.pow(self.x, 2) +
                         math.pow(self.y, 2) + math.pow(self.z, 2))

    def unit(self):
        """
        Return the unit Quaternion

        Returns:
            Unit quaternion
        """
        return self/self.norm()

    def __sub__(self, other):
        """Subtract the argument from the caller

        Args:
            other (Quaternion): The quaternion to subtract by

        Returns:
            The result of subtraction
        """
        if not isinstance(other, Quaternion):
            raise ValueError("Must pass in another Quaternion")
        return Quaternion(self.w - other.w, self.x - other.x,
                          self.y - other.y, self.z - other.z)

    def __div__(self, other):
        """
        Divide a quaternion by a scalar value.

        Args:
            other (number): The number to divide by

        Returns:
            The result of division
        """
        if isinstance(other, (int, float)):
            return Quaternion(self.w / other, self.x / other,
                              self.y / other, self.z / other)
        else:
            return NotImplemented

    __truediv__ = __div__

    def distance(self, other):
        """
        Return the rotational distance between two quaternions.

        Args:
            other (Quaternion): The other

        Returns:
            The angular distance between two quaternions
        """
        if not isinstance(other, Quaternion):
            raise ValueError("Must pass in another Quaternion")
        return (self - other).norm()

    def inverse(self):
        """Return the inverse of the quaternion. Specifically, this is the
        multiplicative inverse, such that multiplying the inverse of q by q
        yields the multiplicative identity.

        Returns:
            The multiplicative inverse of the quaternion
        """
        if abs(self.norm()-1) > self.norm_error:
            raise ValueError("The quaternion must be normalized")
        if self:
            return self.conjugate()/self.norm()
        else:
            raise ZeroDivisionError("The quaternion has no non-zero values")

    def dot(self, other):
        """
        Return the quaternion dot product

        Args:
            other (Quaternion): Quaternion with which to calculate the dot
                                product

        Returns:
            The dot product of the two quaternions
        """
        if not isinstance(other, Quaternion):
            raise ValueError("Must pass in another Quaternion")
        return ((self.w * other.w) + (self.x * other.x) + (self.y * other.y) +
                (self.z * other.z))

    def __nonzero__(self):
        """Determines whether the quaternion is non-zero

        Returns:
            True if the quaternion is non-zero, false otherwise
        """
        return any([self.w, self.x, self.y, self.z])

    def __eq__(self, other):
        """Test if two quaternions are equal

        Args:
            other (Quaternion): The quaternion with which to compare equality

        Returns:
            true if the quaternions are equal. False otherwise
        """
        if not isinstance(other, Quaternion):
            raise ValueError("Must pass in another Quaternion")
        return (self.w == other.w and self.x == other.x and
                self.y == other.y and self.z == other.z)

    def __str__(self):
        """Return a string representation of the quaternion in form:
        <w, x, y, z>

        Returns:
            A string representation of the quaternion
        """
        return ("<" + str(self.w) + ", " + str(self.x) + ", " + str(self.y) +
                ", " + str(self.z) + ">")

    def almost_equal(self, other, delta=.00000001):
        """Determines whether a quaternion is approximately equal to another
        using a naive comparison of the 4 values w, x, y, z

        Args:
            other (Quaternion): The quaternion with which to test equality
            delta (float): The error range in which to define two quaternions
                           equal

        Returns:

        """
        result = self-other
        return (abs(result.w) < delta and abs(result.x) < delta and
                abs(result.y) < delta and abs(result.z) < delta)

    def get_rotation_matrix(self):
        """Return the rotation matrix which this quaternion is equivalent to

        Returns:
            The rotation matrix which this quaternion is equivalent to as a
            list of three lists of three elements each
        """
        rot_matx = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        rot_matx[0][0] = (math.pow(self.w, 2) + math.pow(self.x, 2) -
                          math.pow(self.y, 2) - math.pow(self.z, 2))
        rot_matx[0][1] = (2*self.x*self.y - 2*self.w*self.z)
        rot_matx[0][2] = (2*self.x*self.z + 2*self.w*self.y)
        rot_matx[1][0] = (2*self.x*self.y + 2*self.w*self.z)
        rot_matx[1][1] = (math.pow(self.w, 2) - math.pow(self.x, 2) +
                          math.pow(self.y, 2) - math.pow(self.z, 2))
        rot_matx[1][2] = (2*self.y*self.z - 2*self.w*self.x)
        rot_matx[2][0] = (2*self.x*self.z - 2*self.w*self.y)
        rot_matx[2][1] = (2*self.y*self.z + 2*self.w*self.x)
        rot_matx[2][2] = (math.pow(self.w, 2) - math.pow(self.x, 2) -
                          math.pow(self.y, 2) + math.pow(self.z, 2))
        return rot_matx

    def __neg__(self):
        """ Negate the quaternion

        Returns:
            The negated form of the quaternion
        """
        return Quaternion(-self.w, -self.x, -self.y, -self.z)

    def get_euler(self):
        """Return an euler angle representation of the quaternion. Taken from
        [wikipedia](https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles)

        Returns:
            The list of euler angles x,y,z
        """
        w = self.w
        x = self.x
        y = self.y
        z = self.z

        ysqr = y * y

        t0 = +2.0 * (w * x + y * z)
        t1 = +1.0 - 2.0 * (x * x + ysqr)
        X = math.atan2(t0, t1)

        t2 = +2.0 * (w * y - z * x)
        t2 = +1.0 if t2 > +1.0 else t2
        t2 = -1.0 if t2 < -1.0 else t2
        Y = math.asin(t2)

        t3 = +2.0 * (w * z + x * y)
        t4 = +1.0 - 2.0 * (ysqr + z * z)
        Z = math.atan2(t3, t4)

        return [X, Y, Z]

    def get_rotation_vector(self):
        angle = math.acos(self.w)*2
        if angle == 0:
            return [0, 0, 0]
        x = angle * self.x/math.sin(angle/2)
        y = angle * self.y/math.sin(angle/2)
        z = angle * self.z/math.sin(angle/2)
        return [x, y, z]

    def get_xyz_vector(self):
        return [self.x, self.y, self.z]

    @staticmethod
    def average(quats, init, threshold=0.01):
        qt_bar = init
        dist = 5
        while dist > threshold:
            error = [0, 0, 0]
            for element in quats:
                addition = (element*(qt_bar.inverse())).get_rotation_vector()
                error = [error[idx]+addition[idx] for idx in range(3)]
                # TODO: take to axis angle for averaging then back to
                # quaternion
            error = [element/len(quats) for element in error]
            error = Quaternion.from_rotation_vector(error)
            qt_bar_new = error*qt_bar
            dist = qt_bar.distance(qt_bar_new)
            qt_bar = qt_bar_new
        return qt_bar

    @staticmethod
    def vector_average(quats):
        to_return = [0]*4
        for quat in quats:
            to_return[0] += quat.w
            to_return[1] += quat.x
            to_return[2] += quat.y
            to_return[3] += quat.z
        n = len(quats)
        return [to_return[0]/n, to_return[1]/n, to_return[2]/n, to_return[3]/n]
