"""
Tests for quaternion class
"""

# The MIT License (MIT)
#
# Copyright (c) 2016 GTRC.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


from unittest import TestCase

import math

from quaternion import Quaternion
from random import Random

class TestQuaternion(TestCase):
    def setUp(self):
        random = Random(100)
        N = 50*4
        int_vals = [random.randint(-1000,1000) for _ in range(int(N/2))]
        float_vals = [random.uniform(-1000, 1000) for _ in range(int(N / 2))]
        val_list = int_vals + float_vals
        random.shuffle(val_list)
        # val_list = \
        #     [0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 3, 4.5, 1.5, 2.5, 0, 9.7] + val_list
        self.p = [Quaternion(val_list[4*n], val_list[4*n+1], val_list[4*n+2],
                             val_list[4*n+3], ) for n in range(int(N/8))]
        self.p = [p.unit() for p in self.p]
        self.q = [Quaternion(val_list[4*n], val_list[4*n+1], val_list[4*n+2],
                             val_list[4*n+3], ) for n in range(int(N/8), int(N/4))]
        self.q = [q.unit() for q in self.q]

        self.all = self.p + self.q

    def test_equality(self):
        for quat in self.all:
            self.assertTrue(quat == quat)

    def test_from_Quaternion(self):
        for quat in self.all:
            self.assertEqual(Quaternion.from_quaternion(quat),
                         quat)

    def test_conjugate_a(self):
        """:math: q=(q*)*"""
        for q in self.q:
            self.assertEqual(q, q.conjugate().conjugate())

    def test_conjugate_b(self):
        """:math: (pq)* = q*p*"""
        for q in self.q:
            for p in self.p:
                self.assertTrue((p*q).conjugate().almost_equal(q.conjugate()*p.conjugate()),
                                "Failed conjugation test b:\n\tq: " + str(q) +
                                "\n\tp: " + str(p) +
                                "\n\nLeft:\t" + str((p*q).conjugate()) +
                                "\nRight:\t" + str((q.conjugate())*(p.conjugate())))

    def test_norm_a(self):
        """norm(conjugate(q))=norm(q)"""
        for q in self.all:
            self.assertAlmostEqual(q.conjugate().norm(), q.norm())

    def test_norm_b(self):
        """norm(pq) = norm(p)norm(q)"""
        for q in self.q:
            for p in self.p:
                self.assertAlmostEqual((p*q).norm(), p.norm()*q.norm())

    def test_inverse_a(self):
        """inverse(q) q = Quaternion(1,0,0,0)"""
        for q in self.all:
            self.assertTrue((q.inverse()*q).almost_equal(Quaternion(1, 0, 0, 0)))

    def test_inverse_b(self):
        """q inverse(q) = Quaternion(1,0,0,0)"""
        for q in self.all:
            self.assertTrue((q*q.inverse()).almost_equal(Quaternion(1, 0, 0, 0)))

    def test_inverse_c(self):
        """q inverse(q) = inverse(q) q"""
        for q in self.all:
            self.assertTrue((q * q.inverse()).almost_equal(q.inverse()*q))

    def test_inverse_d(self):
        """inverse(inverse(q)) = q"""
        for q in self.all:
            self.assertTrue(q.inverse().inverse().almost_equal(q))

    def test_inverse_e(self):
        """inverse(pq) = inverse(q)inverse(p)"""
        for q in self.q:
            for p in self.p:
                self.assertTrue((p*q).inverse().almost_equal(q.inverse()*p.inverse()))
