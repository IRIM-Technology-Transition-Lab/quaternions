from unittest import TestCase

import math

from quaternion import Quaternion
from random import Random

class TestQuaternion(TestCase):
    def setUp(self):
        random = Random(100)
        N = 50*4
        int_vals = [random.randint(-1000,1000) for _ in range(N/2)]
        float_vals = [random.uniform(-1000, 1000) for _ in range(N / 2)]
        val_list = int_vals + float_vals
        random.shuffle(val_list)
        # val_list = \
        #     [0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 3, 4.5, 1.5, 2.5, 0, 9.7] + val_list
        self.p = [Quaternion(val_list[4*n], val_list[4*n+1], val_list[4*n+2],
                             val_list[4*n+3], ) for n in range(N/8)]
        self.p = [p.unit() for p in self.p]
        self.q = [Quaternion(val_list[4*n], val_list[4*n+1], val_list[4*n+2],
                             val_list[4*n+3], ) for n in range(N/8, N/4)]
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