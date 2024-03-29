# Copyright (c) 2022-2023 Toposware, Inc.
#
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

"""
Utility module for hashing to curve with the optimized version
of the Simplified Shallue-van de Woestijne-Ulas method
for elliptic curves in short Weierstrass with both A and B coefficients
non-zero.

Adapted from https://github.com/cfrg/draft-irtf-cfrg-hash-to-curve/
"""

from sage.all import *


# Helper functions

def CMOV(x, y, b):
    """
    Returns x if b=False; otherwise returns y
    """
    return int(not bool(b))*x + int(bool(b))*y


def find_z_sswu(F, A, B):
    R = F['xx']
    (xx,) = R._first_ngens(1)  # Polynomial ring over F
    g = xx ** 3 + F(A) * xx + F(B)  # y^2 = g(x) = x^3 + A * x + B
    ctr = F.gen()
    while True:
        for Z_cand in (F(ctr), F(-ctr)):
            # Criterion 1: Z is non-square in F.
            if is_square(Z_cand):
                continue
            # Criterion 2: Z != -1 in F.
            if Z_cand == F(-1):
                continue
            # Criterion 3: g(x) - Z is irreducible over F.
            if not (g - Z_cand).is_irreducible():
                continue
            # Criterion 4: g(B / (Z * A)) is square in F.
            if is_square(g(B / (Z_cand * A))):
                return Z_cand
        ctr += 1


def sgn0(x):
    """
    Returns 1 if x is 'negative' (little-endian sense), else 0.
    """
    ZZR = PolynomialRing(ZZ, name='XX')

    degree = x.parent().degree()
    if degree == 1:
        # not a field extension
        xi_values = (ZZ(x),)
    else:
        # field extension
        # extract vector repr of field element (faster than x._vector_())
        xi_values = ZZR(x)
    sign = 0
    zero = 1
    # compute the sign in constant time
    for i in range(0, degree):
        zz_xi = xi_values[i]
        # sign of this digit
        sign_i = zz_xi % 2
        zero_i = zz_xi == 0
        # update sign and zero
        sign = sign | (zero & sign_i)
        zero = zero & zero_i
    return sign


# cache for per-p values
sqrt_cache = {}


def square_root_random_sign(x):
    a = square_root(x)
    if a is not None and randint(0, 1) == 1:
        return -a
    return a


def square_root(x):
    """
    Returns a square root defined through fixed formulas.
    (non-constant-time)
    """
    F = x.parent()
    p = F.order()

    if p % 16 == 1:
        return tonelli_shanks_ct(x)

    if p % 4 == 3:
        if sqrt_cache.get(p) is None:
            sqrt_cache[p] = (F(1),)
        z = x ** ((p + 1) // 4)

    if p % 8 == 5:
        if sqrt_cache.get(p) is None:
            sqrt_cache[p] = (F(1), F(-1).sqrt())
        z = x ** ((p + 3) // 8)

    elif p % 16 == 9:
        if sqrt_cache.get(p) is None:
            sqrt_m1 = F(-1).sqrt()
            sqrt_sqrt_m1 = sqrt_m1.sqrt()
            sqrt_cache[p] = (F(1), sqrt_m1, sqrt_sqrt_m1,
                             sqrt_sqrt_m1 * sqrt_m1)
        z = x ** ((p + 7) // 16)

    for mul in sqrt_cache[p]:
        sqrt_cand = z * mul
        if sqrt_cand ** 2 == x:
            return sqrt_cand

    return None


# constant-time Tonelli-Shanks
# Adapted from https://github.com/zkcrypto/jubjub/blob/master/src/fq.rs by Michael Scott.
# See also Cohen, "A Course in Computational # Algebraic Number Theory," Algorithm 1.5.1.
def tonelli_shanks_ct(x):
    F = x.parent()
    p = F.order()
    if sqrt_cache.get(p) is None:
        ts_precompute(p, F)
    (c1, c3, c5) = sqrt_cache[p]

    z = x ** c3
    t = z * z
    t = t * x
    z = z * x
    b = t
    c = c5
    for i in range(c1, 1, -1):
        for j in range(1, i - 1):
            b = b * b
        e = b == 1
        zt = z * c
        z = CMOV(zt, z, e)
        c = c * c
        tt = t * c
        t = CMOV(tt, t, e)
        b = t

    if z ** 2 == x:
        return z
    assert not x.is_square()
    return None


# cache pre-computable values -- no need for CT here
def ts_precompute(p, F):
    c2 = p - 1
    c1 = 0
    while c2 % 2 == 0:
        c2 //= 2
        c1 += 1
    assert c2 == (p - 1) // (2 ** c1)
    c4 = F.gen()
    while c4.is_square():
        c4 += 1
    assert p == c2 * 2**c1 + 1
    c3 = (c2 - 1) // 2
    c5 = c4 ** c2
    sqrt_cache[p] = (c1, c3, c5)


def _get_lo(q):
    o = q - 1
    l = 0
    while o % 2 == 0:
        o = o // 2
        l = l + 1
    assert o * 2 ** l == q - 1
    return (l, o)


# This routine identifies a suitable Z in the absence of one provided by the
# caller of sqrt_checked or sqrt_ratio_straightline. When those functions are
# called from elsewhere in the codebase, callers generally provide a Z value;
# this routine's return values are only used for testing.
#
# The values returned by this routine should not be used in hash-to-curve
# implementations. Use the Z value generated by the appropriate function in
# z_selection.sage instead.
def find_Z(F):
    (l, o) = _get_lo(F.order())
    ctr = F.gen()
    ll0 = 2 ** l
    ll1 = 2 ** (l-1)
    while True:
        Z_cand = F(ctr)
        Z_cand_o = Z_cand ** o
        s2l0 = Z_cand_o ** ll0
        s2l1 = Z_cand_o ** ll1
        if not Z_cand.is_square():
            assert s2l1 != F(1)
            assert s2l0 == F(1)
            return Z_cand
        ctr += 1


_s_vals = {}


def _get_Z_val(F):
    global _s_vals
    if F in _s_vals:
        return _s_vals[F]
    pe = find_Z(F)
    _s_vals[F] = pe
    return pe


_consts = {}


def _get_consts(F, Z):
    global _consts
    if (F, Z) in _consts:
        return _consts[(F, Z)]
    q = F.order()
    (l, o) = _get_lo(q)
    # c1, the largest integer such that 2^c1 divides q - 1.
    c1 = l
    c2 = (q - 1) / (2 ** c1)        # Integer arithmetic
    assert c2 == o
    c3 = (c2 - 1) / 2            # Integer arithmetic
    c4 = 2 ** c1 - 1                # Integer arithmetic
    c5 = 2 ** (c1 - 1)              # Integer arithmetic
    c6 = Z ** c2
    c7 = Z ** ((c2 + 1) / 2)
    ret = (c1, c3, c4, c5, c6, c7)
    _consts[(F, Z)] = ret
    return ret


def sqrt_checked(F, x, Z=None):
    x = F(x)
    isQR = True
    order = F.order()
    m = 0
    r = order - 1
    while r % 2 == 0:
        r = r / 2
        m += 1
    assert 2 ** m * r == order-1, "bad initialization"
    if Z is None:
        Z = _get_Z_val(F)
    z = x ** ((r-1)/2)
    t = z * z * x  # x^r
    z = z * x  # x^((r+1)/2)
    c = Z ** r
    inital_tweak_z = Z ** ((r+1)/2)
    if t ** (2 ** (m-1)) != 1:
        isQR = False
        assert not is_square(x), "incorrect determination of squareness"
        z = z*inital_tweak_z
        t = t*c

    for i in range(m, 1, -1):
        if t ** (2 ** (i-2)) != 1:
            z = z * c
            t = t * c * c
        c = c * c
    if isQR:
        assert z*z == x, "incorrect square root: %s squared is not %s" % (z, x)
    if not isQR:
        assert z*z == x * \
            Z, "incorrect tweaked square root: %s squared is not %s" % (z, x*Z)
    return (isQR, z)


def sqrt_ratio_straightline(F, u, v, Z=None):
    u = F(u)
    v = F(v)
    if Z is None:
        Z = _get_Z_val(F)
    (c1, c3, c4, c5, c6, c7) = _get_consts(F, Z)

    tv1 = c6
    tv2 = v ** c4
    tv3 = tv2 ** 2
    tv3 = tv3 * v
    tv5 = u * tv3
    tv5 = tv5 ** c3
    tv5 = tv5 * tv2
    tv2 = tv5 * v
    tv3 = tv5 * u
    tv4 = tv3 * tv2
    tv5 = tv4 ** c5
    isQR = tv5 == 1
    tv2 = tv3 * c7
    tv5 = tv4 * tv1
    tv3 = CMOV(tv2, tv3, isQR)
    tv4 = CMOV(tv5, tv4, isQR)
    for i in range(c1, 1, -1):
        tv5 = i - 2
        tv5 = 2 ** tv5
        tv5 = tv4 ** tv5
        e1 = tv5 == 1
        tv2 = tv3 * tv1
        tv1 = tv1 * tv1
        tv5 = tv4 * tv1
        tv3 = CMOV(tv2, tv3, e1)
        tv4 = CMOV(tv5, tv4, e1)

    assert (isQR, tv3) == sqrt_checked(F, u/v, Z), "incorrect sqrt_ratio"
    return (isQR, tv3)


# Helper class for hashing to curve

class GenericSSWU(object):
    def __init__(self, F, A, B):
        self.name = "SSWU"
        self.F = F
        self.A = F(A)
        self.B = F(B)
        if self.A == 0:
            raise ValueError("S-SWU requires A != 0")
        if self.B == 0:
            raise ValueError("S-SWU requires B != 0")
        self.Z = find_z_sswu(F, F(A), F(B))
        self.E = EllipticCurve(F, [F(A), F(B)])

        # constants for straight-line impl
        self.c1 = -F(B) / F(A)
        self.c2 = -F(1) / self.Z

        # values at which the map is undefined
        # i.e., when Z^2 * u^4 + Z * u^2 = 0
        # which is at u = 0 and when Z * u^2 = -1
        self.undefs = [F(0)]
        if self.c2.is_square():
            ex = self.c2.sqrt()
            self.undefs += [ex, -ex]

        self.sqrt = square_root_random_sign

    def not_straight_line(self, u):
        inv0 = self.inv0
        is_square = self.is_square
        u = self.F(u)
        A = self.A
        B = self.B
        Z = self.Z

        tv1 = inv0(Z**2 * u**4 + Z * u**2)
        x1 = (-B / A) * (1 + tv1)
        if tv1 == 0:
            x1 = B / (Z * A)
        gx1 = x1**3 + A * x1 + B
        x2 = Z * u**2 * x1
        gx2 = x2**3 + A * x2 + B
        if is_square(gx1):
            x = x1
            y = square_root(gx1)
        else:
            x = x2
            y = square_root(gx2)
        if sgn0(u) != sgn0(y):
            y = -y
        return (x, y)

    def straight_line(self, u):
        inv0 = self.inv0
        is_square = self.is_square
        u = self.F(u)
        A = self.A
        B = self.B
        Z = self.Z
        c1 = self.c1
        c2 = self.c2

        tv1 = Z * u**2
        tv2 = tv1**2
        x1 = tv1 + tv2
        x1 = inv0(x1)
        e1 = x1 == 0
        x1 = x1 + 1
        x1 = CMOV(x1, c2, e1)    # If (tv1 + tv2) == 0, set x1 = -1 / Z
        x1 = x1 * c1      # x1 = (-B / A) * (1 + (1 / (Z^2 * u^4 + Z * u^2)))
        gx1 = x1**2
        gx1 = gx1 + A
        gx1 = gx1 * x1
        gx1 = gx1 + B             # gx1 = g(x1) = x1^3 + A * x1 + B
        x2 = tv1 * x1            # x2 = Z * u^2 * x1
        tv2 = tv1 * tv2
        gx2 = gx1 * tv2           # gx2 = (Z * u^2)^3 * gx1
        e2 = is_square(gx1)
        x = CMOV(x2, x1, e2)    # If is_square(gx1), x = x1, else x = x2
        y2 = CMOV(gx2, gx1, e2)  # If is_square(gx1), y2 = gx1, else y2 = gx2
        y = square_root(y2)
        e3 = sgn0(u) == sgn0(y)  # Fix sign of y
        y = CMOV(-y, y, e3)
        return (x, y)

    def map_to_curve(self, u):
        (x1, y1) = self.straight_line(u)
        (x2, y2) = self.not_straight_line(u)
        assert (x1, y1) == (
            x2, y2), "straight-line / non-straight-line mismatch"
        return self.E((x1, y1))

    def is_square(self, x):
        return self.F(x).is_square()

    def inv0(self, x):
        if self.F(x) == self.F(0):
            return self.F(0)
        return self.F(1) / self.F(x)


class OptimizedSSWU:
    def __init__(self, F, A, B):
        self.name = "Optimized SSWU"
        self.F = F
        self.A = F(A)
        self.B = F(B)

        if self.A == 0:
            raise ValueError("S-SWU requires A != 0")
        if self.B == 0:
            raise ValueError("S-SWU requires B != 0")
        self.Z = find_z_sswu(F, F(A), F(B))
        self.E = EllipticCurve(F, [F(A), F(B)])

        self.ref_map = GenericSSWU(F, A, B)
        self.undefs = self.ref_map.undefs

    def map_to_curve(self, u):
        (x1, y1) = self.straight_line(u)
        (x2, y2) = self.not_straight_line(u)
        assert (x1, y1) == (
            x2, y2), "straight-line / non-straight-line mismatch"
        return self.E((x1, y1))

    def is_square(self, x):
        return self.F(x).is_square()

    def inv0(self, x):
        if self.F(x) == self.F(0):
            return self.F(0)
        return self.F(1) / self.F(x)

    def not_straight_line(self, u):
        inv0 = self.inv0
        is_square = self.is_square
        u = self.F(u)
        A = self.A
        B = self.B
        Z = self.Z

        tv1 = inv0(Z ** 2 * u ** 4 + Z * u ** 2)
        x1 = (-B / A) * (1 + tv1)
        if tv1 == 0:
            x1 = B / (Z * A)
        gx1 = x1 ** 3 + A * x1 + B
        x2 = Z * u ** 2 * x1
        gx2 = x2 ** 3 + A * x2 + B
        if is_square(gx1):
            x = x1
            y = square_root(gx1)
        else:
            x = x2
            y = square_root(gx2)
        if sgn0(u) != sgn0(y):
            y = -y

        (xp, yp, zp) = self.ref_map.map_to_curve(u)
        xp = xp / zp
        yp = yp / zp
        assert xp == x
        assert yp == y

        return (x, y)

    def sqrt_ratio(self, u, v):
        x = self.F(u) / self.F(v)
        r1 = sqrt_checked(self.F, x, self.Z)
        r2 = sqrt_ratio_straightline(self.F, u, v, self.Z)
        assert r1 == r2
        return r2

    def straight_line(self, u):
        A = self.A
        B = self.B
        Z = self.Z
        u = self.F(u)
        sqrt_ratio = self.sqrt_ratio

        tv1 = u ** 2
        tv1 = Z * tv1
        tv2 = tv1 ** 2
        tv2 = tv2 + tv1
        tv3 = tv2 + 1
        tv3 = B * tv3
        tv4 = CMOV(Z, -tv2, tv2 != 0)
        tv4 = A * tv4
        tv2 = tv3 ** 2
        tv6 = tv4 ** 2
        tv5 = A * tv6
        tv2 = tv2 + tv5
        tv2 = tv2 * tv3
        tv6 = tv6 * tv4
        tv5 = B * tv6
        tv2 = tv2 + tv5
        x = tv1 * tv3
        (is_gx1_square, y1) = sqrt_ratio(tv2, tv6)
        y = tv1 * u
        y = y * y1
        x = CMOV(x, tv3, is_gx1_square)
        y = CMOV(y, y1, is_gx1_square)
        e1 = sgn0(u) == sgn0(y)
        y = CMOV(-y, y, e1)
        x = x / tv4

        return (x, y)
