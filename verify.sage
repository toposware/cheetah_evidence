# Copyright (c) 2022-2023 Toposware, Inc.
#
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

"""
This module provides evidences of the security of the Cheetah elliptic
curve defined over a sextic extension of the prime field GF(2^64 - 2^32 + 1).
"""

from utils import *
from util_hashtocurve import OptimizedSSWU
from constants import *


def verify():
    k = GF(P)
    kx = k['x']
    k6 = k.extension(kx(POLY), "u")

    p_isprime = P.is_prime(proof=True)
    q_isprime = CURVE_PRIME_ORDER.is_prime(proof=True)

    E = EllipticCurve(k6, [k6(CURVE_COEFF_A), k6(CURVE_COEFF_B)])

    # Enforce that the curve has `CURVE_FULL_ORDER` points, and that
    # `CURVE_FULL_ORDER` equals `CURVE_PRIME_ORDER` times `CURVE_COFACTOR`.
    assert(E.count_points() == CURVE_FULL_ORDER)
    assert(CURVE_FULL_ORDER == CURVE_PRIME_ORDER * CURVE_COFACTOR)

    G = E(k6(CURVE_GENERATOR_X), k6(CURVE_GENERATOR_Y))

    # Enforce that basepoint is of order `CURVE_PRIME_ORDER`
    assert(G * CURVE_COFACTOR != E(0, 1, 0))
    assert(G * CURVE_PRIME_ORDER == E(0, 1, 0))

    # Enforce that basepoint has been generated from the "Cheetah"
    # string and the SSWU hashing-to-curve algorithm
    cheetah_sswu = OptimizedSSWU(k6, k6(CURVE_COEFF_A), k6(CURVE_COEFF_B))
    bin = BinaryStrings()
    sswu_bin_encoding = bin.encoding("Cheetah")
    sswu_int = k6(int(str(sswu_bin_encoding), 2))
    g_sswu = cheetah_sswu.map_to_curve(sswu_int)

    # The obtained point is not yet on the prime-order subgroup
    assert(g_sswu * CURVE_PRIME_ORDER != E(0, 1, 0))
    g_sswu = g_sswu * CURVE_COFACTOR
    # We can now enforce equality with the hardcoded basepoint
    assert(G == g_sswu)

    # Compute Pollard-Rho security and embedding degree for the curve and its twist
    e_security = generic_curve_security(
        P ** 6, CURVE_FULL_ORDER, CURVE_PRIME_ORDER, CURVE_PRIME_ORDER_MINUS_ONE_FACTORS)
    t_security = generic_twist_security(
        P ** 6, CURVE_FULL_ORDER, TWIST_PRIME_ORDER, TWIST_PRIME_ORDER_MINUS_ONE_FACTORS)

    is_pollard_rho_secure = e_security[0] > POLLARD_RHO_SECURITY
    twist_is_pollard_rho_secure = t_security[0] > POLLARD_RHO_TWIST_SECURITY
    is_mov_secure = e_security[1].nbits() > EMBEDDING_DEGREE_SECURITY
    twist_is_mov_secure = t_security[1] > EMBEDDING_DEGREE_SECURITY

    # Check complex discriminant
    # p^6 + 1 - t = #E
    t = P ** 6 + 1 - CURVE_FULL_ORDER
    D = t ** 2 - 4*P ** 6
    assert(prod(x ** y for x, y in DISCRIMINANT_FACTORS)
           == -D)
    for (factor, _) in DISCRIMINANT_FACTORS:
        while D % factor**2 == 0:
            D //= factor**2

    if D % 4 != 1:
        D *= 4
    is_discriminant_large = D < -2 ** DISCRIMINANT_SECURITY

    # Check sextic-extension specific attack security of the curve
    is_genus_2_secure = genus_2_cover_security(E)
    is_genus_3_h_secure = genus_3_hyperelliptic_cover_security(
        E, CURVE_FULL_ORDER)
    is_genus_3_nh_secure = genus_3_nonhyperelliptic_cover_security(E)
    is_ghs_secure = ghs_security(k6(CURVE_COEFF_A), k6(CURVE_COEFF_B), k6)

    # Print final results
    display_result(p_isprime, q_isprime, CURVE_PRIME_ORDER.nbits(), is_pollard_rho_secure, is_mov_secure, e_security, twist_is_pollard_rho_secure,
                   twist_is_mov_secure, t_security, is_discriminant_large, D.nbits(), is_genus_2_secure, is_genus_3_h_secure, is_genus_3_nh_secure, is_ghs_secure)


verify()
