"""
This module provides evidences of the security of the Cheetah elliptic curve defined over
a sextic extension of the prime field Fp with p = 2^64 - 2^32 + 1.
"""

from utils import *
from constants import *


def verify():
    k = GF(P)
    kx = k['x']
    k6 = k.extension(kx(POLY), "u")

    p_isprime = P.is_prime(proof=True)
    q_isprime = CURVE_PRIME_ORDER.is_prime(proof=True)

    E = EllipticCurve(k6, [k6(CURVE_COEFF_A), k6(CURVE_COEFF_B)])
    assert(k6.primitive_element() **
           POWER_PRIMITIVE_ELEMENT == k6(CURVE_COEFF_B))

    # Enforce that the curve has `CURVE_FULL_ORDER` points, and that
    # `CURVE_FULL_ORDER` equals `CURVE_PRIME_ORDER` times `CURVE_COFACTOR`.
    assert(E.count_points() == CURVE_FULL_ORDER)
    assert(CURVE_FULL_ORDER == CURVE_PRIME_ORDER * CURVE_COFACTOR)

    G = E(k6(CURVE_GENERATOR_X), k6(CURVE_GENERATOR_Y))

    # Enforce that basepoint is of order `CURVE_PRIME_ORDER`
    assert(G.order() == CURVE_PRIME_ORDER)
    assert(G * CURVE_COFACTOR != E(0, 1, 0))
    assert(G * CURVE_PRIME_ORDER == E(0, 1, 0))

    # Enforce that hardcoded helper constants are valid
    assert(TWIST_PRIME_ORDER.is_prime(proof=True))
    assert((2*(P ^ 6 + 1) - CURVE_FULL_ORDER) % TWIST_PRIME_ORDER == 0)
    assert(prod(x ^ y for x, y in CURVE_PRIME_ORDER_MINUS_ONE_FACTORS)
           == CURVE_PRIME_ORDER - 1)
    assert(prod(x ^ y for x, y in TWIST_PRIME_ORDER_MINUS_ONE_FACTORS)
           == TWIST_PRIME_ORDER - 1)

    # Compute Pollard-Rho security and embedding degree for the curve and its twist
    e_security = curve_security(
        P ^ 6, CURVE_FULL_ORDER, CURVE_PRIME_ORDER, CURVE_PRIME_ORDER_MINUS_ONE_FACTORS)
    t_security = twist_security(
        P ^ 6, CURVE_FULL_ORDER, TWIST_PRIME_ORDER, TWIST_PRIME_ORDER_MINUS_ONE_FACTORS)

    is_pollard_rho_secure = e_security[0] > RHO_SECURITY
    twist_is_pollard_rho_secure = t_security[0] > TWIST_SECURITY
    is_mov_secure = e_security[1].nbits() > EMBEDDING_DEGREE_SECURITY
    twist_is_mov_secure = t_security[1] > EMBEDDING_DEGREE_SECURITY

    # Check sextic-extension specific attack security of the curve
    is_genus_2_secure = genus_2_cover_security(E)
    is_genus_3_h_secure = genus_3_hyperelliptic_cover_security(
        E, CURVE_FULL_ORDER)
    is_genus_3_nh_secure = genus_3_nonhyperelliptic_cover_security(E)
    is_ghs_secure = ghs_security(k6(CURVE_COEFF_A), k6(CURVE_COEFF_B), k6)

    # Print final results
    display_result(p_isprime, q_isprime, CURVE_PRIME_ORDER.nbits(), is_pollard_rho_secure, is_mov_secure, e_security, twist_is_pollard_rho_secure,
                   twist_is_mov_secure, t_security, is_genus_2_secure, is_genus_3_h_secure, is_genus_3_nh_secure, is_ghs_secure)


verify()
