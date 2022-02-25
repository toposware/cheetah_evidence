# Copyright (c) 2022 Toposware, Inc.
#
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

"""
This module provides all the necessary constants to build and verify the
security of the Cheetah elliptic curve, formed on a sextic extension
of the prime field GF(2^64 - 2^32 + 1).
"""

from sage.all import *

# Parameters defining the extension
P = Integer(18446744069414584321)

# Polynomials defining the extension tower
POLY = "x**6 - 7"


# Parameters defining the Cheetah elliptic curve and its prime subgroup.
CURVE_FULL_ORDER = Integer(
    39402006141350512473373051550956862238057064441100873236175653175294391087266801759100731275759467840677123787392470)
CURVE_PRIME_ORDER = Integer(
    55610362957290864006699123731285679659474893560816383126640993521607086746831)
CURVE_COFACTOR = Integer(708537115134665106932687062569690615370)

CURVE_COEFF_A = 1
CURVE_COEFF_B = "u + 395"

# Prime-order subgroup generator obtained from the Simplified Shallue-van de Woestijne-Ulas method
CURVE_GENERATOR_X = "12938930721685970739*u**5 + 375185138577093320*u**4 + 4830863958577994148*u**3 + 10526511002404673680*u**2 + 8599518745794843693*u + 2754611494552410273"
CURVE_GENERATOR_Y = "9990732138772505951*u**5 + 13187678623570541764*u**4 + 10708493419890101954*u**3 + 14375303400746062753*u**2 + 2774812795997841935*u + 15384029202802550068"

TWIST_PRIME_ORDER = Integer(
    223601976666913148532398089120886101575454553929079671735989704519146384604414660366213)


# Helper constants for speeding-up embedding degree calculations
CURVE_PRIME_ORDER_MINUS_ONE_FACTORS = [(2, 1),
                                       (3, 2),
                                       (5, 1),
                                       (617892921747676266741101374792063107327498817342404256962677705795634297187,
                                        1)]

TWIST_PRIME_ORDER_MINUS_ONE_FACTORS = [(2, 2),
                                       (7, 2),
                                       (1883330971, 1),
                                       (3474824543011033669, 1),
                                       (174325146328985902076700361979837689750196287435948514503, 1)]

DISCRIMINANT_FACTORS = [(2, 2),
                        (3, 1),
                        (5, 1),
                        (41, 1),
                        (139, 1),
                        (1877, 1),
                        (135277, 1),
                        (14221247, 1),
                        (62876612242092976553, 1),
                        (258086030176912481681567, 1),
                        (6680106624555714761968517908536894284995258092121529, 1)]
