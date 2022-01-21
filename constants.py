"""
This module provides all the necessary constants to build and verify the
security of the Cheetah elliptic curve, formed on a sextic extension
of the prime field GF(2^64 - 2^32 + 1).
"""

from sage.all import *

# Parameters defining the extension
P = Integer(18446744069414584321)
DEGREE = 6
POWER_PRIMITIVE_ELEMENT = 311

# Polynomials defining the extension tower
POLY = "x^6 - 7"


# Parameters defining the Cheetah elliptic curve and its prime subgroup.
CURVE_FULL_ORDER = Integer(
    39402006141350512473373051550956862238057064441100873236169646240204104116866551216944808951746172933195042092384306)
CURVE_PRIME_ORDER = Integer(
    10230194559610033405039867617070259612247645045591847851798073552054039295467)
CURVE_COFACTOR = Integer(3851540252901356796008861125671805917718)

CURVE_COEFF_A = 1
CURVE_COEFF_B = "8751648879042009540*u^5 + 9152663345541292493*u^4 + 11461655964319546493*u^3 + 18139339604299112321*u^2 + 11484739320139982777*u + 13239264519837388909"

CURVE_GENERATOR_X = "14525701791208974225*u^5 + 9317930114507436616*u^4 + 15366750451009235558*u^3 + 1487465174323083563*u^2 + 7015637262057020447*u + 7245987309612903154"
CURVE_GENERATOR_Y = "15793083952783722763*u^5 + 5153214154709431636*u^4 + 17741001138250422543*u^3 + 5112079552571193492*u^2 + 14203600126873428360*u + 5846570486458364631"

TWIST_PRIME_ORDER = Integer(
    2777277120284427568023445991898494413128327016410668251003909727886171705563353)


# Helper constants for speeding-up embedding degree calculations
CURVE_PRIME_ORDER_MINUS_ONE_FACTORS = [(2, 1),
                                       (3, 2),
                                       (7, 1),
                                       (13, 1),
                                       (31, 1),
                                       (1777, 1),
                                       (113375933054658078597970753655786304020500605907940181244640080727161, 1)]

TWIST_PRIME_ORDER_MINUS_ONE_FACTORS = [(2, 3),
                                       (3, 2),
                                       (2428087638393691, 1),
                                       (2659884960708926281, 1),
                                       (5972546065456619376387003508506761614559921, 1)]
