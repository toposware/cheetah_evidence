"""
This module provides all the necessary constants to build and verify the
security of the Cheetah elliptic curve, formed on a sextic extension
of the prime field GF(2^64 - 2^32 + 1).
"""

from sage.all import *

# Parameters defining the extension
P = Integer(18446744069414584321)
DEGREE = 6
POWER_PRIMITIVE_ELEMENT = 294

# Polynomials defining the extension tower
POLY = "x^6 - x - 1"


# Parameters defining the Cheetah elliptic curve and its prime subgroup.
CURVE_FULL_ORDER = Integer(
    39402006141350512473373051550956862238057064441100873236167583409106056278293958754762526768910210209442281490706790)
CURVE_PRIME_ORDER = Integer(
    47901524144134838524951555146609336724916023524245614293878204099229447205501)
CURVE_COFACTOR = Integer(822562681362508909289329747633262861790)

CURVE_COEFF_A = 1
CURVE_COEFF_B = "6035437078588083401*u^5 + 2857736027621020694*u^4 + 14652430667295839026*u^3 + 9953599263706139154*u^2 + 4878031766911544164*u + 1789756228530856276"

CURVE_GENERATOR_X = "7794598347260010540*u^5 + 10618622314870859645*u^4 + 14820135150250105254*u^3 + 4679529240359230485*u^2 + 16711805998942165686*u + 8337925652927234608"
CURVE_GENERATOR_Y = "3536746701953062266*u^5 + 1688150812217466895*u^4 + 13679753399522663948*u^3 + 14150323329710684588*u^2 + 3361356312544043514*u + 1202746335776013372"

TWIST_PRIME_ORDER = Integer(
    856565350898924184203761990238192657349066618284801592090738568505074030580094487355638833440455248956461653523849)


# Helper constants for speeding-up embedding degree calculations
CURVE_PRIME_ORDER_MINUS_ONE_FACTORS = [(2, 2),
                                       (3, 5),
                                       (5, 3),
                                       (23, 1),
                                       (127, 1),
                                       (88547, 1),
                                       (817183, 1),
                                       (1865298419040846456677757911581088148744131952002113469237, 1)]

TWIST_PRIME_ORDER_MINUS_ONE_FACTORS = [(2, 3),
                                       (3, 1),
                                       (613, 1),
                                       (23560183264403, 1),
                                       (2471212671089594991510266465766648346427860392607611500021183414432851518179366995174720853085893,
                                        1)]
