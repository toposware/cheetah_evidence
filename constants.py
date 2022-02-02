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
POLY = "x**6 - 7"


# Parameters defining the Cheetah elliptic curve and its prime subgroup.
CURVE_FULL_ORDER = Integer(
    39402006141350512473373051550956862238057064441100873236169646240204104116866551216944808951746172933195042092384306)
CURVE_PRIME_ORDER = Integer(
    10230194559610033405039867617070259612247645045591847851798073552054039295467)
CURVE_COFACTOR = Integer(3851540252901356796008861125671805917718)

CURVE_COEFF_A = 1
CURVE_COEFF_B = "8751648879042009540*u**5 + 9152663345541292493*u**4 + 11461655964319546493*u**3 + 18139339604299112321*u**2 + 11484739320139982777*u + 13239264519837388909"

# Prime-order subgroup generator obtained from the Simplified Shallue-van de Woestijne-Ulas method
CURVE_GENERATOR_X = "3424925235903636224*u**5 + 4220885415633744070*u**4 + 8114383727112424035*u**3 + 2279564461935938533*u**2 + 2626899613631150971*u + 8695031540959941377"
CURVE_GENERATOR_Y = "14965683120853365983*u**5 + 14304619975303237162*u**4 + 744748944745194181*u**3 + 14146767844428843935*u**2 + 8141780031807903371*u + 9566031633564717113"

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

DISCRIMINANT_FACTORS = [(2, 2),
                        (3, 1),
                        (19, 1),
                        (109, 1),
                        (229, 1),
                        (719, 1),
                        (479761, 1),
                        (703127, 1),
                        (31586170873, 1),
                        (89657950031647, 1),
                        (11566669382877326859111792446843, 1),
                        (3457410493340060929559242158224744420789, 1)]
