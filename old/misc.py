#!/usr/bin/env python

#Copyright (C) 2008-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt

import sys
import os
import re
import math

#########################################################
#########################################################
#########################################################
#misc functions
#########################################################
#########################################################
#########################################################

def linOriginRegression(points):
    """
    computes a linear regression starting at zero
    """
    j = sum([ i[0] for i in points ])
    k = sum([ i[1] for i in points ])
    if j != 0:
        return k/j, j, k
    return 1, j, k

def close(i, j, tolerance):
    """
    check two float values are within a bound of one another
    """
    return i <= j + tolerance and i >= j - tolerance

#########################################################
#########################################################
#########################################################
#continuous math functions
#########################################################
#########################################################
#########################################################

LOG_ZERO_PROB = -1e30000
LOG_ONE_PROB = 0.0
ZERO_PROB = 0.0

def logAdd(x, y):
    if x < y:
        if x <= LOG_ZERO_PROB:
            return y
        return math.log(math.exp(x - y) + 1) + y
    if y <= LOG_ZERO_PROB:
        return x
    return math.log(math.exp(y - x) + 1) + x

#########################################################
#########################################################
#########################################################
#bio functions
#########################################################
#########################################################
#########################################################

dNAMap_reverseComp_IUPAC = { 'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'n':'n' }

dNAMap_reverseComp_Int = { 0:3, 1:2, 2:1, 3:0, 4:4 }

def reverseComplement(seq, rCM=dNAMap_reverseComp_Int):
    seq.reverse()
    i = [ rCM[j] for j in seq ]
    seq.reverse()
    return i


def main():
    pass

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()
