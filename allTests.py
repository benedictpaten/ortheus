#!/usr/bin/env python

#Copyright (C) 2008-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest

import OrtheusTests
import TransducerCompilerTest
import TransducerComposerTest

from sonLib.bioio import parseSuiteTestOptions

def allSuites():
    OrtheusSuite = unittest.makeSuite(OrtheusTests.TestCase, 'test')
    TransducerCompilerSuite = unittest.makeSuite(TransducerCompilerTest.TestCase, 'test')
    TransducerComposerSuite = unittest.makeSuite(TransducerComposerTest.TestCase, 'test')
    
    allTests = unittest.TestSuite((OrtheusSuite, TransducerCompilerSuite, TransducerComposerSuite))
    return allTests
        
def main():
    parseSuiteTestOptions()
    
    suite = allSuites()
    runner = unittest.TextTestRunner()
    runner.run(suite)
        
if __name__ == '__main__':
    main()