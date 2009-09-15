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