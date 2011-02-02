#!/usr/bin/env python

#Copyright (C) 2008-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt

import unittest
import os

from TransducerCompiler import compileTransducerToCCode

from ortheus.common import ortheusRootPath
 
class TestCase(unittest.TestCase):
    
    def setUp(self):
        self.rootPath = os.path.join(ortheusRootPath(), "models")
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        unittest.TestCase.tearDown(self)
        
    def testCompileTransducer(self):
        inputFile = self.rootPath + "/affineModel.otra"
        outputCFile = self.rootPath + "/affineModel.c"
        outputHFile = self.rootPath + "/affineModel.h"
        compileTransducerToCCode(inputFile, outputCFile, outputHFile)
        #Remove outputFile
        os.remove(outputCFile)
        os.remove(outputHFile)
        
if __name__ == '__main__':
    unittest.main()