import unittest
import os

from TransducerCompiler import compileTransducerToCCode

from sonLib.misc import sonTraceRootPath
 
class TestCase(unittest.TestCase):
    
    def setUp(self):
        self.rootPath = sonTraceRootPath() + "/src/ortheus/models"
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