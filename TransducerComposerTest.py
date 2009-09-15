import unittest
import os

from TransducerComposer import parseInputTransducerFile
from TransducerComposer import mainScript

from sonLib.misc import sonTraceRootPath

class TestCase(unittest.TestCase):
    
    def setUp(self):
        self.rootPath = sonTraceRootPath() + "/src/ortheus/models"
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        unittest.TestCase.tearDown(self)
        
    def testParseInputTransducerFile(self):
        inputFile = self.rootPath + "/affineModel.sxpr"
        branchTransducerX, branchTransducerZ, rootTransducer = parseInputTransducerFile(inputFile)
        branchTransducerX.report()
        branchTransducerZ.report()
        rootTransducer.report()
        
    def testMainScript(self):
        inputFile1 = self.rootPath + "/affineModel.sxpr"
        outputFile = self.rootPath + "/temp.otra"
        outputDotFile = self.rootPath + "/temp.dot"
        mainScript(inputFile1, outputFile, outputDotFile)
        
        inputFile1 = self.rootPath + "/geometric_parent.sxpr"
        mainScript(inputFile1, outputFile, outputDotFile)
        
        inputFile1 = self.rootPath + "/doubleAffineModel.sxpr"
        mainScript(inputFile1, outputFile, outputDotFile)
        
        #Remove outputFile
        os.remove(outputFile)
        os.remove(outputDotFile)
     
if __name__ == '__main__':
    unittest.main()