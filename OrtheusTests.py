"""Tests ortheus_core and the old Ortheus python scripts.
"""

import os
import random
import sys

import unittest

from sonLib.bioio import printBinaryTree
from sonLib.bioio import fastaAlignmentRead
from sonLib.bioio import fastaWrite
from sonLib.bioio import getTempFile

from sonLib.tree import BinaryTree

from sonLib.bioio import system
from sonLib.bioio import getRandomSequence
from sonLib.bioio import mutateSequence

from sonLib.bioio import TestStatus
from sonLib.bioio import parseSuiteTestOptions

class TestCase(unittest.TestCase):
    
    def setUp(self):
        self.testNo = TestStatus.getTestSetup()
        self.tempFiles = []
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        for tempFile in self.tempFiles:
            os.remove(tempFile)
        unittest.TestCase.tearDown(self)
        
        
    def testENm001(self):
        if TestStatus.getTestStatus() == TestStatus.TEST_VERY_LONG:
            encodePath = TestStatus.getPathToDataSets() + "/MAY-2005/ENm001"
            outputPath = TestStatus.getPathToDataSets() + "/ortheus/encodeTest"
            #treeString = '(((((((((((((human:0.006969,chimp:0.009727):0.025291,((baboon:0.008968):0.011019):0.024581):0.023649):0.066673):0.018405,((rat:0.081244,mouse:0.072818):0.238435):0.021892):0.02326,(((cow:0.164728,(cat:0.109852,dog:0.107805):0.049576):0.004663):0.010883):0.033242):0.028346):0.016015):0.226853):0.063898):0.126639):0.119814):0.16696);'
            treeString = '((((human:0.006969,chimp:0.009727):0.025291,baboon:0.044568):0.108727,(rat:0.081244,mouse:0.072818):0.260327):0.02326,(cow:0.164728,(cat:0.109852,dog:0.107805):0.049576):0.048788):0.749525;'
            seqFiles = [ "human.ENm001.fa", "chimp.ENm001.fa", "baboon.ENm001.fa", "rat.ENm001.fa", "mouse.ENm001.fa", "cow.ENm001.fa", "cat.ENm001.fa", "dog.ENm001.fa" ]
            seqFiles = [ encodePath + "/" + i for i in seqFiles ]
            outputFile = outputPath + "/outputENm001.mfa"
            command = "Ortheus.py -e %s -d '%s' -f %s -j -a -b" % \
            (" ".join(seqFiles), treeString, outputFile)
            print "running command", command
            system(command)
        
    def testSimulation(self):
        if TestStatus.getTestStatus() == TestStatus.TEST_LONG:
            blanchettePath = TestStatus.getPathToDataSets() + "/blanchettesSimulation/00.job"
            outputPath = TestStatus.getPathToDataSets() + "/ortheus/blanchettesSimulationTest"
            treeString = '(((((((((((((human:0.006969,chimp:0.009727):0.025291,((baboon:0.008968):0.011019):0.024581):0.023649):0.066673):0.018405,((rat:0.081244,mouse:0.072818):0.238435):0.021892):0.02326,(((cow:0.164728,(cat:0.109852,dog:0.107805):0.049576):0.004663):0.010883):0.033242):0.028346):0.016015):0.226853):0.063898):0.126639):0.119814):0.16696);'
            seqFiles = [ "HUMAN", "CHIMP", "BABOON", "RAT", "MOUSE", "COW", "CAT", "DOG" ]
            seqFiles = [ blanchettePath + "/" + i for i in seqFiles ]
            outputFile = outputPath + "/outputJob1.mfa"
            command = "Ortheus.py -e %s -d '%s' -f %s -j -a -b" % \
            (" ".join(seqFiles), treeString, outputFile)
            print "running command", command
            system(command)
            
    def testAndyYatesFirstExample(self):
        if TestStatus.getTestStatus() == TestStatus.TEST_LONG:
            filePath = TestStatus.getPathToDataSets() + "/ortheus/andyYatesExample1"
            seqs = "seq1.fa seq2.fa seq3.fa seq4.fa seq5.fa seq6.fa seq7.fa seq8.fa seq9.fa seq10.fa seq11.fa \
            seq12.fa seq13.fa seq14.fa seq15.fa seq16.fa seq17.fa seq18.fa seq19.fa seq20.fa seq21.fa seq22.fa seq23.fa seq24.fa seq25.fa seq26.fa \
            seq27.fa seq28.fa seq29.fa seq30.fa seq31.fa seq32.fa seq33.fa seq34.fa seq35.fa seq36.fa"
            seqs = " ".join([ "%s/%s" % (filePath, i) for i in seqs.split() ])
            command = 'Ortheus.py -l "#-j 0 -e" -e %s -z \
            "(((1012:0.0112,1051:0.0119):0.0026,(1055:0.0015,1052:0.0018):0.0370):0.0022,1054:0.0108,1053:0.0116);" \
            -A 1054 1051 1054 1054 1053 1012 1054 1054 1053 1054 1051 1054 1051 1051 1053 1051 1051 1012 1051 1054 1012 1054 1053 1051 1053 \
            1054 1054 1051 1012 1012 1054 1053 1053 1012 1054 1051 -f %s/output.16163.mfa -g %s/output.16163.tree-a -k "# -A" -m "java -Xmx1800m -Xms1800m" -a -b' % \
            (seqs, filePath, filePath)
            print "running command", command
            system(command)
        
    def testRandom(self):
        """Makes random sequences and tests that Ortheus can align them and produce a valid output.
        """
        outputFile = getTempFile()
        self.tempFiles.append(outputFile)
        
        MAX_SEQS = 20
        
        for i in xrange(MAX_SEQS):
            self.tempFiles.append(getTempFile())
        
        for test in xrange(0, self.testNo):
            print "test no : %i " % test
            #seqNo
            binaryTree = randomTree()
            middleSeq = getRandomSequence(250)[1]
            seqs = []
            getTreeSeqs(binaryTree, middleSeq, seqs)
           
            if len(seqs) <= MAX_SEQS and len(seqs) > 2:
                seqFiles = []
                for i in xrange(0, len(seqs)):
                    seqFiles.append(self.tempFiles[1+i])
                    fileHandle = open(seqFiles[i], 'w')
                    fastaWrite(fileHandle, "%i" % i, seqs[i])
                    fileHandle.close()
                print "Have seq files ", seqFiles
            
                treeString = printBinaryTree(binaryTree, True)
                print "For tree ", treeString
                
                #align seqs and check no failure
                command = "ortheus_core -a %s -b '%s' -d %s -e" % (" ".join(seqFiles), treeString, outputFile)
                print "command to call", command
                system(command)
                
                #check alignment is complete
                alignment = [ i[:] for i in fastaAlignmentRead(outputFile) ]
                #print "alignment", alignment
                checkAlignment(alignment, seqs)
                
                print "test no is finished : %i " % test

def randomTree():
    leafNo = [-1]
    def fn():
        if random.random() > 0.6:
            return BinaryTree(random.random()*0.8, True, fn(), fn(), None)
        else:
            leafNo[0] += 1
            return BinaryTree(random.random()*0.8, False, None, None, str(leafNo[0]))
    return BinaryTree(random.random(), True, fn(), fn(), None)

def getTreeSeqs(binaryTree, seq, l):
    seq = mutateSequence(seq, binaryTree.distance)
    if binaryTree.internal:
        getTreeSeqs(binaryTree.left, seq, l)
        getTreeSeqs(binaryTree.right, seq, l)
    else:
        l.append(seq)

def checkAlignment(align, seqs):
    i = [0]*len(seqs)
    for j in align:
        for k in xrange(0, len(seqs)):
            if j[k*2] != '-':
                assert j[k*2] == seqs[k][i[k]]
                i[k] += 1
    for j in xrange(0, len(seqs)):
        assert i[j] == len(seqs[j])

def main():
    parseSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()

if __name__ == '__main__':
    main()