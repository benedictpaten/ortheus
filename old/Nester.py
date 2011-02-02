#!/usr/bin/env python

#Copyright (C) 2008-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt

"""
A script suitable for aligning large numbers of sequences, in a progressive, nested fashion using ancestor inferral. It sub-divides the 
tree into an optimal number of sub-trees and uses Pecan and Ortheus to progressively create the final tree-alignment.
"""

import sys
import os
import time
import math

from ortheus.old.bioio import multiFastaRead
from ortheus.old.bioio import printBinaryTree
from ortheus.old.bioio import getTempFile as getTempFile_Global
from ortheus.old.bioio import logger
from ortheus.old.bioio import printMod
from ortheus.old.bioio import getOpenSeqFiles
from ortheus.old.bioio import closeSeqIterators
from ortheus.old.bioio import removeSeqFiles
from ortheus.old.bioio import concatanateSeqFiles

from ortheus.old.tree import binaryTree_nodeNames
from ortheus.old.tree import binaryTree_depthFirstNumbers
from ortheus.old.tree import felsensteins
from ortheus.old.tree import normaliseWV
from ortheus.old.tree import getBinaryTreeNodes
from ortheus.old.Stitcher import stitchAlignAndReconstruct

def addDefaultNesterArgs(alignerArgs):
    alignerArgs.MAX_SEQ_NO = 30
    alignerArgs.MAX_NODE_NO = alignerArgs.MAX_SEQ_NO*2 - 1

def parseMods(mods, alignerArgs, indices, skipped):
    mods.reverse()
    while mods:
        mod = mods.pop()
        if mod == '-' + indices[0]:
            alignerArgs.MAX_SEQ_NO = int(mods.pop())
            alignerArgs.MAX_NODE_NO = alignerArgs.MAX_SEQ_NO*2 - 1
            continue
        skipped.append(mod)
    return indices[1:]

def printMods(alignerArgs, indices):
    printMod(indices[0], '[INTEGER] maximum number of sequences to align in nested clique, default %s ' % alignerArgs.MAX_SEQ_NO)
    return indices[1:]

def makeAlignment(binaryTree, seqFiles, alignmentFile, outputScoreFile, alignerArgs):
    seqFiles = getChildSeqs(binaryTree, seqFiles)
    return stitchAlignAndReconstruct(len(seqFiles), seqFiles, printBinaryTree(binaryTree, True, False), alignmentFile, outputScoreFile, alignerArgs)

def getTempFile():
    return getTempFile_Global(".nest")

def extractSubAlignment(alignmentFile, startSeq, endSeq, newAlignmentFile):
    seqNo = endSeq-startSeq
    outputFiles, outputIters = getOpenSeqFiles(seqNo, getTempFile)
    for column in multiFastaRead(alignmentFile):
        for i in xrange(startSeq, endSeq):
            if column[i] != '-':
                for j in xrange(startSeq, endSeq):
                    outputIters[j-startSeq].write(column[j])
                break
    closeSeqIterators(outputIters, seqNo)
    concatanateSeqFiles(outputFiles, newAlignmentFile, seqNo, [ str(i) for i in xrange(startSeq, endSeq) ])
    removeSeqFiles(outputFiles, seqNo)

def outputMergedColumn(outputIters, i, j, k):
    m = 0
    for l in i:
        outputIters[m].write(l)
        m += 1
    for l in j:
        outputIters[m].write(l)
        m += 1
    for l in k:
        outputIters[m].write(l)
        m += 1

def mergeAlignments(daddyAlignmentFile, childAlignmentFile, childSeq, childSeq2, daddySeqNo, childSeqNo, newAlignmentFile, nodeLabels):
    childGapColumn = ['-']*childSeqNo
    topDaddyGapColumn = ['-']*childSeq
    bottomDaddyGapColumn = ['-']*(daddySeqNo-(childSeq+1))
    
    totalSeqNo = daddySeqNo -1 + childSeqNo
    outputFiles, outputIters = getOpenSeqFiles(totalSeqNo, getTempFile)
    
    childIter = multiFastaRead(childAlignmentFile)
    def nextChild():
        try:
            return childIter.next()
        except StopIteration:
            return None
    childColumn = nextChild()
    for daddyColumn in multiFastaRead(daddyAlignmentFile):
        if daddyColumn[childSeq] != '-':
            while childColumn[childSeq2] == '-':
                outputMergedColumn(outputIters, topDaddyGapColumn, childColumn, bottomDaddyGapColumn)
                childColumn = nextChild()
            outputMergedColumn(outputIters, daddyColumn[:childSeq], childColumn, daddyColumn[childSeq+1:])
            childColumn = nextChild()
        else:
            outputMergedColumn(outputIters, daddyColumn[:childSeq], childGapColumn, daddyColumn[childSeq+1:])
    while childColumn != None:
        outputMergedColumn(outputIters, topDaddyGapColumn, childColumn, bottomDaddyGapColumn)
        childColumn = nextChild()
    closeSeqIterators(outputIters, totalSeqNo)
    concatanateSeqFiles(outputFiles, newAlignmentFile, totalSeqNo, nodeLabels[0:totalSeqNo])
    removeSeqFiles(outputFiles, totalSeqNo)

def mergeTogetherAllAlignments(binaryTree, alignmentFiles, nodeLabels, subTreeCounter):
    subTreeCounter[0] += 1
    if binaryTree.internal:
        alignmentFile = alignmentFiles[binaryTree.traversalID.mid]
        if alignmentFile != None:
            subTreeCounter=[1]
            childAlignments = mergeTogetherAllAlignments(binaryTree.left, alignmentFiles, nodeLabels, subTreeCounter) + mergeTogetherAllAlignments(binaryTree.right, alignmentFiles, nodeLabels, subTreeCounter)
            for childAlignmentFile, childTree in childAlignments:
                mergeAlignments(alignmentFile, childAlignmentFile, childTree.traversalID.midStart - binaryTree.traversalID.midStart, \
                                childTree.traversalID.mid - childTree.traversalID.midStart, \
                                subTreeCounter[0], childTree.traversalID.midEnd - childTree.traversalID.midStart, alignmentFile, nodeLabels)
                subTreeCounter[0] += childTree.traversalID.midEnd - childTree.traversalID.midStart - 1
                os.remove(childAlignmentFile)
            return ((alignmentFile, binaryTree),)
        return mergeTogetherAllAlignments(binaryTree.left, alignmentFiles, nodeLabels, subTreeCounter) + mergeTogetherAllAlignments(binaryTree.right, alignmentFiles, nodeLabels, subTreeCounter)
    return ()
    
def getChildSeqs(binaryTree, seqFiles):
    if binaryTree == None:
        return ()
    seqFile = seqFiles[binaryTree.traversalID.mid]
    if seqFile != None:
        return (seqFile,)
    return getChildSeqs(binaryTree.left, seqFiles) + getChildSeqs(binaryTree.right, seqFiles)

def calculateTreeNodeCosts(binaryTree, alpha=1.5):
    def avg(binaryTree, i, j):
        if binaryTree.internal:
            return (i[binaryTree.left.traversalID.mid] - normaliseWV(j[binaryTree.left.traversalID.mid])[0]) + \
            (i[binaryTree.right.traversalID.mid] - normaliseWV(j[binaryTree.right.traversalID.mid])[0])
        return 0.0
    def subMatrix(distance):
        i = 0.5 + 0.5*math.exp(-2.0*distance*alpha)
        return ( (i, 1.0 - i), (1.0 - i, i) )
    binaryTreeNodes = []
    getBinaryTreeNodes(binaryTree, binaryTreeNodes)
    subMatrices = [ subMatrix(i.distance) for i in binaryTreeNodes ]
    leaves = [ (1.0, 0.0) for i in xrange(0, binaryTree.traversalID.midEnd) ]
    j = felsensteins(binaryTree, subMatrices, (0.5, 0.5), leaves, 2)
    for i in j.keys():
        j[i] = normaliseWV(j[i])[0]
    finalCosts = [ avg(i, j, felsensteins(i, subMatrices, (0.5, 0.5), leaves, 2)) for i in binaryTreeNodes ]
    return finalCosts

def calculatePath(binaryTree, costs, maxSeqNo):
    pointers = {}
    matrix = {}
    def fn(binaryTree):
        iD = binaryTree.traversalID.mid
        if binaryTree.internal:
            fn(binaryTree.left)
            fn(binaryTree.right)
            #the action
            matrixLeft = matrix[binaryTree.left.traversalID.mid]
            matrixRight = matrix[binaryTree.right.traversalID.mid]
            pointersLeft = pointers[binaryTree.left.traversalID.mid]
            pointersRight = pointers[binaryTree.right.traversalID.mid]
            matrixNode = {}
            pointersNode = {}
           
            for i in matrixLeft.keys():
                for j in matrixRight.keys():
                    k = matrixLeft[i] + matrixRight[j]
                    if (not matrixNode.has_key(i + j + 1)) or matrixNode[i + j + 1] > k:
                        matrixNode[i + j + 1] = k
                        pointersNode[i + j + 1] = pointersLeft[i] + pointersRight[j]
            j = sys.maxint
            k = None
            for i in xrange(3, maxSeqNo+1):
                if matrixNode.has_key(i) and matrixNode[i] < j:
                    j = matrixNode[i]
                    k = pointersNode[i]
            j += costs[iD]
            if (not matrixNode.has_key(3)) or j < matrixNode[3]: #don't replace 3 way internal nodes with two leaves attached
                matrixNode[3] = j
                pointersNode[3] = k + (binaryTree,)
            matrix[iD] = matrixNode
            pointers[iD] = pointersNode
        else:
            matrix[iD] = { 1:0.0 }
            pointers[iD] = { 1:() }
    fn(binaryTree)
    j = sys.maxint
    k = sys.maxint
    iD = binaryTree.traversalID.mid
    matrixNode = matrix[iD]
    for i in matrixNode.keys():
        if i > 3 and i < maxSeqNo+1 and matrixNode[i] < j:
            j = matrixNode[i]
            k = i
    if j != sys.maxint:
        return (j, pointers[iD][k])
    return (0.0, ())

def removeInternalIDs(binaryTree):
    """Removes the ids from internal nodes.
    """
    if binaryTree.internal == True:
        binaryTree.iD = None
    if binaryTree.left:
        removeInternalIDs(binaryTree.left)
        removeInternalIDs(binaryTree.right)

def nestAlign(binaryTree, leafSeqFiles, outputFile, outputScoreFile, alignerArgs):
    logger.info("Starting Nester")
    maxNodeNo = alignerArgs.MAX_NODE_NO
    
    removeInternalIDs(binaryTree)
    
    logger.info("Binary tree : %s " % printBinaryTree(binaryTree, True, False))
    binaryTree_depthFirstNumbers(binaryTree)
    nodeNo = binaryTree.traversalID.midEnd
    logger.info("Labelled tree with numbers ")
    
    seqNo = len(leafSeqFiles)
    logger.info(" Sequence files : %s" % " ".join(leafSeqFiles))
    #assert seqNo*2 - 1 == nodeNo
    
    logger.info("Output file %s " % outputFile)
    
    labels = binaryTree_nodeNames(binaryTree)
    costs = calculateTreeNodeCosts(binaryTree)
    logger.info("Calculated node costs")
    for node in xrange(0, nodeNo):
        logger.info("Node : %s , reconstruction value : %f , %f" % (labels[node], costs[node], 1.0 - costs[node]))
    pathCost, treePath = calculatePath(binaryTree, costs, maxNodeNo)
    logger.info(" Calculated nested path. Cost : %f , Path : %s" % (pathCost, " ".join([ labels[i.traversalID.mid] for i in treePath ])))
    assert len(leafSeqFiles) == seqNo
    alignmentFiles = [None] * nodeNo
    seqFiles = [None] * nodeNo
    for i in xrange(0, seqNo):
        seqFiles[i*2] = leafSeqFiles[i]
    logger.debug("About to start main nested loop")
    for subTree in treePath:
        assert subTree != binaryTree
        logger.info("Chosen sub tree to align : %s " % printBinaryTree(subTree, True, False))
        alignmentFile = getTempFile()
        startTime = time.time()
        makeAlignment(subTree, seqFiles, alignmentFile, outputScoreFile, alignerArgs)
        logger.info("Made alignment of subtree, time taken : %s (seconds)" % (time.time()-startTime))
        #get the two ancestors
        subTreeTraversalIDs = binaryTree_depthFirstNumbers(subTree, labelTree=False, dontStopAtID=False)
        
        if subTree.left.internal:
            offset = subTreeTraversalIDs[subTree].midStart
            childXAlignmentFile = getTempFile()
            extractSubAlignment(alignmentFile, 0, subTreeTraversalIDs[subTree].mid-offset, childXAlignmentFile)
            alignmentFiles[subTree.left.traversalID.mid] = childXAlignmentFile
            logger.info("Extracted alignment of left child : %s " % printBinaryTree(subTree.left, True, False))
            
            assert offset == subTreeTraversalIDs[subTree.left].midStart
            childXSeqFile = getTempFile()
            extractSubAlignment(childXAlignmentFile, subTreeTraversalIDs[subTree.left].mid - offset, subTreeTraversalIDs[subTree.left].mid - offset + 1, childXSeqFile)
            seqFiles[subTree.left.traversalID.mid] = childXSeqFile
            logger.info("Extracted sequence of left child : %s " % printBinaryTree(subTree.left, True, False))
        
        if subTree.right.internal:
            offset = subTreeTraversalIDs[subTree].midStart
            childYAlignmentFile = getTempFile()
            extractSubAlignment(alignmentFile, subTreeTraversalIDs[subTree].mid + 1 - offset, subTreeTraversalIDs[subTree].midEnd - offset, childYAlignmentFile)
            alignmentFiles[subTree.right.traversalID.mid] = childYAlignmentFile  
            logger.info("Extracted alignment of right child : %s " % printBinaryTree(subTree.right, True, False))
            
            offset = subTreeTraversalIDs[subTree.right].midStart
            childYSeqFile = getTempFile()
            extractSubAlignment(childYAlignmentFile, subTreeTraversalIDs[subTree.right].mid - offset, subTreeTraversalIDs[subTree.right].mid - offset + 1, childYSeqFile)
            seqFiles[subTree.right.traversalID.mid] = childYSeqFile  
            logger.info("Extracted sequence of right child : %s " % printBinaryTree(subTree.right, True, False))
        
        subTree.left.iD = labels[subTree.left.traversalID.mid] #labels tree, so we only print relevant bits
        subTree.right.iD = labels[subTree.right.traversalID.mid]
        os.remove(alignmentFile)
        logger.info("Finished loop and reduced tree to : %s " % printBinaryTree(subTree, True, False))
    startTime = time.time()
    makeAlignment(binaryTree, seqFiles, outputFile, outputScoreFile, alignerArgs)
    logger.info("Finished final nested alignment, time taken : %s (seconds)" % (time.time()-startTime))
    alignmentFiles[binaryTree.traversalID.mid] = outputFile
    mergeTogetherAllAlignments(binaryTree, alignmentFiles, labels, [0])
    logger.info("Merged together all alignments")
    for i in xrange(1, nodeNo, 2):
        if seqFiles[i] != None:
            os.remove(seqFiles[i])
    removeInternalIDs(binaryTree)
    logger.info("Have cleaned up, and am returning")
    
def main():
    pass

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()
