#!/usr/bin/env python

#Copyright (C) 2008-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt

import sys
import os

from ortheus.old.tree import BinaryTree
from ortheus.old.tree import binaryTree_depthFirstNumbers
from ortheus.old.tree import correctTreeDistances
from ortheus.old.tree import calculateProbableRootOfGeneTree

from ortheus.old.bioio import multiFastaRead, printBinaryTree
from ortheus.old.bioio import newickTreeParser
from ortheus.old.bioio import getTempFile as getTempFile_Global
from ortheus.old.bioio import logger
from ortheus.old.bioio import printMod
from ortheus.old.bioio import parseVarArgs

from ortheus.old.bioio import getOpenSeqFiles
from ortheus.old.bioio import closeSeqIterators
from ortheus.old.bioio import removeSeqFiles
from ortheus.old.bioio import concatanateSeqFiles

from ortheus.old.Nester import nestAlign
from ortheus.old.Stitcher import makePecanAlignment

"""
A script to estimate a phylogenetic tree for a set of large co-linear sequences.

Logic of script:

    Generate random binary tree
    for specified iteration number
        Call Gamma-Pecan with tree and sequences, placing output in output.mfa
        Count number of gapless columns in output.mfa
        if number of gapless columns greater > threshold:
            estimate new tree using only gapless columns and standard neighbour joining
        else:
            print warning
            estimate new tree using entire alignment and standard neighbour joining
    if branch-scaling set:
        call estimate-tree script with small sub-trees
        use average ratio to linearly restimate scaling of branch lengths.
    print tree in output file
"""

def addDefaultEstimateTreeArgs(treeArgs):
    ###user editable
    #default distance used on star tree, pretty unimportant for current Pecan (0.6)
    treeArgs.DEFAULT_DISTANCE = 0.1
    #number of iteration to do tree re-estimation
    treeArgs.ITERATION_NUMBER = 2
    treeArgs.SEMPHY_PATH="semphy"
    #used for initial tree estimation
    treeArgs.SEMPHY_ARGS_TREE="-a 4 --hky --verbose=-1 -J -H" #--posteriorDTME -O -S" #"
    treeArgs.SEMPHY_ARGS_PAIRS="-a 4 --hky --verbose=-1 -J -H -S"
    #number of alignment columns to pass to semphy for tree estimation
    treeArgs.COLUMN_MIN_GAPLESS_NO = 1000
    #branch re-estimation
    treeArgs.DO_SUBTREE_BRANCH_LENGTH_ESTIMATION = True
    #size of subtree allowed for overall rate branch lengthestimation
    treeArgs.BRANCH_LENGTH_ESTIMATION_SUBTREE_DISTANCE = 0.4
    #the species tree with which to reconcile the root of the tree
    treeArgs.SPECIES_TREE_STRING = None
    #the in-order species of the input sequences (corresponding to the names in the species tree
    treeArgs.LEAF_SPECIES = None
    
def makeAlignment(seqFiles, tree, outputAlignmentFile, alignerArgs):
    if len(seqFiles) < 30:
        makePecanAlignment(seqFiles, printBinaryTree(tree, True), outputAlignmentFile, alignerArgs)
    else:
        alignmentFile = getTempFile()
        outputScoreFile = getTempFile()
        nestAlign(tree, seqFiles, alignmentFile, outputScoreFile, alignerArgs)
        splitOutAncestors(alignmentFile,outputAlignmentFile)
        os.remove(alignmentFile)
        os.remove(outputScoreFile)
    
def getTempFile():
    return getTempFile_Global(".tree_estimate")

###probably don't look down
def parseEstimateTreeMods(mods, treeArgs, indices, skipped):
    mods.reverse()
    while mods:
        mod = mods.pop()
        if mod == '-' + indices[0]:
            treeArgs.DEFAULT_DISTANCE = float(mods.pop())
            continue
        if mod == '-' + indices[1]:
            treeArgs.ITERATION_NUMBER = int(mods.pop())
            assert treeArgs.ITERATION_NUMBER > 0
            continue
        if mod == '-' + indices[2]:
            treeArgs.SEMPHY_PATH = mods.pop()
            continue
        if mod == '-' + indices[3]:
            treeArgs.SEMPHY_ARGS_TREE = mods.pop()[1:]
            continue
        if mod == '-' + indices[4]:
            treeArgs.SEMPHY_ARGS_PAIRS = mods.pop()[1:]
            continue
        if mod == '-' + indices[5]:
            treeArgs.COLUMN_MIN_GAPLESS_NO = int(mods.pop())
            continue
        if mod == '-' + indices[6]:
            treeArgs.DO_SUBTREE_BRANCH_LENGTH_ESTIMATION = False
            continue
        if mod == '-' + indices[7]:
            treeArgs.BRANCH_LENGTH_ESTIMATION_SUBTREE_DISTANCE = float(mods.pop())
            continue
        if mod == '-' + indices[8]:
            treeArgs.SPECIES_TREE_STRING = mods.pop()
            continue
        if mod == '-' + indices[9]:
            treeArgs.LEAF_SPECIES = parseVarArgs(mods)
            continue
        skipped.append(mod)
    return indices[10:]

def printEstimateTreeMods(alignerArgs, indices):
    printMod(indices[0], '[FLOAT] default distance to use for branches in initial star tree')
    printMod(indices[1], '[INTEGER] number of iterations of tree estimation to run')
    printMod(indices[2], '[STRING] path to Semphy executable')
    printMod(indices[3], '[BRACKETED STRING STARTING WITH "#" CHARACTER] parameters to pass to Semphy for tree building')
    printMod(indices[4], '[BRACKETED STRING STARTING WITH "#" CHARACTER] parameters to pass to Semphy for pairs of sequences')
    printMod(indices[5], '[INTEGER] minimum number of columns to allow to allow for gapless tree estimations ')
    printMod(indices[6], 'Do sub-tree estimation of the rate, default: %s' % alignerArgs.DO_SUBTREE_BRANCH_LENGTH_ESTIMATION)
    printMod(indices[7], '[FLOAT] total sub tree distance to allow in subtrees used for branch length re-estimation')
    printMod(indices[8], '[STRING] species tree, if provided will guess root of tree buy minimising the number of duplications')
    printMod(indices[9], '[STRING]xN in order leaf species of the leaf sequences, needed if reconciliation-prediction of the root with the species tree is performed')
    return indices[10:]

def makeStarTree(seqNo, counter, defaultDistance):
    """
    makes binary tree, using the example of MAVIDs tree construction
    """
    if seqNo >= 2:
        return BinaryTree(defaultDistance, True, makeStarTree(seqNo/2, counter, defaultDistance), makeStarTree((seqNo/2) + (seqNo%2), counter+(seqNo/2), defaultDistance), None)
    return BinaryTree(defaultDistance, False, None, None, None)
        
def splitOutAncestors(alignmentFile, outputAlignmentFile):
    alignment = open(alignmentFile, 'r')
    out = open(outputAlignmentFile, 'w')
    counter = 0
    for line in alignment:
        if line[0] == '>':
            counter += 1
            if counter % 2 == 1:
                out.write(line)
        else:
            if counter % 2 == 1:
                out.write(line)
    out.close()
    alignment.close()
    
def formatForSemphy(alignmentFile):
    alignment = open(alignmentFile, 'r')
    outputAlignmentFile = getTempFile()
    out = open(outputAlignmentFile, 'w')
    counter = 0
    for line in alignment:
        if line[0] == '>':
            out.write(">%i\n" % counter)
            counter += 1
        else:
            i = 0
            while i < len(line)-1:
                j = 0
                while j < 100 and i < len(line)-1:
                    out.write(line[i])  
                    i += 1
                    j += 1
                out.write('\n')
    alignment.close()
    out.close()
    return outputAlignmentFile

def calculateSemphyTreeEstimate(alignmentFile, treeArgs, seqNo):
    if seqNo == 2:
        semphyArgs = treeArgs.SEMPHY_ARGS_PAIRS
    else:    
        semphyArgs = treeArgs.SEMPHY_ARGS_TREE
    semphyAlignmentFile = formatForSemphy(alignmentFile)
    outputTreeFile = getTempFile()
    characterFrequencies = " --ACGprob=%f,%f,%f" % tuple(treeArgs.EXPECTED_CHARACTER_FREQUENCIES[:-1])
    command = "%s --treeoutputfile=%s %s %s --sequence=%s " % (treeArgs.SEMPHY_PATH, outputTreeFile, semphyArgs, characterFrequencies, semphyAlignmentFile)
    #if existingTreeFile != None: #just optimise branch lengths
    #    command += " --bbl --tree=%s " % existingTreeFile
    logger.info("Calling Semphy with %s ", command)
    pipe = os.popen(command)
    if pipe.close():
        logger.info("tree building failed, so must exit")
        sys.exit(1)
    fileHandle = open(outputTreeFile, 'r')
    treeString = fileHandle.readlines()[0]
    fileHandle.close()
    binaryTree = newickTreeParser(treeString, False)
    binaryTree_depthFirstNumbers(binaryTree)
    #clean up
    os.remove(semphyAlignmentFile)
    os.remove(outputTreeFile)
    correctTreeDistances(binaryTree)
    return binaryTree

def countGaplessColumns(alignment):
    gapless = 0
    total = 0
    for column in multiFastaRead(alignment):
        if '-' not in column:
            gapless += 1
        total += 1
    return gapless, total
    
def getGaplessAlignment(alignment, seqNo):
    outputAlignment = getTempFile()
    outputFiles, outputIters = getOpenSeqFiles(seqNo, getTempFile)
    for column in multiFastaRead(alignment):
        if '-' not in column:
            for i in xrange(0, seqNo):
                outputIters[i].write(column[i])
    closeSeqIterators(outputIters, seqNo)
    concatanateSeqFiles(outputFiles, outputAlignment, seqNo, [ str(i) for i in xrange(0, seqNo) ])
    removeSeqFiles(outputFiles, seqNo)
    return outputAlignment
        
def getSubtrees(tree, totalDistance):
    subtrees = []
    def fn(tree):
        if tree.internal:
            i = fn(tree.left) + fn(tree.right)
            if i <= totalDistance:
                subtrees.append(tree)
            return i + tree.distance
        else:
            return tree.distance
    fn(tree)
    return subtrees

def getSubtreeSeqs(seqs, subTree):
    if subTree.internal:
        return getSubtreeSeqs(seqs, subTree.left) + getSubtreeSeqs(seqs, subTree.right)
    return (seqs[int(subTree.iD)],)

def labelTree(tree, counter):
    if tree.internal:
        labelTree(tree.left, counter)
        labelTree(tree.right, counter)
    else:
        tree.iD = counter()

def adjustTreeRates(tree, rateCorrection):
    tree.distance *= rateCorrection
    if tree.internal:
        adjustTreeRates(tree.left, rateCorrection)
        adjustTreeRates(tree.right, rateCorrection)

def calculateRateCorrection(tree, rateTree):
    def fn(tree):
        if tree.internal:
            i = fn(tree.left) + fn(tree.right)
            return i + tree.distance
        else:
            return tree.distance
    i = fn(rateTree.left) + fn(rateTree.right)
    j = fn(tree.left) + fn(tree.right)
    return i/j

def strCounter(i):
    def fn():
        i[0] += 1
        return str(i[0])
    return fn

def estimateTree(seqFiles, tree, iterations, doSubTreeBranchEstimation, treeArgs):
    #get sequence files
    seqNo = len(seqFiles)
    #run alignment
    treeStrings = [ printBinaryTree(tree, False) + " " + " ".join(seqFiles) ]
    for iteration in xrange(0, iterations):
        ####edit this line to set
        outputAlignment = getTempFile()
        makeAlignment(seqFiles, tree, outputAlignment, treeArgs)
        gaplessColumnNo, totalColumnNo = countGaplessColumns(outputAlignment)
        logger.info("Total number of gapless columns: %s " % gaplessColumnNo)
        if gaplessColumnNo > treeArgs.COLUMN_MIN_GAPLESS_NO: #total number of columns exceeds minimum required to do tree estimation
            gaplessOutputAlignment = getGaplessAlignment(outputAlignment, seqNo)
            tree = calculateSemphyTreeEstimate(gaplessOutputAlignment, treeArgs, seqNo)
            os.remove(gaplessOutputAlignment)
        elif totalColumnNo > 0:
            logger.info("Warning, insufficient columns to estimate tree using only gapless columns")
            tree = calculateSemphyTreeEstimate(outputAlignment, treeArgs, seqNo)
        else:
            logger.info("Warning, no alignment from which to estimate tree!!")
        logger.info("Found tree topology : %s " % printBinaryTree(tree, True))
        seqFiles = getSubtreeSeqs(seqFiles, tree)
        labelTree(tree, strCounter([-1]))
        treeString = printBinaryTree(tree, False) + " " + " ".join(seqFiles)
        logger.info("On iteration : %i , found tree and seq files (ordered) : %s " % (iteration, treeString))
        if treeString in treeStrings:
            logger.info("Topology of tree is equal to one previously seen, so exiting")
            break
        if iteration+1 < iterations:
            os.remove(outputAlignment)
        #now scale by global estimates of branch length
    if doSubTreeBranchEstimation:
        subTrees = getSubtrees(tree, treeArgs.BRANCH_LENGTH_ESTIMATION_SUBTREE_DISTANCE)
        if len(subTrees) > 0:
            rateCorrections = []
            for subTree in subTrees:
                subTree2, seqFiles2, outputAlignment2 = estimateTree(getSubtreeSeqs(seqFiles, subTree), subTree, 1, False, treeArgs)
                os.remove(outputAlignment2)
                rateCorrections.append(calculateRateCorrection(subTree, subTree2))
            for i in xrange(0, len(subTrees)):
                logger.info("Rate correction for subtree: %s %s , is calculated as : %f ", \
                            printBinaryTree(subTrees[i], True), \
                            " ".join(getSubtreeSeqs(seqFiles, subTrees[i])), rateCorrections[i])
            rateCorrection = sum(rateCorrections)/len(rateCorrections)
            logger.info("Average rate correction is calculated as : %f ", rateCorrection)
            adjustTreeRates(tree, rateCorrection)
        else:
            logger.info("No suitable branches found for rate re-estimation")    
    return tree, seqFiles, outputAlignment

def estimateTreeAlign(seqFiles, outputTreeFile, treeArgs):
    origSeqFileOrder = seqFiles[:]
    tree = makeStarTree(len(seqFiles), 0, treeArgs.DEFAULT_DISTANCE)
    binaryTree_depthFirstNumbers(tree)
    labelTree(tree, strCounter([-1]))
    tree, seqFiles, outputAlignment = estimateTree(seqFiles, tree, treeArgs.ITERATION_NUMBER, \
                                                    treeArgs.DO_SUBTREE_BRANCH_LENGTH_ESTIMATION, treeArgs)
    seqFiles = list(seqFiles)
    if treeArgs.SPECIES_TREE_STRING != None:
        logger.info("Predicting root of tree using species tree")
        speciesTree = newickTreeParser(treeArgs.SPECIES_TREE_STRING)
        binaryTree_depthFirstNumbers(speciesTree)
        logger.info("Parsed species tree: %s" % printBinaryTree(speciesTree, True))
        i = [-1]
        def fn():
            i[0] += 1
            j = origSeqFileOrder.index(seqFiles[i[0]])
            return "%s_%s" % (treeArgs.LEAF_SPECIES[j], str(i[0]))
        labelTree(tree, fn)
        tree, dupCount, lossCount = calculateProbableRootOfGeneTree(speciesTree, tree, processID=lambda x : x.split("_")[0])
        def fn2(tree):
            if tree.internal:
                fn2(tree.left)
                fn2(tree.right)
            else:
                tree.iD = tree.iD.split('_')[1]
        fn2(tree)
        seqFiles = getSubtreeSeqs(seqFiles, tree)
        logger.info("Reconciled tree with root : %s %s " % (printBinaryTree(tree, True), " ".join(seqFiles)))
        logger.info("Number of dups needed for reconcilliations : %s " % dupCount)
        logger.info("Number of losses needed for reconcilliations : %s " % lossCount)
    seqFiles = list(seqFiles)
    out = open(outputTreeFile, 'w')
    out.write("%s\n" % printBinaryTree(tree, True))
    out.write("%s\n" % " ".join(seqFiles))
    out.close()
    logger.info("Finished estimate tree")
    return tree, seqFiles, outputAlignment

def main():
    pass

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()

