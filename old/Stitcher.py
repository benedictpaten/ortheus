#!/usr/local/bin/python2.3

"""
Runs Pecan, then Ortheus upon output, and creates ancestor sequences
"""

import sys
import os
import os.path
import time

from ortheus.old.bioio import multiFastaRead
from ortheus.old.bioio import getMultiFastaOffsets
from ortheus.old.bioio import getTempFile as getTempFile_Global
from ortheus.old.bioio import concatanateSeqFiles
from ortheus.old.bioio import logger
from ortheus.old.bioio import printMod
from ortheus.old.bioio import newickTreeParser
from ortheus.old.bioio import printBinaryTree
from ortheus.old.bioio import getOpenSeqFiles
from ortheus.old.bioio import closeSeqIterators
from ortheus.old.bioio import removeSeqFiles

from ortheus.old.tree import binaryTree_depthFirstNumbers
from ortheus.old.tree import binaryTree_nodeNames

#########################################################
#########################################################
#########################################################
#global output settings
#########################################################
#########################################################
#########################################################

#aligner
def addDefaultStitcherArgs(alignerArgs):
    alignerArgs.JAVA_PREFIX = "java -server "
    alignerArgs.ALIGNER_PREFIX =  " bp.pecan.Pecan "
    #alignerArgs.ALIGNER_PREFIX =  " bp.pecan.Pecan"
    alignerArgs.RECONSTRUCTION_PREFIX = "ortheus_core"
    alignerArgs.ALIGNMENT_ARGS = " " #-X -d -q -r 1.0 "
    alignerArgs.RECONSTRUCTION_ARGS = " "
    alignerArgs.ALIGNMENT_ARGS_FAST = " " #-X -d -q -r 1.0 "
    alignerArgs.RECONSTRUCTION_ARGS_FAST = " -i 10 -j 0 "
    alignerArgs.CAUTIOUS_ARGS = " -j 0 "
    alignerArgs.ALIGNMENT_CHUNK_MAX_COLUMN_SIZE = 2000
    alignerArgs.VITERBI_ALIGNMENT_COLUMN_GAP = 200
    alignerArgs.FAST_SETTING = False
    
def parseMods(mods, alignerArgs, indices, skipped):
    mods.reverse()
    while mods:
        mod = mods.pop()
        if mod == '-' + indices[0]:
            s = " " + mods.pop()[1:] + " "
            alignerArgs.ALIGNMENT_ARGS += s
            alignerArgs.ALIGNMENT_ARGS_FAST += s
            continue
        if mod == '-' + indices[1]:
            s = " " + mods.pop()[1:] + " "
            alignerArgs.RECONSTRUCTION_ARGS += s
            alignerArgs.RECONSTRUCTION_ARGS_FAST += s
            continue
        if mod == '-' + indices[2]:
            alignerArgs.JAVA_PREFIX = mods.pop()
            continue
        if mod == '-' + indices[3]:
            alignerArgs.ALIGNMENT_CHUNK_MAX_COLUMN_SIZE = int(mods.pop())
            continue
        if mod == '-' + indices[4]:
            alignerArgs.VITERBI_ALIGNMENT_COLUMN_GAP = int(mods.pop())
            continue
        skipped.append(mod)
    return indices[5:]

def printMods(alignerArgs, indices):
    printMod(indices[0], '[BRACKETED STRING STARTING WITH "#" CHARACTER] Extra command line arguments to pass to Pecan (SETTING SOME OPTIONS WILL BE INCOMPATIBLE WITH THESE SCRIPTS)')
    printMod(indices[1], '[BRACKETED STRING STARTING WITH "#" CHARACTER] Extra command line arguments to pass to Ortheus (SETTING SOME OPTIONS WILL BE INCOMPATIBLE WITH THESE SCRIPTS)')
    printMod(indices[2], '[BRACKETED STRING] Java string (i.e. java -server)')
    printMod(indices[3], '[INTEGER] maximum number of alignment columns in chunk for Ortheus reconstruction')
    printMod(indices[4], '[INTEGER] number of columns to leave between end of reconstruction chunk and final reconstruction')
    return indices[5:]

def getTempFile():
    return getTempFile_Global(".stitch")

def getNextAlignmentChunk(previousAlignment, alignment, size, seqNo, labels):
    """
    outputs fragment of multiple alignment, and individual sequence files
    """
    seqFiles, seqIterators = getOpenSeqFiles(seqNo, getTempFile)
    for seq in xrange(0, seqNo):
        seqIterators[seq].write(">\n")
    alignmentFiles, alignmentIterators = getOpenSeqFiles(seqNo, getTempFile)
    columnCount = 0
    end = False
    for column in previousAlignment:
        assert len(column) == seqNo
        assert len(column) != ['-']*seqNo
        if columnCount >= size:
            break
        for seq in xrange(0, seqNo):
            residue = column[seq]
            alignmentIterators[seq].write(residue)
            if column[seq] != '-':
                seqIterators[seq].write(residue)
        columnCount += 1
    else:
        for column in alignment:
            assert column != None
            assert len(column) == seqNo
            assert len(column) != ['-']*seqNo
            previousAlignment.append(column[:])
            if columnCount >= size:
                break
            for seq in xrange(0, seqNo):
                residue = column[seq]
                alignmentIterators[seq].write(residue)
                if column[seq] != '-':
                    seqIterators[seq].write(residue)
            columnCount += 1
        else:
            end = True
    closeSeqIterators(alignmentIterators, seqNo)
    closeSeqIterators(seqIterators, seqNo)
    tempAlignmentFile = getTempFile()
    concatanateSeqFiles(alignmentFiles, tempAlignmentFile, seqNo, labels)
    removeSeqFiles(alignmentFiles, seqNo)
    if columnCount > 0:
        return seqFiles, tempAlignmentFile, end
    removeSeqFiles(seqFiles, seqNo)
    os.remove(tempAlignmentFile)
    return None, None, None

def removeFromLeft(completedAlignment, alignment, nodeNo, seqNo):
    indices = [0]*seqNo
    for column in completedAlignment:
        assert len(column) == nodeNo
        for i in xrange(0, seqNo):
            if column[i*2] != '-':
                indices[i] += 1
    #logger.debug("Indices of sequences aligned : %s ", " ".join([ str(i) for i in indices ]))
    if indices == [0]*seqNo:
        print "nnnnnnnoooooo"
        sys.exit(1)
    l = []
    for i in xrange(0, len(alignment)):
        column = alignment[i]
        assert len(column) == seqNo
        gapCount = 0
        for j in xrange(0, seqNo):
            if column[j] != '-' :
                if indices[j] > 0:
                    gapCount += 1
                    column[j] = '-'
                    indices[j] -= 1
            else:
                gapCount += 1
        if gapCount != seqNo:
            l.append(column)
    return l
    
def appendToAlignment(alignmentIter, outputIter, seqNo):
    for column in alignmentIter:
        assert len(column) == seqNo
        for seq in xrange(0, seqNo):
            outputIter[seq].write(column[seq])
            
def appendScore(scoreFile, previousScoreFile):
    i = open(scoreFile, 'r')
    j = float(i.readline())
    i.close()
    try:
        i = open(previousScoreFile, 'r')
        line = i.readline()
        if line != '':
            k = float(line)
        else:
            k = 0.0
        i.close()
    except IOError:
        k = 0.0
    i = open(previousScoreFile, 'w')
    i.write("%f\n" % (j + k))
    i.close()
    
def stitchAlignAndReconstruct(seqNo, inputSeqFiles, treeString, outputFile, outputScoreFile, alignerArgs):
    alignmentFile = getTempFile()
    makePecanAlignment(inputSeqFiles, treeString, alignmentFile, alignerArgs)
    stitchReconstruct(seqNo, inputSeqFiles, treeString, outputFile, outputScoreFile, alignmentFile, alignerArgs)
    os.remove(alignmentFile)
    
def makePecanAlignment(inputSeqFiles, treeString, alignmentFile, alignerArgs):
    if alignerArgs.FAST_SETTING:
        alignmentArgs = alignerArgs.ALIGNMENT_ARGS_FAST
    else:
        alignmentArgs = alignerArgs.ALIGNMENT_ARGS
    pecanTime = time.time()
    command = "%s %s -F %s -E '%s' -G %s %s " % (alignerArgs.JAVA_PREFIX, alignerArgs.ALIGNER_PREFIX, " ".join(inputSeqFiles), treeString, alignmentFile, alignmentArgs)
    logger.info("Calling Pecan with : %s", command)
    if os.system(command):
        print "Something went wrong calling aligner, so I've got to go"
        sys.exit(1)
    logger.info("Completed alignment in : %s (seconds)" % (time.time()-pecanTime))

def stitchReconstruct(seqNo, inputSeqFiles, treeString, outputFile, outputScoreFile, inputAlignmentFile, alignerArgs):
    startTime = time.time() #epoch time in seconds
    
    logger.info("Starting Stitcher")
    reconstructionPrefix = alignerArgs.RECONSTRUCTION_PREFIX
    if alignerArgs.FAST_SETTING:
        reconstructionArgs = alignerArgs.RECONSTRUCTION_ARGS_FAST
    else:
        reconstructionArgs = alignerArgs.RECONSTRUCTION_ARGS
    cautiousArgs = alignerArgs.CAUTIOUS_ARGS
    alignmentChunkMaxSeqSize = alignerArgs.ALIGNMENT_CHUNK_MAX_COLUMN_SIZE
    viterbiAlignmentColumnGap = alignerArgs.VITERBI_ALIGNMENT_COLUMN_GAP
    #parse tree 
    binaryTree = newickTreeParser(treeString)
    binaryTree_depthFirstNumbers(binaryTree)
    logger.info("Newick tree read : %s " % printBinaryTree(binaryTree, True))
    labels = binaryTree_nodeNames(binaryTree)
    leafLabels = [ labels[i] for i in xrange(0, len(labels)) if (i%2) == 0]
    #load alignment iterator
    alignmentReader = multiFastaRead(inputAlignmentFile, lambda x : x)
    #number of sequences, including ancestors
    nodeNumber = binaryTree.traversalID.midEnd
    assert nodeNumber == seqNo * 2 - 1
    #create output files
    outputFiles, outputIterators = getOpenSeqFiles(nodeNumber, getTempFile)
    #while has chunk
    previousAlignment = []
    alignmentSeqs, alignmentFile, end = getNextAlignmentChunk(previousAlignment, alignmentReader, alignmentChunkMaxSeqSize, seqNo, leafLabels)
    tempTreeStatesFile = getTempFile()
    loopOptions = " "  
    logger.info("Starting main loop")
    characterFrequenciesString = " ".join([ str(i) for i in alignerArgs.EXPECTED_CHARACTER_FREQUENCIES ])
    while alignmentSeqs != None:
        if(end):
            viterbiAlignmentColumnGap = 0
        tempAncestorFile = getTempFile()
        tempScoreFile = getTempFile()
        command = "%s -b '%s' -c %s -a %s -u %s -s %s %s %s -d %s -n %s -x %s " % (reconstructionPrefix, treeString, alignmentFile, \
                                                                       " ".join(alignmentSeqs), tempTreeStatesFile, \
                                                                       viterbiAlignmentColumnGap, loopOptions, reconstructionArgs, tempAncestorFile, characterFrequenciesString, tempScoreFile)
        logger.info("Calling Ortheus with : %s", command)
        exitValue = os.system(command)
        if exitValue != 0:
            logger.info("Something went wrong calling Ortheus : %i ", exitValue)
            #if exitValue != 73:
            #    logger.info("Unrecognised issue, so am exiting to be cautious")
            #    sys.exit(1)
            logger.info("Going to retry with caution settings")
            command = "%s -b '%s' -c %s -a %s -u %s -s %s %s %s -d %s -x %s" % (reconstructionPrefix, treeString, alignmentFile, \
                                                                       " ".join(alignmentSeqs), tempTreeStatesFile, \
                                                                       viterbiAlignmentColumnGap, loopOptions, cautiousArgs, tempAncestorFile, tempScoreFile)
            logger.info("Calling Ortheus with : %s", command)
            if os.system(command):
                logger.info("Already tried caution, so have to go")
                sys.exit(1)
        logger.info("Completed reconstruction of chunk")
        appendScore(tempScoreFile, outputScoreFile)
        os.remove(tempScoreFile)
        loopOptions = " -t " + tempTreeStatesFile
        tempAncestorFastaOffsets = getMultiFastaOffsets(tempAncestorFile)
        previousAlignment = removeFromLeft(multiFastaRead(tempAncestorFile, lambda x : x, tempAncestorFastaOffsets), previousAlignment, nodeNumber, seqNo)
        appendToAlignment(multiFastaRead(tempAncestorFile, lambda x : x, tempAncestorFastaOffsets), outputIterators, nodeNumber)
        logger.info("Added reconstructed chunk to complete alignment")
        os.remove(tempAncestorFile)
        removeSeqFiles(alignmentSeqs, seqNo)
        os.remove(alignmentFile)
        logger.info("Cleaned up at end of loop")
        alignmentSeqs, alignmentFile, end = getNextAlignmentChunk(previousAlignment, alignmentReader, alignmentChunkMaxSeqSize, seqNo, leafLabels)
    logger.info("Finished main loop")
    #load into single output file
    closeSeqIterators(outputIterators, nodeNumber)
    concatanateSeqFiles(outputFiles, outputFile, nodeNumber, labels)
    logger.info("Written out alignment to single file")
    #clean up
    os.remove(tempTreeStatesFile)
    removeSeqFiles(outputFiles, nodeNumber)
    logger.info("Cleaned up final files")
    logger.info("Finished, total time taken for stitcher: %s (seconds)" % (time.time()-startTime))

def main():
    pass

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()