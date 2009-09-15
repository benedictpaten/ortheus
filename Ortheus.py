#!/usr/bin/env python

import sys
import os
import re
import time

from ortheus.old.Nester import nestAlign
from ortheus.old.Nester import addDefaultNesterArgs
from ortheus.old.Nester import parseMods as parseModsNester
from ortheus.old.Nester import printMods as printModsNester

from ortheus.old.Stitcher import addDefaultStitcherArgs
from ortheus.old.Stitcher import parseMods as parseModsStitcher
from ortheus.old.Stitcher import printMods as printModsStitcher

from ortheus.old.EstimateTree import addDefaultEstimateTreeArgs
from ortheus.old.EstimateTree import parseEstimateTreeMods
from ortheus.old.EstimateTree import printEstimateTreeMods
from ortheus.old.EstimateTree import estimateTreeAlign

from ortheus.old.bioio import removeMod
from ortheus.old.bioio import getDefaultLogger 
from ortheus.old.bioio import newickTreeParser
from ortheus.old.bioio import getDefaultArgs
from ortheus.old.bioio import parseVarArgs
from ortheus.old.bioio import printBinaryTree
from ortheus.old.bioio import logger
from ortheus.old.bioio import loggerIndices
from ortheus.old.bioio import printFirstMods
from ortheus.old.bioio import parseFirstMods
from ortheus.old.bioio import printMod
from ortheus.old.bioio import dNAMap_IUPACToInt
from ortheus.old.bioio import fastaRead

from ortheus.old.tree import calculateCharacterFrequencies
from ortheus.old.tree import normaliseWV
from ortheus.old.tree import sumWVA

VERSION_NO="0.5.0"

def addDefaultArgs(alignerArgs):
    alignerArgs.OUTPUT_FILE = "output.mfa" 
    alignerArgs.OUTPUT_TREE_FILE = "output.newick"
    alignerArgs.OUTPUT_SCORE_FILE = "output.score" 
    alignerArgs.NEWICK_TREE_STRING = None  
    alignerArgs.EXPECTED_CHARACTER_FREQUENCIES = [0.3, 0.2, 0.2, 0.3]
    alignerArgs.EMPIRICALLY_ESTIMATE_CHARACTER_FREQUENCIES = False
    alignerArgs.MAKE_FINAL_ALIGNMENT = True

def removeReservedIndices(indices, alignerArgs):
    alignerArgs.MAKE_FINAL_ALIGNMENT_INDICE = 'y'
    removeMod(indices, alignerArgs.MAKE_FINAL_ALIGNMENT_INDICE)

def parseMods(mods, alignerArgs, indices, skipped):
    mods.reverse()
    while mods:
        mod = mods.pop()
        if mod == '-' + indices[0]:
            alignerArgs.NEWICK_TREE_STRING = mods.pop()
            continue
        if mod == '-' + indices[1]:
            alignerArgs.SEQUENCE_FILES = parseVarArgs(mods)
            continue
        if mod == '-' + indices[2]:
            alignerArgs.OUTPUT_FILE = mods.pop()
            continue
        if mod == '-' + indices[3]:
            alignerArgs.OUTPUT_TREE_FILE = mods.pop()
            continue
        if mod == '-' + indices[4]:
            alignerArgs.OUTPUT_SCORE_FILE = mods.pop()
            continue
        if mod == '-' + indices[5]:
            A = float(mods.pop())
            C = float(mods.pop())
            G = float(mods.pop())
            T = float(mods.pop())
            alignerArgs.EXPECTED_CHARACTER_FREQUENCIES = normaliseWV([ A, C, G, T ])
            continue
        if mod == '-' + indices[6]:
            alignerArgs.EMPIRICALLY_ESTIMATE_CHARACTER_FREQUENCIES = True
            continue
        if mod == '-' + alignerArgs.MAKE_FINAL_ALIGNMENT_INDICE:
            alignerArgs.MAKE_FINAL_ALIGNMENT = not alignerArgs.MAKE_FINAL_ALIGNMENT
            continue
        skipped.append(mod)
    return indices[7:]

def printMods(alignerArgs, indices):
    printMod(indices[0], '[BRACKETED STRING] set newick tree string ')
    printMod(indices[1], '[STRING]xN set sequence files (this is a mandatory argument)')
    printMod(indices[2], '[STRING] set output file (else output.mfa)')
    printMod(indices[3], '[STRING] set output file for tree (else output.newick)')
    printMod(indices[4], '[STRING] set output score file (else output.score)')
    printMod(indices[5], '[FLOAT]x[ACGT] use following nucleotide frequencies (otherwise expected frequencies have a stationary GC of 40 percent -- this is passed to all the underlying relevant methods)')
    printMod(indices[6], 'Empirically estimate stationary nucleotide frequencies from input sequences')
    printMod(alignerArgs.MAKE_FINAL_ALIGNMENT_INDICE, 'Make alignment (i.e. don\'t just do tree estimation or whatever) default : %s' % alignerArgs.MAKE_FINAL_ALIGNMENT)
    return indices[7:]

def empiricallyEstimateNucleotideFrequencies(seqFiles):
    return normaliseWV(sumWVA([ calculateCharacterFrequencies(fastaRead(seqFile), dNAMap_IUPACToInt, 4) for seqFile in seqFiles ], 4))

def main():
    sys.stderr.write("Arguments received : %s \n" % "_".join(sys.argv))
    startTime = time.time()
    alignerArgs = getDefaultArgs()
    addDefaultArgs(alignerArgs)
    addDefaultStitcherArgs(alignerArgs)
    addDefaultNesterArgs(alignerArgs)
    addDefaultEstimateTreeArgs(alignerArgs)
    i = loggerIndices
    removeReservedIndices(i, alignerArgs)
    if len(sys.argv) < 3:
        print "Ortheus.py [MODIFIER_ARGUMENTS]"
        print "Version: ", VERSION_NO
        print "A top level script for running Ortheus and Pecan to produce substitution and indel aware reconstructed chunks of genome"
        print "If you would like to contribute to this program's development please contact me at bjp (AT) ebi (DOT) ac (DOT) uk "
        print "Arguments:"
        i = printFirstMods(alignerArgs, i)
        i = printMods(alignerArgs, i)
        i = printModsStitcher(alignerArgs, i)
        i = printModsNester(alignerArgs, i)
        i = printEstimateTreeMods(alignerArgs, i)
        print "-------------Ortheus help string as follows (Changing these arguments may break the script)-------------"
        os.system("ortheus_core")
        print "-------------End Ortheus help string-------------"
        print "-------------Pecan help string as follows (Changing these arguments may break the script)-------------"
        os.system("%s bp.pecan.Pecan -help" % (alignerArgs.JAVA_PREFIX,))
        print "-------------End Pecan help string-------------"
        sys.exit(0)
        
    mods = sys.argv[1:]
    l = []
    i = parseFirstMods(mods, alignerArgs, i, l)
    i = parseMods(l, alignerArgs, i, mods)
    i = parseModsStitcher(mods, alignerArgs, i, l)
    i = parseModsNester(l, alignerArgs, i, mods)
    i = parseEstimateTreeMods(mods, alignerArgs, i, l)
    if len(l) != 0:
        logger.info("Ooops, remaining arguments %s ", " ".join(l))
        assert False  
    logger.info("Arguments received : %s " % " ".join(sys.argv))
    logger.info("Sequence files : %s " % " ".join(alignerArgs.SEQUENCE_FILES))
    if alignerArgs.EMPIRICALLY_ESTIMATE_CHARACTER_FREQUENCIES:
        alignerArgs.EXPECTED_CHARACTER_FREQUENCIES = empiricallyEstimateNucleotideFrequencies(alignerArgs.SEQUENCE_FILES)
        logger.info("Empirically estimated character frequencies : %s " % " ".join([ str(i) for i in alignerArgs.EXPECTED_CHARACTER_FREQUENCIES ]))
    try:
        os.remove(alignerArgs.OUTPUT_SCORE_FILE)
    except OSError:
        pass
    if alignerArgs.NEWICK_TREE_STRING != None:
        binaryTree = newickTreeParser(alignerArgs.NEWICK_TREE_STRING)  
        logger.info("Newick tree read : %s " % printBinaryTree(binaryTree, True))
    else:
        binaryTree, seqFiles, outputAlignment = estimateTreeAlign(alignerArgs.SEQUENCE_FILES, alignerArgs.OUTPUT_TREE_FILE, alignerArgs)
        os.remove(outputAlignment) #for now, this should be
        alignerArgs.SEQUENCE_FILES = seqFiles
    if alignerArgs.MAKE_FINAL_ALIGNMENT:
        nestAlign(binaryTree, alignerArgs.SEQUENCE_FILES, alignerArgs.OUTPUT_FILE, alignerArgs.OUTPUT_SCORE_FILE, alignerArgs)        
    #logger.info("Finished, total time taken : %s (seconds)" % (time.time()-startTime))
    print "total_time %s " % (time.time()-startTime)

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()
