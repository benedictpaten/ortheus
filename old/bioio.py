#!/usr/bin/env python

#Copyright (C) 2008-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt

import sys
import os
import re
import logging
import tempfile

#from tree import BinaryTree
import ortheus.old.tree

DEFAULT_DISTANCE = 0.001

#########################################################
#########################################################
#########################################################
#global logging settings
#########################################################
#########################################################
#########################################################

loggingFormatter = logging.Formatter('%(asctime)s %(levelname)s %(lineno)s %(message)s')

def getDefaultLogger(level=logging.INFO):
    logger = logging.getLogger()
    logger.setLevel(level)
    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(logging.CRITICAL) #null logger, to stop annoying error message
    logger.addHandler(handler)
    return logger

logger = getDefaultLogger()

loggerIndices = [ 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
                  'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z' ]

LOG_FILE = "Ortheus.log"


def getDefaultArgs():
    class a:
        pass
    return a()

def printMod(indice, string):
    print '\t-' + indice + ' ' + string

def removeMod(indices, indice):
    indices.remove(indice)

def parseFirstMods(mods, alignerArgs, indices, skipped):
    mods.reverse()
    while mods:
        mod = mods.pop()
        if mod == '-' + indices[0]:
            logger.setLevel(logging.DEBUG)
            continue
        if mod == '-' + indices[1]:
            handler = logging.StreamHandler(sys.stderr)
            handler.setFormatter(loggingFormatter)
            logger.addHandler(handler)
            continue
        if mod == '-' + indices[2]:
            handler = logging.FileHandler(LOG_FILE)
            handler.setFormatter(loggingFormatter)
            logger.addHandler(handler)
            continue
        skipped.append(mod)
    return indices[3:]

def printFirstMods(alignerArgs, indices):
    printMod(indices[0], 'set logging level to DEBUG (default is INFO, which does not include low level debug stuff)')
    printMod(indices[1], 'log to standard error')
    printMod(indices[2], 'log to an output file in the current directory (ortheus.log)')
    return indices[3:]

def parseVarArgs(mods):
    """
    gets variable numbers of arguments before '-'
    """
    i = len(mods)-1
    while i >= 0 and mods[i][0] != '-':
        i -= 1
    l = mods[i+1:]
    l.reverse()
    for i in xrange(0, len(l)):
        mods.pop()
    return l

#########################################################
#########################################################
#########################################################
#misc parsing functions
#########################################################
#########################################################
#########################################################

def getNextNonCommentLine(file):
    line = file.readline()
    while line[0] == '#':
        line = file.readline()
    return line

def parseIntLine(line):
    return [ int(i) for i in line.split() ]

#########################################################
#########################################################
#########################################################
#fasta functions
#########################################################
#########################################################
#########################################################

dNAMap_IUPACToInt_Local = { 'A':0, 'C':1, 'G':2, 'T':3, 'a':0, 'c':1, 'g':2, 't':3, 'N':4, 'n':4 }

dNAMap_IntToIUPAC_Local = { 0:'A', 1:'C', 2:'G', 3:'T', 4:'N' }

def dNAMap_IUPACToInt(i):
    return dNAMap_IUPACToInt_Local[i]

def dNAMap_IntToIUPAC(i):
    return dNAMap_IntToIUPAC_Local[i]

def dNAMap_IUPACToWVFn(char):
    dNAMap_IntToWV = [ 
         ( 1.0, 0.0, 0.0, 0.0 ), 
         ( 0.0, 1.0, 0.0, 0.0 ), 
         ( 0.0, 0.0, 1.0, 0.0 ), 
         ( 0.0, 0.0, 0.0, 1.0 ),
         ( 0.25, 0.25, 0.25, 0.25 ) 
         ]
    dNAMap_IUPACToInt = { 'A':0, 'C':1, 'G':2, 'T':3, \
           'a':0, 'c':1, 'g':2, 't':3, \
           'N':4, 'n':4 }
    return dNAMap_IntToWV[dNAMap_IUPACToInt[char]]

def fastaRead_MultipleSequences(fasta, map=lambda x:x):
    """
    ignores all '>' lines and sucks up everything else, bar newline characters
    """
    j = open(fasta, 'r')
    seqs = []
    seq = None
    for i in j:
        if i[0] != '>':
            seq = seq + [ map(k) for k in i[:-1] ]
        else:
            if seq != None:
                seqs.append(seq)
            seq = []
    if seq != None:
        seqs.append(seq)
    j.close()
    return seqs

def fastaRead(fasta, map=lambda x:x):
    """
    ignores all '>' lines and sucks up everything else, bar newline characters
    """
    j = open(fasta, 'r')
    seq = [ map(i) for i in "".join([ i[:-1] for i in j if i[0] != '>' ]) ]
    j.close()
    return seq

def getMultiFastaOffsets(fasta):
    """
    reads in columns of multiple alignment and returns them iteratively
    """
    f = open(fasta, 'r')
    i = 0
    j = f.read(1)
    l = []
    while j != '':
        i += 1
        if j == '>':
            i += 1
            while f.read(1) != '\n':
                i += 1
            l.append(i)
        j = f.read(1)
    f.close()
    return l

def multiFastaRead(fasta, map=lambda x : x, l=None):
    """
    reads in columns of multiple alignment and returns them iteratively
    """
    if l == None:
        l = getMultiFastaOffsets(fasta)
    else:
        l = l[:]
    seqNo = len(l)
    for i in xrange(0, seqNo):
        j = open(fasta, 'r')
        j.seek(l[i])
        l[i] = j
    column = [sys.maxint]*seqNo
    if seqNo != 0:
        while True:
            for j in xrange(0, seqNo):
                i = l[j].read(1)
                while i == '\n':
                    i = l[j].read(1)
                column[j] = i
            if column[0] == '>' or column[0] == '':
                for j in xrange(1, seqNo):
                    assert column[j] == '>' or column[j] == ''
                break
            for j in xrange(1, seqNo):
                 assert column[j] != '>' and column[j] != ''
                 column[j] = map(column[j])
            yield column[:]
    for i in l:
        i.close()
        
def writeFastaFile(seqs, names, seqNo, fastaFile, 
                   filter=lambda x : x):
    """
    Writes out column alignment to given file multi-fasta format
    """
    fastaFile = open(fastaFile, 'w')
    for seq in xrange(0, seqNo):
        fastaFile.write(">%s\n" % names[seq])
        for base in seqs[seq]:
            fastaFile.write(filter(base))
        fastaFile.write("\n")
    fastaFile.close()

def writeFastaAlignment(columnAlignment, names, seqNo, fastaFile, 
                        filter=lambda x : True):
    """
    Writes out column alignment to given file multi-fasta format
    """
    fastaFile = open(fastaFile, 'w')
    columnAlignment = [ i for i in columnAlignment if filter(i) ]
    for seq in xrange(0, seqNo):
        fastaFile.write(">%s\n" % names[seq])
        for column in columnAlignment:
            fastaFile.write(column[seq])
        fastaFile.write("\n")
    fastaFile.close()
    
def getTempFile(suffix=""):
    """
    returns a string representing a temporary file, that must be manually deleted
    """
    handle, file = tempfile.mkstemp(suffix)
    os.close(handle)
    return file

def concatanateSeqFiles(seqFiles, outputFile, seqNo, nodeNames):
    """
    concatenate sequence files into one file
    """
    #turn into one alignment
    outputFile = open(outputFile, 'w')
    for seq in xrange(0, seqNo):
        outputFile.write(">%s\n" % nodeNames[seq])
        j = open(seqFiles[seq], 'r')
        for i in j:
            outputFile.write(i)
        j.close()
        outputFile.write("\n")
    outputFile.close()
    
def getOpenSeqFiles(seqNo, getTempFile):
    seqFiles = [ getTempFile() for seq in xrange(0, seqNo) ]
    seqIterators = [ open(seqFiles[seq], 'w') for seq in xrange(0, seqNo) ]
    return (seqFiles, seqIterators)

def closeSeqIterators(seqIterators, seqNo):
    for seq in xrange(0, seqNo):
        seqIterators[seq].close();
    
def removeSeqFiles(seqFiles, seqNo):
    for seq in xrange(0, seqNo):
        os.remove(seqFiles[seq])
        
#########################################################
#########################################################
#########################################################
#newick tree functions
#########################################################
#########################################################
#########################################################
       
def newickTreeParser(newickTree, defaultDistance=DEFAULT_DISTANCE, \
                     sortNonBinaryNodes=False, reportUnaryNodes=False):
    """
    lax newick tree parser
    """
    newickTree = newickTree.replace("(", " ( ")
    newickTree = newickTree.replace(")", " ) ")
    newickTree = newickTree.replace(":", " : ")
    newickTree = newickTree.replace(";", "")
    newickTree = newickTree.replace(",", " , ")
    
    newickTree = re.compile("[\s]*").split(newickTree)
    while "" in newickTree:
        newickTree.remove("")
    def fn(newickTree, i):
        if i[0] < len(newickTree):
            if newickTree[i[0]] == ':':
                d = float(newickTree[i[0]+1])
                i[0] += 2
                return d
        return defaultDistance
    def fn2(newickTree, i):
        if i[0] < len(newickTree):
            j = newickTree[i[0]]
            if j != ':' and j != ')' and j != ',':
                i[0] += 1
                return j
        return None
    def fn3(newickTree, i):
        if newickTree[i[0]] == '(':
            #subTree1 = None
            subTreeList = []
            i[0] += 1
            k = []
            while newickTree[i[0]] != ')':
                if newickTree[i[0]] == ',':
                    i[0] += 1
                subTreeList.append(fn3(newickTree, i))
            i[0] += 1
            def cmp(i, j):
                if i.distance < j.distance:
                    return -1
                if i.distance > j.distance:
                    return 1
                return 0
            if sortNonBinaryNodes:
                subTreeList.sort(cmp)
            subTree1 = subTreeList[0]
            if len(subTreeList) > 1:
                for subTree2 in subTreeList[1:]:
                    subTree1 = ortheus.old.tree.BinaryTree(0.0, True, subTree1, subTree2, None)
                subTree1.iD = fn2(newickTree, i)
                subTree1.distance += fn(newickTree, i)
            elif reportUnaryNodes:
                subTree1 = ortheus.old.tree.BinaryTree(0.0, True, subTree1, None, None)
                subTree1.iD = fn2(newickTree, i)
                subTree1.distance += fn(newickTree, i)
            else:
                fn2(newickTree, i)
                subTree1.distance += fn(newickTree, i)
            return subTree1
        leafID = fn2(newickTree, i)
        return ortheus.old.tree.BinaryTree(fn(newickTree, i), False, None, None, leafID)
    return fn3(newickTree, [0])

def printBinaryTree(binaryTree, includeDistances, dontStopAtID=True):
    def fn(binaryTree):
        #print " tree Node ", binaryTree.left, binaryTree.right, binaryTree.distance, binaryTree.internal, binaryTree.iD 
        if binaryTree.iD != None:
            iD = str(binaryTree.iD)
        else:
            iD = ''
        if binaryTree.internal and (dontStopAtID or binaryTree.iD == None):
            if binaryTree.right != None:
                s = '(' + fn(binaryTree.left) + ',' + fn(binaryTree.right) + ')' + iD
            else:
                s = '(' + fn(binaryTree.left) + ')' + iD
        else:
            s = iD
        if includeDistances:
            return s + (":%f" % binaryTree.distance)
        return s
    return fn(binaryTree) + ';'

#########################################################
#########################################################
#########################################################
#functions for postion weight matrices
#########################################################
#########################################################
#########################################################

def pWMParser(pWMFile, alphabetSize=4):
    """
    reads in standard position weight matrix format,
    rows are different types of base, columns are individual residues
    """
    lines = open(pWMFile, 'r').readlines()
    assert len(lines) == alphabetSize
    l = [ [ float(i) ] for i in lines[0].split() ]
    for line in lines[1:]:
        l2 = [ float(i) for i in line.split() ]
        assert len(l) == len(l2)
        for i in xrange(0, len(l)):
            l[i].append(l2[i])
    for i in xrange(0, len(l)):
        j = sum(l[i]) + 0.0
        l[i] = [ k/j for k in l[i] ]
    return l

def writePWM(pWM, out, alphabetSize=4):
    """
    Writes file in standard PWM format, is reverse of pWMParser
    """
    for i in xrange(0, alphabetSize):
        out.write("%s\n" % ' '.join([ str(pWM[j][i]) for j in xrange(0, len(pWM)) ]))

def writePWMFile(pWM, outFile, alphabetSize=4):
    i = open(outFile, 'w')
    writePWM(pWM, i, alphabetSize)
    i.close()

def printPWM(pWM):
    """
    print a pWM, convenience function
    """
    writePWM(pWM, sys.stdout)

#########################################################
#########################################################
#########################################################
#exonerate functions
#########################################################
#########################################################
#########################################################

CIGAR_MATCH = 0
CIGAR_DELETE = 1
CIGAR_INSERT = 2
CIGAR_PLUS = 3
CIGAR_MINUS = 4
CIGAR_NA = 5

class Cigar:
    def __init__(self, queryID, queryStart, queryEnd, queryStrand, 
                 targetID, targetStart, targetEnd, targetStrand, score, oL):
        self.queryID = queryID
        self.queryStart = queryStart
        self.queryEnd = queryEnd
        self.queryStrand = queryStrand
        self.targetID = targetID
        self.targetStart = targetStart
        self.targetEnd = targetEnd
        self.targetStrand = targetStrand
        self.score = score
        self.oL = oL

def parseCigars(cigarFile):
    """
    converts file into list of cigars
    """
    cigarFile = open(cigarFile, 'r')
    cigars = []
    p = re.compile("cigar:\\s+(.+)\\s+([0-9]+)\\s+([0-9]+)\\s+([\\+\\-\\.])\\s+(.+)\\s+([0-9]+)\\s+([0-9]+)\\s+([\\+\\-\\.])\\s+([0-9]+)(\\s+(.*)\\s*)*")
    for line in cigarFile.readlines():
        i = p.match(line)
        if i != None:
            m = i.groups()
            if len(m) == 11:
                l = m[10].split(" ")
                ops = []
                assert len(l) % 2 == 0
                for j in xrange(0, len(l), 2):
                    if l[j] == 'M':
                        ops.append((CIGAR_MATCH, int(l[j+1])))
                    elif l[j] == 'D':
                        ops.append((CIGAR_DELETE, int(l[j+1]))) #a gap in the query
                    else: 
                        assert l[j] == 'I'
                        ops.append((CIGAR_INSERT, int(l[j+1]))) #a gap in the target
            else:
                ops = []
            def getStrand(i):
                if i == '+':
                    return CIGAR_PLUS
                if i == '-':
                    return CIGAR_MINUS
                assert i == '.'
                return CIGAR_NA
            cigars.append(Cigar(m[0], int(m[1]), int(m[2]), getStrand(m[3]), m[4], int(m[5]), int(m[6]), getStrand(m[7]), int(m[8]), ops))
    cigarFile.close()
    return cigars

RUN_EXONERATE_DNA=" --querytype dna --targettype dna "
RUN_EXONERATE_SOFTMASKED=" --softmaskquery true --softmasktarget true "
RUN_EXONERATE_AFFINE_LOCAL=" --model affine:local "

RUN_EXONERATE_DNA_SOFTMASKED_AFFINE_LOCAL=RUN_EXONERATE_DNA + RUN_EXONERATE_SOFTMASKED + RUN_EXONERATE_AFFINE_LOCAL
RUN_EXONERATE_DNA_SOFTMASKED_AFFINE_GLOBAL=RUN_EXONERATE_DNA_SOFTMASKED_AFFINE_LOCAL + " --bestn 1 "

def runExonerate(querySeqFile, targetSeqFile, path='exonerate', stringOptions=RUN_EXONERATE_DNA_SOFTMASKED_AFFINE_LOCAL):
    """
    runs exonerate, with given options, returns a list of cigars
    """
    f = getTempFile("bioio")
    com = "%s --showcigar true --showvulgar false --showalignment false %s %s %s > %s" % (path, stringOptions, querySeqFile, targetSeqFile, f)
    print "commmand :", com
    exitValue = os.system(com)
    if exitValue != 0:
        logger.info("Something went wrong calling Exoncerate : %i ", exitValue)
        sys.exit(1)
    c = parseCigars(f)
    os.remove(f)
    return c
    
def main():
    pass

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()