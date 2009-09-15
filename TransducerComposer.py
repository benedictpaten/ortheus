#!/usr/bin/env python

import sys
import os
import re

########################################################
########################################################
#This is a script to compose a 3-leaf transducer  (i.e. an internal of a tree) 
#using two three branch transducers, one for the root (y) and one for each 
#descendant branch (x and z).
#It takes as input a file using a loose version of the Phylocomposer format 
#(lisp inspired)
#(see Holmes 2007 Bioinformatics paper) and affineModel.sexpr for example.
#It spits out a file in a special formal read only by TransducerCompiler.py.
########################################################
########################################################

def getNextNonCommentLine(file):
    line = file.readline()
    while line != '' and line[0] == '#':
        line = file.readline()
    return line

#####These are the state types of a transducer, they also include values only
#used in the three leaf model.
START = "START"
END = "END"
SILENT = "SILENT"
INSERT_X = "INSERT_X"
INSERT_Z = "INSERT_Z"
DELETE_X = "DELETE_X"
DELETE_Z = "DELETE_Z"
DELETE_XZ = "DELETE_XZ"
DELETE = "DELETE"
INSERT = "INSERT"
MATCH = "MATCH"

class Transducer:
    """These is a convenient object for holding the states and transitions
    of a transducer. Something similar is in TransducerCompiler.py, but
    I've basically cut and pasted for now.
    """
    def __init__(self):
        self.matrix = {}
        self.inverseMatrix = {}
        self.stateTypes = {}
        self.statesList = []
        
    def addState(self, name, type):
        assert name not in self.matrix
        self.matrix[name] = {}
        self.inverseMatrix[name] = {}
        self.stateTypes[name] = type
        self.statesList.append(name)
    
    def addTransition(self, stateNameFrom, stateNameTo, value):
        self.matrix[stateNameFrom][stateNameTo] = value
        self.inverseMatrix[stateNameTo][stateNameFrom] = value
        
    def getTransitionsFrom(self, state):
        return self.matrix[state].copy()
    
    def getTransitionsTo(self, state):
        return self.inverseMatrix[state].copy()
        
    def getStateNames(self, stateType=None):
        if stateType == None:
            return list(self.statesList)
        return [ i for i in self.statesList if self.getStateType(i) == stateType ]
    
    def getStateNumber(self):
        return len(self.getStateNames())
    
    def getStateType(self, state): 
        return self.stateTypes[state]
    
    def report(self):
        print "Reporting on transducer"
        print "State names", " ".join(self.getStateNames())
        for fS in self.getStateNames():
            i = self.getTransitionsFrom(fS)
            for tS in i.keys():
                print "From state %s, to state %s, value %s" % (fS, tS, i[tS])
                
########################################################
########################################################
#IO Functions
########################################################
########################################################

def parseStateLine(line):
    name = re.compile("\(name ([^)]*)\)").search(line).group(1)
    type = re.compile("\(type ([^)]*)\)").search(line).group(1)
    hash = { 'start':START, 'end':END, 'wait':SILENT, 
             'silent':SILENT, 'insert':INSERT, 'match':MATCH,
             'delete':DELETE }
    return name, hash[type]

def parseTransitionLine(line):
    fState = re.compile("\(from ([^)]*)\)").search(line).group(1)
    tState = re.compile("\(to ([^)]*)\)").search(line).group(1)
    s = re.compile("\(label ([^)]*)\)").search(line)
    value = []
    if s != None:    
        value = [ s.group(1) ]
    return fState, tState, value

def parseInputTransducerFile(inputFile):
    """This function reads the lisp like transducer file.
    The transducer must contain three models, in order the branch-X model,
    the branch-Z model and finally the root model.
    """
    inputFile = open(inputFile, 'r')
    def parseTransducer(inputFile):
        line = getNextNonCommentLine(inputFile)
        transducer = Transducer()
        while line != '':
            if "transducer" in line:
                line = getNextNonCommentLine(inputFile)
                while "transition" not in line:
                    if "state" in line:
                        state, type = parseStateLine(line)
                        transducer.addState(state, type)
                    line = getNextNonCommentLine(inputFile)
                while line[0] != ')':
                    if "transition" in line:
                        fState, tState, value = parseTransitionLine(line)
                        transducer.addTransition(fState, tState, value)
                    line = getNextNonCommentLine(inputFile)
                return transducer
            line = getNextNonCommentLine(inputFile)
        return None
    branchTransducerX = parseTransducer(inputFile)
    branchTransducerZ = parseTransducer(inputFile)
    rootTransducer = parseTransducer(inputFile)
    inputFile.close()
    return branchTransducerX, branchTransducerZ, rootTransducer

def writeModel(states, transitions, outputFile):
    """Writes out a file for TransducerCompiler.py to read.
    """
    def fn(stateName):
        return stateName.replace('/', '')
    outputFile = open(outputFile, 'w')
    outputFile.write("\n# States: %s Transitions: %s \n" % (len(states), len(transitions)))
    for state in states.keys():
        outputFile.write("S %s = %s\n" % (fn(state), states[state]))
    for fS, tS in transitions.keys():
        outputFile.write("T %s --> %s = %s\n" % (fn(fS), fn(tS), " * ".join(transitions[(fS, tS)])))
    outputFile.close()
    
def writeDotFile(states, transitions, outputFile):
    """Writes out a graph-viz file describing the graph.
    """
    stateList = list(states.keys())
    outputFile = open(outputFile, 'w')
    outputFile.write("graph G {\n")
    outputFile.write("overlap=false\n")
    for state in states:
        outputFile.write("node[width=0.3,height=0.3,shape=box,style=filled,color=red,fontsize=14];\n")
        outputFile.write('n%in [label="%s"];\n' % (stateList.index(state), state))
    for transition in transitions:
        i, j = transition
        outputFile.write("edge[color=green,len=0.6,weight=100,dir=forward];\n")
        outputFile.write('n%in -- n%in [label="%s"];\n' % (stateList.index(i), stateList.index(j), "*".join(transitions[transition])))
    outputFile.write("}\n")    
    outputFile.close() 

########################################################
########################################################
#The actual algorithms
########################################################
########################################################

def composeInflatedTransducer(branchTransducerX, branchTransducerZ, rootTransducer):
    """This builds an 'inflated' version of the combined three leaf
    state graph. Subsequently we collapse out some of the silent states to
    make the graph more compact. The arguments are the three transducers
    previously parsed.
    """
    #Composite states.
    states = {}
    #State types.
    stateTypes = {}
    #Transitions between states.
    transitions = {}
    
    def getStartState(transducer):
        i = None
        for state in transducer.getStateNames():
            if transducer.getStateType(state) == START:
                assert i == None
                i = state
        assert i != None
        return i
    
    #Get start states for each branche's transducer.
    sX = getStartState(branchTransducerX)
    sY = getStartState(rootTransducer)
    sZ = getStartState(branchTransducerZ)
    
    #Construct composite start state.
    def createCompositeState(mX, sX, mY, sY, mZ, sZ):
        return '%s_%s/%s_%s/%s_%s' % (sX, mX, sY, mY, sZ, mZ)
    
    s = createCompositeState(0, sX, 0, sY, 0, sZ)
    
    #Decompose states
    def decomposeState(s):
        i = s.split('/')
        assert len(i) == 3
        return [ j.split('_')[0] for j in i ]
    
    def addToStates(s, sType):
        if s in states:
            assert states[s] == sType
            return False
        states[s] = sType
        return True
    
    addToStates(s, START)
    
    #Create stack
    stack = []
    #Push on next states
    def pushOnStates(s):
        def fn(t, s):
            i = t.getStateType(s)
            if i in [ MATCH, INSERT, DELETE ]:
                return 1
            return 0
        sX, sY, sZ = decomposeState(s)
        for sX2 in branchTransducerX.getTransitionsFrom(sX):
            for sZ2 in branchTransducerZ.getTransitionsFrom(sZ):
                for sY2 in rootTransducer.getTransitionsFrom(sY):
                    s2 = createCompositeState(fn(branchTransducerX, sX2), sX2, fn(rootTransducer, sY2), sY2, fn(branchTransducerZ, sZ2), sZ2)
                    if (s, s2) not in transitions and (s, s2) not in stack:
                        stack.append((s, s2))
    
    pushOnStates(s)
    #while stack not empty:
    while len(stack) != 0:
        pS, s = stack.pop()
        
        #decompose states
        pSX, pSY, pSZ = decomposeState(pS)
        sX, sY, sZ = decomposeState(s)
        
        def addTransition(pS, s, cost):
            if (pS, s) in transitions:
                assert transitions[(pS, s)] == cost
            transitions[(pS, s)] = cost
    
        if not branchTransducerX.getStateType(sX) in [ DELETE, MATCH, END ]:
            tCost = branchTransducerX.getTransitionsFrom(pSX)[sX]
            if branchTransducerX.getStateType(sX) == INSERT:
                s2 = createCompositeState(1, sX, 0, pSY, 0, pSZ)
                isNewState = addToStates(s2, INSERT_X)
            else:
                s2 = createCompositeState(0, sX, 0, pSY, 0, pSZ)
                isNewState = addToStates(s2, SILENT)
            addTransition(pS, s2, tCost)
            if isNewState:
                #Push on new states
                pushOnStates(s2)
                        
        elif not branchTransducerZ.getStateType(sZ) in [ DELETE, MATCH, END ]:
            tCost = branchTransducerZ.getTransitionsFrom(pSZ)[sZ]
            if branchTransducerZ.getStateType(sZ) == INSERT:
                s2 = createCompositeState(0, pSX, 0, pSY, 1, sZ)
                isNewState = addToStates(s2, INSERT_Z)
            else:
                s2 = createCompositeState(0, pSX, 0, pSY, 0, sZ)
                isNewState = addToStates(s2, SILENT)
            addTransition(pS, s2, tCost)
            if isNewState:
                #Push on new states
                pushOnStates(s2)
                        
        elif rootTransducer.getStateType(sY) == SILENT:
            tCost = rootTransducer.getTransitionsFrom(pSY)[sY]
            s2 = createCompositeState(0, pSX, 0, sY, 0, pSZ)
            isNewState = addToStates(s2, SILENT)
            addTransition(pS, s2, tCost)
            if isNewState:
                pushOnStates(s2)
                
        elif rootTransducer.getStateType(sY) == INSERT and \
        branchTransducerX.getStateType(sX) in [ DELETE, MATCH ] and \
        branchTransducerZ.getStateType(sZ) in [ DELETE, MATCH ]:
            #find new state type
            if branchTransducerX.getStateType(sX) == MATCH:
                if branchTransducerZ.getStateType(sZ) == MATCH :
                    sType = MATCH
                else:
                    assert branchTransducerZ.getStateType(sZ) == DELETE
                    sType = DELETE_Z
            else:
                assert branchTransducerX.getStateType(sX) == DELETE
                if branchTransducerZ.getStateType(sZ) == MATCH:
                    sType = DELETE_X
                else:
                    assert branchTransducerZ.getStateType(sZ) == DELETE
                    sType = DELETE_XZ
            #find transition cost
            tCost = branchTransducerX.getTransitionsFrom(pSX)[sX]
            
            def appendCost(tCost1, tCost2):
                return tCost1 + tCost2
            
            tCost = appendCost(tCost, branchTransducerZ.getTransitionsFrom(pSZ)[sZ])
            tCost = appendCost(tCost, rootTransducer.getTransitionsFrom(pSY)[sY])
            #Add to states, determining if new
            isNewState = addToStates(s, sType)
            #Add transition
            addTransition(pS, s, tCost)
            #iterate through new states
            if True: #isNewState:
                pushOnStates(s)
        elif rootTransducer.getStateType(sY) == END and \
            branchTransducerX.getStateType(sX) == END and \
            branchTransducerZ.getStateType(sZ) == END:
            isNewState = addToStates(s, END)
            #Add transition
            addTransition(pS, s, rootTransducer.getTransitionsFrom(pSY)[sY])
        else:
            #These are states which are impossible to enter/pass through.
            pass

    return states, transitions

def compactModel(states, transitions, collapseCoefficient):
    """This takes the 'inflated' model and collapses out silent states
    to make the graph more compact.
    """
    #return states, transitions   
    for state in states.keys():
        if states[state] == SILENT:
            inTrans = [ i for i in transitions if i[1] == state ]
            outTrans = [ i for i in transitions if i[0] == state ]
            if len(inTrans) > 0 and len(outTrans) > 0 and \
            (len(inTrans)*len(outTrans) <= collapseCoefficient + len(inTrans)+len(outTrans)):
                #print "removing state", state
                for inT in inTrans:
                    for outT in outTrans:
                        i, j = inT
                        k, l = outT
                        transitions[(i, l)] = transitions[inT] + transitions[outT]
                for transition in inTrans + outTrans:
                    transitions.pop(transition)
                states.pop(state)
    return states, transitions

def checkModel(states, transitions):
    """Checks a graph is well formed.
    
    The collapse coefficient determines how 'collapsed' the resulting graph is.
    By default it is set to one, but if increased will increase the degree of 
    collapse.
    """
    assert len([ i for i in states if states[i] == START ]), "Just one start state please"
    assert len([ i for i in states if states[i] == END ]), "Just one end state please"
    for state in states:
        inTrans = [ i for i in transitions if i[1] == state ]
        outTrans = [ i for i in transitions if i[0] == state ]
        if len(inTrans) == 0:
            assert states[state] == START, "State is not start: %s" % state
        if len(outTrans) == 0:
            assert states[state] == END, "State is not end: %s" % state

########################################################
########################################################
#Control scripts.
########################################################
########################################################


def mainScript(phyloComposerInputFile, outputFile, outputDotFile, collapseCoefficient=1):
    branchTransducerX, branchTransducerZ, rootTransducer = parseInputTransducerFile(phyloComposerInputFile)
    states, transitions = composeInflatedTransducer(branchTransducerX, branchTransducerZ, rootTransducer)
    states, transitions = compactModel(states, transitions, collapseCoefficient)
    checkModel(states, transitions)
    writeModel(states, transitions, outputFile)
    writeDotFile(states, transitions, outputDotFile)

def main():
    phyloComposerInputFile =  sys.argv[1]
    outputFile = sys.argv[2]
    outputDotFile = sys.argv[3]
    collapseCoefficient = int(sys.argv[4])
    mainScript(phyloComposerInputFile, outputFile, 
               outputDotFile, collapseCoefficient)

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()
