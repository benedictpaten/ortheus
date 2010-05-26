#!/usr/bin/env python

import sys

########################################################
########################################################
#This is a script to compile a 3-leaf transducer produced using the
#TransducerComposer.py script to C code, as can then used by Ortheus.
#In addition to the 3-leaf transducer file produced by TransducerComposer.py
#it also takes a .param file which is a procedurally arranged file
#containing C expressions for the transducer model, as well as arguments
#to the substitution model.
########################################################
########################################################
  

########################################################
########################################################
#Data-structures
########################################################
########################################################

#State types, as used in the three leaf model.
START = "START"
END = "END"
SILENT = "SILENT"
INSERT_X = "INSERT_X"
INSERT_Z = "INSERT_Z"
DELETE_X = "DELETE_X"
DELETE_Z = "DELETE_Z"
DELETE_XZ = "DELETE_XZ"
MATCH = "MATCH"

stateTypes = set((START, END, SILENT, 
                  INSERT_X, INSERT_Z, 
                  DELETE_X, DELETE_Z,
                  DELETE_XZ, MATCH))

class Model:
    """Data-structure to describe the three leaf transducer.
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
    
    def isStart(self, state): 
        return self.getStateType(state) == START
    
    def sortStates(self):
        """Function sort states into an order so that they can
        be computed for DP. Bit dense, probably the wrong place for this
        function.
        """
        l = [0]*len(self.statesList)
        for i in xrange(len(self.statesList)+10):
            for stateTo in self.statesList:
                toType = self.getStateType(stateTo)
                toIndex = self.statesList.index(stateTo)
                for stateFrom in self.getTransitionsTo(stateTo):
                    fromType = self.getStateType(stateFrom)
                    fromIndex = self.statesList.index(stateFrom)
                    if l[fromIndex] >= l[toIndex]:
                        if toType == SILENT:
                            l[toIndex] = l[fromIndex]+1
                        elif toType == DELETE_XZ and fromType != SILENT and \
                                                     fromType != DELETE_XZ:
                            #print "better thab", stateFrom, stateTo
                            l[toIndex] = l[fromIndex]+1
        l = [ (l[i], self.statesList[i]) for i in xrange(len(l)) ]
        l.sort()
        self.statesList = [ i[1] for i in l ]
        
########################################################
########################################################
#Input IO Functions
########################################################
########################################################

def getNextNonCommentLine(file):
    line = file.readline()
    while line != '' and line[0] == '#':
        line = file.readline()
    return line
        
def parseParameterLine(line):
    assert line[0] == 'P' or line[0] == 'C' or line[0] == 'O' or line[0] == 'B'
    i = line.split()
    parameter = i[1]
    assert parameter != 'DX', "This is the reserved distance for the X-branch parameter, no need to specify it"
    assert parameter != 'DY', "This is the reserved distance for the X-branch parameter, no need to specify it"
    assert i[2] == '='
    return parameter, "".join(i[3:])

def parseState(line):
    assert line[0] == 'S'
    i = line.split()
    assert len(i) == 4
    state = i[1]
    assert i[2] == '='
    assert i[3] in stateTypes
    return state, i[3]

def parseTransition(line):
    assert line[0] == 'T'
    i = line.split()
    stateFrom = i[1]
    assert i[2] == '-->'
    stateTo = i[3]
    assert i[4] == '='
    return stateFrom, stateTo, "".join(i[5:])

def parseFile(inputFile):
    """This function reads the input file using the accessory functions above.
    """
    primaryParameterList = []
    parameterList = []
    cParameters = {}
    forwardModel = Model()
    silentModel = Model()
    subModel = None
        
    inputFile = open(inputFile, 'r')
    line = getNextNonCommentLine(inputFile)
    while line != '':
        if line[0] == 'P' or line[0] == 'C' or line[0] == 'O' or line[0] == 'B':
            name, value = parseParameterLine(line)
            parameterList.append((name, value))
            if line[0] == 'O' or line[0] == 'B':
                primaryParameterList.append((name, value))
            if line[0] == 'C' or line[0] == 'B':
                cParameters[name] = value
        elif line[0] == 'S':
            state, type = parseState(line)
            forwardModel.addState(state, type)
            if type == SILENT or type == DELETE_XZ:
                silentModel.addState(state, type)
        elif line[0] == 'T':
            fromState, toState, value  = parseTransition(line)
            if fromState in silentModel.getStateNames() and \
            toState in silentModel.getStateNames():
                silentModel.addTransition(fromState, toState, value)
            forwardModel.addTransition(fromState, toState, value)
        elif line[0] == 'M':
            i = line.split()[1:]
            #Write out sub-model
            if i[0] == 'HKY':
                subModel = HKYSubModel(i[1:])
            else:
                raise RuntimeError("Don't recognise substitution model")
        line = getNextNonCommentLine(inputFile)
        
    #Clean up
    inputFile.close()
    assert len(forwardModel.getStateNames(START)) == 1
    assert len(forwardModel.getStateNames(END)) == 1
    forwardModel.sortStates()
    silentModel.sortStates()
    return primaryParameterList, parameterList, cParameters, forwardModel, silentModel, subModel    

def replaceCParameters(string, cParameters):
    #Computes expression, returns float
    i = [ (len(i), i) for i in cParameters.keys() ]
    i.sort()
    i.reverse()
    for length, parameter in i:
        string = string.replace(parameter, "1.0")
        #string.replace(parameter, "1.0")
    return string

########################################################
########################################################
#Output IO Functions
########################################################
########################################################

class HKYSubModel:
    """This class is used to contain some logic to write out a subsitution
    model for Ortheus.
    """
    def __init__(self, params):
        self.EXPECTED_FREQUENCY_A = float(params[0])
        self.EXPECTED_FREQUENCY_C = float(params[1])
        self.EXPECTED_FREQUENCY_G = float(params[2])
        self.EXPECTED_FREQUENCY_T = float(params[3])
        self.TRANSITION_TRANSVERSION_RATIO = float(params[4])
    
    def writeCFile(self, outputFile):
        outputFile.write("\ttemp->subModelX = constructHKYSubModel(DX,\n")
        outputFile.write("\t\t%f, %f,\n" % (self.EXPECTED_FREQUENCY_A, self.EXPECTED_FREQUENCY_C))
        outputFile.write("\t\t%f, %f,\n" % (self.EXPECTED_FREQUENCY_G, self.EXPECTED_FREQUENCY_T))
        outputFile.write("\t\t%f);\n\n" % (self.TRANSITION_TRANSVERSION_RATIO))
        outputFile.write("\ttemp->subModelY = constructHKYSubModel(DZ,\n")
        outputFile.write("\t\t%f, %f,\n" % (self.EXPECTED_FREQUENCY_A, self.EXPECTED_FREQUENCY_C))
        outputFile.write("\t\t%f, %f,\n" % (self.EXPECTED_FREQUENCY_G, self.EXPECTED_FREQUENCY_T))
        outputFile.write("\t\t%f);\n\n" % (self.TRANSITION_TRANSVERSION_RATIO))
        outputFile.write("\ttemp->ancestorProbs = st_malloc(sizeof(float)*ALPHABET_SIZE);\n")
        outputFile.write("\ttemp->ancestorProbs[0] = %f;\n" % self.EXPECTED_FREQUENCY_A) 
        outputFile.write("\ttemp->ancestorProbs[1] = %f;\n" % self.EXPECTED_FREQUENCY_C) 
        outputFile.write("\ttemp->ancestorProbs[2] = %f;\n" % self.EXPECTED_FREQUENCY_G) 
        outputFile.write("\ttemp->ancestorProbs[3] = %f;\n\n" % self.EXPECTED_FREQUENCY_T) 

def writeHModel(outputFile, primaryParameterList, parameterList, cParameters, forwardModel, silentModel, subModel):
    """This function writes the header file.
    """
    outputFile = open(outputFile, 'w')
    
    #Write header lines
    outputFile.write("#ifndef XYZMODELC_H_\n")
    outputFile.write("#define XYZMODELC_H_\n\n")
    
    outputFile.write('#include <inttypes.h>\n')
    outputFile.write('#include "sonLib.h"\n')
    
     
    #Number of states
    outputFile.write("struct CombinedTransitionModel {\n")
    
    outputFile.write("\tstruct SubModel *subModelX;\n")
    outputFile.write("\tstruct SubModel *subModelY;\n")
    outputFile.write("\tfloat *ancestorProbs;\n")
    
    #Include root parameters
    outputFile.write("\tint32_t includeRoot;\n")
    
    #Now write out all transitions
    for toState in forwardModel.getStateNames():
        i = forwardModel.getTransitionsTo(toState)
        for fromState in i.keys():
            #print "This is", fromState, toState, i[fromState]
            outputFile.write("\tfloat ft%s_%s;\n" % (fromState, toState))
            outputFile.write("\tfloat tbt%s_%s;\n" % (fromState, toState))
    outputFile.write("\n")
    
    for toState in silentModel.getStateNames():
        for fromState in silentModel.getStateNames():
            #print "This is", fromState, toState, i[fromState]
            outputFile.write("\tfloat lt%s_%s;\n" % (fromState, toState))
    
    outputFile.write("};\n")
    
    outputFile.write("struct ParameterStruct {\n")
    for parameterName, value in primaryParameterList:
        outputFile.write("\tfloat %s;\n" % parameterName)
    outputFile.write("};\n")
    
    outputFile.write("#endif /*XYZMODELC_H_*/\n")
    outputFile.close()

def writeCModel(outputFile, primaryParameterList, parameterList, cParameters, forwardModel, silentModel, subModel):
    """This function writes the .c file.
    """
    outputFile = open(outputFile, 'w')
    
    #Write header lines
    outputFile.write("#include <stdio.h>\n")
    outputFile.write("#include <assert.h>\n")
    outputFile.write('#include <inttypes.h>\n')
    outputFile.write('#include "sonLib.h"\n')
    outputFile.write('#include "substitutionC.h"\n')
    outputFile.write('#include "sequenceGraphC.h"\n')
    outputFile.write('#include "xyzModelC.h"\n\n')
    
    #Number of states
    outputFile.write("#define STATE_NO %i\n\n" % (len(forwardModel.getStateNames())-2))
    
    #Now write out states.
    stateIndex = 0
    for state in forwardModel.getStateNames():
        i = forwardModel.getStateType(state) 
        if i != START and i != END:
            outputFile.write("const int32_t %s = %i;\n" %(state, stateIndex))
            stateIndex += 1  
    
    outputFile.write("\n")

    #Now write state check functions.
    l = [ (SILENT, "Silent"),  (DELETE_XZ, "XYDelete"),\
          (INSERT_X, "XInsert"), (INSERT_Z, "YInsert"),\
          (DELETE_X, "XDelete"), (DELETE_Z, "YDelete"),\
          (MATCH, "Match") ]
    for stateType, functionName in l:
        outputFile.write("inline int32_t is%s(int32_t state) {\n" % functionName)
        for state in forwardModel.getStateNames(stateType):
            outputFile.write("\tif(state == %s) return TRUE;\n" % state)
        outputFile.write("\treturn FALSE; \n}\n\n")
    
    outputFile.write("struct CombinedTransitionModel *constructCombinedTransitionModel(float DX, float DZ, int32_t includeRoot, struct ParameterStruct *pS) {\n") 
    outputFile.write("\tstruct CombinedTransitionModel *temp = st_malloc(sizeof(struct CombinedTransitionModel));\n\n")
    
    outputFile.write('\tst_logInfo("Building combined transition model, DX: %%f, DY %%f\\n", DX, DZ);\n')
    
    #Are we including the root in these probs.
    outputFile.write("\ttemp->includeRoot = includeRoot;\n")
    outputFile.write('\tst_logInfo("Is root? " INT_STRING " \\n", temp->includeRoot);\n')
    
    #Write out submodel
    subModel.writeCFile(outputFile)
    
    #Write out C parameters.
    outputFile.write('\tint32_t i;\n\tint32_t j;\n')
    
    outputFile.write('\tfor(i=0; i<ALPHABET_SIZE; i++) {\n')
    outputFile.write('\t\tst_logInfo("Stationary frequency %i value %%f \\n", i, temp->subModelX->stationaryDistribution[i]);\n\t}\n')
    
    outputFile.write('\tfor(i=0; i<ALPHABET_SIZE; i++) {\n')
    outputFile.write('\t\tfor(j=0; j<ALPHABET_SIZE; j++) {\n')
    outputFile.write('\t\t\tst_logInfo("Substitition probability (X Branch) %i %i, forward value %%f, backward value %%f\\n", i, j, temp->subModelX->forward[i*ALPHABET_SIZE + j], temp->subModelX->backward[i*ALPHABET_SIZE + j]);\n\t\t}\n\t}\n')
    
    outputFile.write('\tfor(i=0; i<ALPHABET_SIZE; i++) {\n')
    outputFile.write('\t\tfor(j=0; j<ALPHABET_SIZE; j++) {\n')
    outputFile.write('\t\t\tst_logInfo("Substitition probability (Z Branch) %i %i, forward value %%f, backward value %%f\\n", i, j, temp->subModelY->forward[i*ALPHABET_SIZE + j], temp->subModelY->backward[i*ALPHABET_SIZE + j]);\n\t\t}\n\t}\n')
     
    for parameterName, value in primaryParameterList:
        outputFile.write("\tfloat %s = pS->%s;\n" % (parameterName, parameterName))
        outputFile.write('\tst_logInfo("Parameter %s, value %%f\\n", %s);\n' % (parameterName, parameterName))
     
    #Now write out parameters.
    for parameterName, value in parameterList:
        if (parameterName, value) not in primaryParameterList:
            outputFile.write("\tfloat %s = %s;\n" % (parameterName, value))
            outputFile.write('\tst_logInfo("Parameter %s, value %%f\\n", %s);\n' % (parameterName, parameterName))
        
    
    #Now write out all transitions
    for toState in forwardModel.getStateNames():
        i = forwardModel.getTransitionsTo(toState)
        for fromState in i.keys():
            #print "This is", fromState, toState, i[fromState]
            outputFile.write("\ttemp->ft%s_%s = LOG(%s);\n" % (fromState, toState, i[fromState]))
            outputFile.write("\ttemp->tbt%s_%s = LOG(%s);\n" % (fromState, toState, replaceCParameters(i[fromState], cParameters)))
            outputFile.write('\tst_logInfo("From state %s, To state %s, Forward Parameter %%f Traceback parameter %%f \\n", temp->ft%s_%s, temp->tbt%s_%s);\n' % (fromState, toState, fromState, toState, fromState, toState))
    outputFile.write("\n")
    
    #States already organised into correct order.
    #Write out non-log parameters.
    for toState in silentModel.getStateNames():
        i = silentModel.getTransitionsTo(toState)
        for fromState in i.keys():
            #print "This is", fromState, toState, i[fromState]
            outputFile.write("\tfloat tlt%s_%s = %s;\n" % (fromState, toState, i[fromState]))
    outputFile.write("\t float fA[STATE_NO];\n") 
    outputFile.write("\t float fA1[STATE_NO];\n") 
    outputFile.write("\t float fA2[STATE_NO];\n") 
    outputFile.write("\tfloat *l1;\n\tfloat *l2;\n\tfloat *l3;\n") 
    outputFile.write("\tl1 = fA1;\n")
    outputFile.write("\tl2 = fA2;\n")
    for fromState in silentModel.getStateNames():
        #Initialise level 0 and 1 and output-matrix
        outputFile.write("\tfor(i=0; i<STATE_NO; i++) {\n\t\tl1[i] = 0.0;\n\t\tl2[i] = 0.0;\n\t\tfA[i] = 0.0;\n\t}\n")
        i = silentModel.getStateType(fromState)
        if i == SILENT:
            outputFile.write("\tl1[%s] = 1.0;\n" % fromState)
        else:
            outputFile.write("\tl2[%s] = 1.0;\n" % fromState)
        #Do levels
        outputFile.write("\tfor(i=0; i<1000; i++) {\n")
        #Do silent transitions
        for toState in silentModel.getStateNames(SILENT):
            for fromState2 in silentModel.getTransitionsTo(toState):
                outputFile.write("\t\tl1[%s] += l1[%s] * tlt%s_%s;\n" % (toState, fromState2, fromState2, toState))
        #Add level to output
        for state in silentModel.getStateNames():
            outputFile.write("\t\tfA[%s] += l1[%s];\n" % (state, state))
        #Do delete-XZ transitions
        for toState in silentModel.getStateNames(DELETE_XZ):
            for fromState2 in silentModel.getTransitionsTo(toState):
                outputFile.write("\t\tl2[%s] += l1[%s] * tlt%s_%s;\n" % (toState, fromState2, fromState2, toState))
        #Make level 2 level 1 and initialise level 2 with 0.0s
        outputFile.write("\t\tl3 = l1;\n")
        outputFile.write("\t\tl1 = l2;\n")
        outputFile.write("\t\tl2 = l3;\n")
        outputFile.write("\t\tfor(j=0; j<STATE_NO; j++)\n\t\t\tl2[j] = 0.0;\n\t}\n")
        #Now write out these computed transitions
        for toState in silentModel.getStateNames():
            outputFile.write("\ttemp->lt%s_%s = LOG(fA[%s]);\n" % (fromState, toState, toState))
            outputFile.write('\tst_logInfo("From state %s, To state %s, LOOP PARAMETER %%f \\n", temp->lt%s_%s);\n' % (fromState, toState, fromState, toState))
    outputFile.write("\treturn temp; \n}\n\n")
    
    l = [ (INSERT_X, "insertX"), (INSERT_Z, "insertY"), (DELETE_X, "deleteX"), (DELETE_Z, "deleteY") ]
    for stateType, functionName in l:
        #Forward function
        outputFile.write("inline void %sFn(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, struct Edge *edge," % functionName) 
        outputFile.write("void (*assignFn)(struct AlignmentDataStructures *aDS, int32_t, int32_t, float)) {\n")
        for toState in forwardModel.getStateNames(stateType):
            for fromState in forwardModel.getTransitionsTo(toState).keys():
                if not forwardModel.isStart(fromState):
                    outputFile.write("\tassignFn(aDS, %s, %s, model->ft%s_%s + edge->edgeScore + edge->subScore);\n" \
                                     % (fromState, toState, fromState, toState))
        outputFile.write("}\n\n")
        #Traceback function
        outputFile.write("inline void %sFn_TraceBack(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, struct Edge *edge," % functionName) 
        outputFile.write("void (*assignFn)(struct AlignmentDataStructures *aDS, int32_t, int32_t, float, float)) {\n")
        
        outputFile.write("\tif (model->includeRoot) {\n\t\tfloat i;\n")
        for toState in forwardModel.getStateNames(stateType):
            for fromState in forwardModel.getTransitionsTo(toState).keys():
                if not forwardModel.isStart(fromState):
                    outputFile.write("\t\ti = model->ft%s_%s + edge->edgeScore + edge->subScore;\n" % (fromState, toState))
                    outputFile.write("\t\tassignFn(aDS, %s, %s, i, i);\n" % (fromState, toState))
        outputFile.write("\t}\n")  
        outputFile.write("\telse {\n")
        for toState in forwardModel.getStateNames(stateType):
            for fromState in forwardModel.getTransitionsTo(toState).keys():
                if not forwardModel.isStart(fromState):
                    outputFile.write("\t\tassignFn(aDS, %s, %s, model->ft%s_%s + edge->edgeScore + edge->subScore, model->tbt%s_%s + edge->edgeScore);\n" \
                                     % (fromState, toState, fromState, toState, fromState, toState))
        outputFile.write("\t}\n")
                    
        outputFile.write("}\n\n")
    
    #Forward match function
    outputFile.write("inline void matchFn(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model,") 
    outputFile.write("struct Edge *edgeX, struct Edge *edgeY, void (*assignFn)(struct AlignmentDataStructures *aDS, int32_t, int32_t, float)) {\n")
    outputFile.write("\t float m[ALPHABET_SIZE];\n")
    outputFile.write("\t float i;\n")
    outputFile.write("\t float j;\n")
    outputFile.write("\tmultiplyWV(edgeX->wV, edgeY->wV, m, ALPHABET_SIZE);\n")
    outputFile.write("\ti = LOG(combineWV(m, model->ancestorProbs, ALPHABET_SIZE));\n")
    for toState in forwardModel.getStateNames(MATCH):
        for fromState in forwardModel.getTransitionsTo(toState).keys():
            if not forwardModel.isStart(fromState):
                outputFile.write("\tj = model->ft%s_%s + edgeX->edgeScore + edgeY->edgeScore + i;\n" % (fromState, toState))
                outputFile.write("\tassignFn(aDS, %s, %s, j);\n"  % (fromState, toState))
    outputFile.write("}\n\n")
    
    #Traceback match function
    outputFile.write("inline void matchFn_TraceBack(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model,") 
    outputFile.write("struct Edge *edgeX, struct Edge *edgeY, void (*assignFn)(struct AlignmentDataStructures *aDS, int32_t, int32_t, float, float)) {\n")
    outputFile.write("\t float m[ALPHABET_SIZE];\n")
    outputFile.write("\t float i;\n")
    outputFile.write("\t float j;\n")
    outputFile.write("\t float k;\n")
    outputFile.write("\t float l;\n")
    outputFile.write("\tmultiplyWV(edgeX->wV, edgeY->wV, m, ALPHABET_SIZE);\n")
    outputFile.write("\ti = LOG(sumWV(m, ALPHABET_SIZE));\n")
    outputFile.write("\tj = LOG(combineWV(m, model->ancestorProbs, ALPHABET_SIZE));\n")
    
    outputFile.write("\tif (model->includeRoot) {\n")
    for toState in forwardModel.getStateNames(MATCH): 
        for fromState in forwardModel.getTransitionsTo(toState).keys():
            if not forwardModel.isStart(fromState):
                outputFile.write("\t\tl = model->ft%s_%s + edgeX->edgeScore + edgeY->edgeScore + j;\n" % (fromState, toState))
                outputFile.write("\t\tassignFn(aDS, %s, %s, l, l);\n"  % (fromState, toState))
    outputFile.write("\t}\n") 
    outputFile.write("\telse {\n")
    for toState in forwardModel.getStateNames(MATCH): 
        for fromState in forwardModel.getTransitionsTo(toState).keys():
            if not forwardModel.isStart(fromState):
                outputFile.write("\t\tk = model->tbt%s_%s + edgeX->edgeScore + edgeY->edgeScore + i;\n" % (fromState, toState))
                outputFile.write("\t\tl = model->ft%s_%s + edgeX->edgeScore + edgeY->edgeScore + j;\n" % (fromState, toState))
                outputFile.write("\t\tassignFn(aDS, %s, %s, l, k);\n"  % (fromState, toState))
    outputFile.write("\t}\n")
    outputFile.write("}\n\n")
    
    #Forward silent/delete_xz function
    outputFile.write("inline void silentFn(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, float *cell,")
    outputFile.write("void (*assignFn)(struct AlignmentDataStructures *aDS, int32_t, int32_t, float)) {\n")
    
    for toState in silentModel.getStateNames():
        for fromState in forwardModel.getTransitionsTo(toState).keys():
            if not forwardModel.isStart(fromState):
                i = forwardModel.getStateType(fromState)
                if i != SILENT and i != DELETE_XZ:
                    outputFile.write("\tassignFn(aDS, %s, %s, model->ft%s_%s);\n"  % (fromState, toState, fromState, toState))
    
    outputFile.write("\t float fA[STATE_NO];\n\tfloat i;\n")
    for state in silentModel.getStateNames():
        outputFile.write("\tfA[%s] = cell[%s];\n" % (state, state))
    for toState in silentModel.getStateNames():
        outputFile.write("\ti = LOG_ZERO;\n") 
        for fromState in silentModel.getStateNames():
            outputFile.write("\tLOG_PLUS_EQUALS(&i, fA[%s] + model->lt%s_%s);\n" % (fromState, fromState, toState))
        outputFile.write("\tcell[%s] = i;\n" % toState)
    outputFile.write("}\n\n") 
    
    #Traceback silent/delete_xz functions
    l = [ (SILENT, "silent"), (DELETE_XZ, "delete") ]
    for stateType, functionName in l:
        outputFile.write("inline void %sFn_TraceBack(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, void (*assignFn)(struct AlignmentDataStructures *aDS, int32_t, int32_t, float, float)) {\n" % functionName)
        outputFile.write("\tif (model->includeRoot) {\n")
        for toState in forwardModel.getStateNames(stateType):
            for fromState in forwardModel.getTransitionsTo(toState).keys():
                if not forwardModel.isStart(fromState):
                    outputFile.write("\t\tassignFn(aDS, %s, %s, model->ft%s_%s, model->ft%s_%s);\n" \
                                     % (fromState, toState, fromState, toState, fromState, toState))
        outputFile.write("\t}\n") 
        outputFile.write("\telse {\n")
        for toState in forwardModel.getStateNames(stateType):
            for fromState in forwardModel.getTransitionsTo(toState).keys():
                if not forwardModel.isStart(fromState):
                    outputFile.write("\t\tassignFn(aDS, %s, %s, model->ft%s_%s, model->tbt%s_%s);\n" \
                                     % (fromState, toState, fromState, toState, fromState, toState))        
        outputFile.write("\t}\n}\n\n")
    
    #Destructor
    outputFile.write("void destructCombinedTransitionModel(struct CombinedTransitionModel *model) { free(model); }\n\n")
    
    #State no
    outputFile.write("int32_t stateNo() { return STATE_NO; }\n\n")
        
    #Start states
    outputFile.write("float *startStates(struct CombinedTransitionModel *model) {\n")
    outputFile.write("\tfloat *i; int32_t j; i = st_malloc(sizeof(float)*STATE_NO);\n")
    outputFile.write("\tfor (j=0;j<STATE_NO; j++) i[j] = LOG_ZERO;\n")
    i = forwardModel.getStateNames(START)
    assert len(i) == 1, "More than one start state"
    startState = list(i)[0]
    for toState in forwardModel.getTransitionsFrom(startState).keys():
        if forwardModel.getStateType(toState) != END:
            outputFile.write("\ti[%s] = model->ft%s_%s;\n" % (toState, startState, toState))
    outputFile.write("\treturn i; \n}\n\n")
    
    #End states
    outputFile.write("float *endStates(struct CombinedTransitionModel *model) {\n")
    outputFile.write("\tfloat *i; int32_t j; i = st_malloc(sizeof(float)*STATE_NO);\n")
    outputFile.write("\tfor (j=0;j<STATE_NO; j++) i[j] = LOG_ZERO;\n\n")
    i = forwardModel.getStateNames(END)
    assert len(i) == 1, "More than one end state"
    endState = list(i)[0]
    for fromState in forwardModel.getTransitionsTo(endState).keys():
        if forwardModel.getStateType(fromState) != START:
            outputFile.write("\ti[%s] = model->ft%s_%s;\n" % (fromState, fromState, endState))
    outputFile.write("\treturn i; \n}\n\n")
    
    #Paramaeter struct constructor
    outputFile.write("struct ParameterStruct *constructParamStruct(int argc, char *argv[]) {\n")
    outputFile.write("\tint32_t i;\n\tchar *mod;\n\tstruct ParameterStruct *pM = st_malloc(sizeof(struct ParameterStruct));\n")
    
    outputFile.write("\tfloat floatParser;\n")
    for parameterName, value in primaryParameterList:
        outputFile.write("\tpM->%s = %s;\n" % (parameterName, value))
    outputFile.write("\tfor(i=1; i<argc; i++) {\n")
    outputFile.write("\t\tmod = argv[i];\n")
    outputFile.write("\t\tif(mod[0] != '-') {\n\t\t\tcontinue;\n\t\t}\n\t\tassert(mod[0] == '-');\n")
    outputFile.write("\t\tswitch(mod[1]) {\n")
    characterNo = 65
    for parameterName, value in primaryParameterList:
        outputFile.write("\t\t\tcase '%s':\n" % chr(characterNo))
        characterNo += 1
        outputFile.write('\t\t\t\tsscanf(argv[++i], "%f", &floatParser);\n')
        outputFile.write("\t\t\t\tpM->%s = floatParser;\n\t\t\tbreak;\n" % parameterName)
    outputFile.write("\t\t}\n\t}\n")
    outputFile.write("\treturn pM; \n}\n\n")
    
    #Parameter struct destructor.
    outputFile.write("void destructParamStruct(struct ParameterStruct *pM) {\n\tfree(pM);\n}\n\n")
    
    #Parameter reporter
    outputFile.write("void printParamStruct() {\n")
    characterNo = 65
    for parameterName, value in primaryParameterList:
        outputFile.write('\tfprintf(stderr, "\t-%s [FLOAT] value of parameter %s , default: %s \\n");\n' % (chr(characterNo), parameterName, value))
        characterNo += 1
    outputFile.write("}\n\n")
    
    outputFile.close()
    
########################################################
########################################################
#Control scripts.
########################################################
########################################################

def compileTransducerToCCode(inputFile, outputCFile, outputHFile):
    #First oarse input-file
    primaryParameterList, parameterList, cParameters, forwardModel, silentModel, subModel = parseFile(inputFile)
    #Writes c file
    writeCModel(outputCFile, primaryParameterList, parameterList, cParameters, forwardModel, silentModel, subModel)
    writeHModel(outputHFile, primaryParameterList, parameterList, cParameters, forwardModel, silentModel, subModel)

def main():
    compileTransducerToCCode(sys.argv[1], sys.argv[2], sys.argv[3])

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()
