/*
 * Copyright (C) 2008-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <ctype.h>

#include "fastCMaths.h"
#include "bioioC.h"
#include "commonC.h"
#include "substitutionC.h"
#include "commonC.h"
#include "hashTableC.h"
#include "heapC.h"


#include "sequenceGraphC.h"
#include "xyzModelC.h"
#include "constraintsC.h"

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//parameters
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

int64_t NUMBER_OF_SAMPLES = 100;
int64_t CONSTRAINT_HALF_ANTI_DIAGONAL_LOOSEN_SIZE = 10;
int64_t TOTAL_CONSTRAINTS = TRUE;
float DEFAULT_BINARY_TREE_DISTANCE = 0.0001f;
int64_t WRITE_MFA_ALIGNMENT = TRUE;

char *OUTPUT_ALIGNMENT_FILE = NULL;
int64_t VITERBI_COLUMN_GAP = 0;
char *TREE_STATES_INPUT_FILE = NULL;
char *TREE_STATES_OUTPUT_FILE = NULL;
char *SCORE_OUTPUT_FILE = NULL;

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//misc functions that are needed for the main script,
//but otherwise don't have a home
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

struct TreeNodeHolder {
    struct TreeNode *treeNode;
    struct Edge *edge;
};

struct TreeNodeHolder *constructTreeNodeHolder(struct TreeNode *treeNode, struct Edge *edge) {
    struct TreeNodeHolder *treeNodeHolder;

    treeNodeHolder = st_malloc(sizeof(struct TreeNodeHolder));
    treeNodeHolder->treeNode = treeNode;
    treeNodeHolder->edge = edge;
    return treeNodeHolder;
}

void destructTreeNodeHolder(struct TreeNodeHolder *treeNodeHolder) {
    //destructTreeNode(treeNodeHolder->treeNode);
    free(treeNodeHolder);
}

float *convertToWVSeq(char *seq, int64_t length, float *(*map)(char i)) {
    int64_t i;
    float *wV;
    float *seqWV;

    seqWV = st_malloc(sizeof(float)*length*ALPHABET_SIZE);
    for(i=0; i<length; i++) {
        wV = map(seq[i]);
        memcpy(seqWV + ALPHABET_SIZE*i, wV, ALPHABET_SIZE*sizeof(float));
    }
    return seqWV;
}

float *dNAMap_CharToWVFn(char i) {
    static float A[] = { 1.0f, 0.0f, 0.0f, 0.0f };
    static float C[] = { 0.0f, 1.0f, 0.0f, 0.0f };
    static float G[] = { 0.0f, 0.0f, 1.0f, 0.0f };
    static float T[] = { 0.0f, 0.0f, 0.0f, 1.0f };
    static float N[] = { 1.0f, 1.0f, 1.0f, 1.0f };

    switch(i) {
        case 'A':
        case 'a':
            return A;
        case 'C':
        case 'c':
            return C;
        case 'G':
        case 'g':
            return G;
        case 'T':
        case 't':
            return T;
        default:
            return N;
    }
}

char dNAMap_WVToACTGFn(float *wV) {
    int64_t i;
    float j;
    int64_t k;
    //converts char vectors to an ACTG
    if(wV[0] < -0.5f) { //is gap
        return '-';
        for(i=1; i<ALPHABET_SIZE; i++) {
            assert(wV[i] < 0);
        }
    }
    j = wV[0];
    assert(j >= 0);
    k = 0;
    for(i=1; i<ALPHABET_SIZE; i++) {
        assert(wV[i] >= 0);
        if(wV[i] > j) {
            j = wV[i];
            k = i;
        }
    }
    switch(k) {
        case 0:
                for(i=1; i<ALPHABET_SIZE; i++) {
                    if(wV[i] != j) {
                        return 'A';
                    }
                }
                return 'N';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        case 4: return 'N';
        default:
        assert(FALSE);
    }
    return -1;
}

void writeFastaAlignment(FILE *file, float **columnAlignment, int64_t columnNo, char **names, int64_t startSeq, int64_t seqInc, int64_t seqNo, char (*map)(float *), char *repeatBases) {
    //Writes out column alignment to given file multi-fasta format
    int64_t seq;
    int64_t i;
    int64_t j;
    char k;
    float *column;

    for(seq=startSeq; seq<seqNo; seq+=seqInc) {
        fprintf(file, ">%s\n", names[seq]);
        j = 0;
        for(i=0; i<columnNo; i++) {
            column = columnAlignment[i];
            k = map(column + seq*ALPHABET_SIZE);
            fprintf(file, "%c", repeatBases[j + seq] ? tolower(k) : k);
            j += seqNo;
        }
        fprintf(file, "\n");
    }
}

struct SequenceGraph *convertSeqToSeqGraph(float *seq, int64_t seqLength, struct TraversalID *traversalID) {
    struct List *edges;
    struct TreeNode *treeNode;
    float *wV;
    int64_t i;

    edges = constructEmptyList(seqLength, (void (*)(void *))destructEdge);
    //converts a list chars to a sequence graph
    for(i=0; i<seqLength; i++) {
        wV = seq + (i*ALPHABET_SIZE);
        treeNode = copyConstructTreeNode(TREE_NODE_LEAF, INT32_MAX, traversalID, NULL, NULL, wV);
        edges->list[i] = copyConstructEdge(i, i+1, LOG_ONE, LOG_ZERO, wV, FALSE, treeNode, i);
    }
    return constructSequenceGraph(edges, seqLength+1);
}

void lineariseTreeNodes(struct TreeNode *treeNode, int64_t addToList, struct Edge *edge, struct List *list) {
    if(treeNode->left != NULL) {
        lineariseTreeNodes(treeNode->left, TRUE, edge, list);
    }
    if(treeNode->treeNodeX != NULL) {
        lineariseTreeNodes(treeNode->treeNodeX, FALSE, edge, list);
    }
    if(treeNode->treeNodeY != NULL) {
        lineariseTreeNodes(treeNode->treeNodeY, FALSE, edge, list);
    }
    if(addToList) {
        listAppend(list, constructTreeNodeHolder(treeNode, edge));
    }
}

struct List *lineariseAlignment(struct List *alignment) {
    int64_t i;
    struct Edge *edge;
    struct List *list;
    //struct TreeNode *treeNode;

    //get in order non-normalised column probs for tree nodes
    list = constructEmptyList(alignment->length+1, (void (*)(void *))destructTreeNodeHolder);
    list->length = 0;
    //for edge in alignment:
    for(i=0; i<alignment->length; i++) {
        edge = alignment->list[i];
        lineariseTreeNodes(edge->treeNode, TRUE, edge, list);
    }
    return list;
}

struct TreeNode *getFirstNonSilentVertex(struct TreeNode *treeNode) {
    while(treeNode->type & TREE_NODE_EFFECTIVELY_SILENT) {
        treeNode = treeNode->treeNodeX;
        if(treeNode == NULL) {
            return NULL;
        }
    }
    return treeNode;
}

struct List *convertAlignmentToColumns(struct List *alignment, int64_t nodeNumber, struct SubModel **subModels, float *ancestorProbs) {
    struct List *list;
    int64_t i;
    struct TreeNodeHolder *treeNodeHolder;
    struct TreeNode *treeNode;

    list = constructEmptyList(alignment->length, free);
    list->length = 0;
    for(i=0; i<alignment->length; i++) {
        //uglyf(" hi %" PRIi64 " \n", i);
        treeNodeHolder = alignment->list[i];
        treeNode = getFirstNonSilentVertex(treeNodeHolder->treeNode);
        if(treeNode != NULL) {
            listAppend(list, felsensteins(treeNode, subModels, nodeNumber, ancestorProbs));
            /*INT_32 xx;
            FLOAT_32 *fA;
            fA = list->list[list->length-1];
            for(xx=0; xx<nodeNumber; xx++) {
                uglyf(" %c ", dNAMap_WVToACTGFn(fA + xx*ALPHABET_SIZE));
            }
            uglyf(" %f %" PRIi64 " %" PRIi64 " %" PRIi64 " \n ", treeNodeHolder->edge->edgeScore, treeNodeHolder->edge->iD, treeNodeHolder->treeNode->left, treeNodeHolder->treeNode->transitionID); //, ((struct TreeNode *)alignment->list[i])->);*/
        }
    }
    return list;
}

char *treeNodeNames(struct BinaryTree *binaryTree, char **labels, char **leafLabels) {
    char *i;
    char *j;
    char *k;

    if(binaryTree->internal) {
        i = treeNodeNames(binaryTree->left, labels, leafLabels);
        j = treeNodeNames(binaryTree->right, labels, leafLabels);
        k = st_malloc(sizeof(char) * (strlen(i) + strlen(j) + 2));
        return labels[binaryTree->traversalID->mid] = strcat(strcat(strcpy(k, i), "_"), j);
    }
    else {
        i = leafLabels[binaryTree->traversalID->leafNo];
        labels[binaryTree->traversalID->mid] = i;
        return i;
    }
}

struct Constraints ***buildPrimeConstraints(char *alignmentFile, int64_t seqNo, int64_t *seqLengths, struct BinaryTree *binaryTree) {
    struct Constraints ***primeConstraints;
    FILE *alignmentStream;
    struct List *alignmentSeqs;
    struct List *alignmentSeqLengths;
    struct List *fastaNames;

    //constraints from previously constructed alignment
    st_logInfo("Alignment file : %s\n", alignmentFile);
    if (alignmentFile != NULL) {
        alignmentSeqs = constructEmptyList(0, free);
        alignmentSeqLengths = constructEmptyList(0, free);
        fastaNames = constructEmptyList(0, free);

        alignmentStream = fopen(alignmentFile, "r");
        fastaRead(alignmentStream, alignmentSeqs, alignmentSeqLengths, fastaNames);
        fclose(alignmentStream);

        assert(seqNo == alignmentSeqs->length);
        primeConstraints = buildAllConstraints_FromAlignment((char **)alignmentSeqs->list, *((int64_t *)alignmentSeqLengths->list[0]), seqNo, seqLengths, CONSTRAINT_HALF_ANTI_DIAGONAL_LOOSEN_SIZE, '-', binaryTree, TOTAL_CONSTRAINTS);
        //free memory used
        destructList(alignmentSeqs);
        destructList(alignmentSeqLengths);
        destructList(fastaNames);
        //end
    }
    else {
        primeConstraints = buildAllConstraints_FromAlignment(NULL, INT32_MIN, seqNo, (int64_t *)seqLengths, CONSTRAINT_HALF_ANTI_DIAGONAL_LOOSEN_SIZE, '-', binaryTree, TOTAL_CONSTRAINTS);
    }
    st_logInfo("Finished parsing alignment file\n");
    st_logInfo("Built constraints lists\n");
    return primeConstraints;
}

struct List *parseSequences(char **argv, int64_t inputSequenceNumber, int64_t *seqNo, int64_t **seqLengthsA) {
    int64_t i;
    struct List *seqs;
    struct List *seqLengths;
    struct List *fastaNames;
    FILE *fastaStream;

    //the sequence files
    seqLengths = constructEmptyList(0, free);
    seqs = constructEmptyList(0, free);
    fastaNames = constructEmptyList(0, free);
    for(i=0; i<inputSequenceNumber; i++) {
        st_logInfo("Sequence file : %s\n", argv[i]);
        fastaStream = fopen(argv[i], "r");
        fastaRead(fastaStream, seqs, seqLengths, fastaNames);
        fclose(fastaStream);
    }
    *seqNo = seqs->length;
    *seqLengthsA = st_malloc(sizeof(int64_t)*(*seqNo));
    for(i=0; i<*seqNo; i++) {
        (*seqLengthsA)[i] = *((int64_t *)seqLengths->list[i]);
    }
    //free memory
    destructList(seqLengths);
    destructList(fastaNames);
    return seqs;
}

struct SequenceGraph **buildLeafSequenceGraphs(struct List *seqs, int64_t seqNo, int64_t *seqLengths, int64_t nodeNumber, struct BinaryTree **binaryTreeNodes) {
    int64_t i;
    int64_t j;
    float *wV;
    struct SequenceGraph **inputSequenceGraphs;

    inputSequenceGraphs = st_malloc(sizeof(void *)*seqNo);
    i=0;
    for(j=0; j<nodeNumber; j++) {
        if(!binaryTreeNodes[j]->internal) {
            wV = convertToWVSeq((char *)seqs->list[i], seqLengths[i], dNAMap_CharToWVFn);
            inputSequenceGraphs[i] = convertSeqToSeqGraph(wV, seqLengths[i], binaryTreeNodes[j]->traversalID);
            free(wV);
            i++;
        }
    }
    //no current memory clean up
    return inputSequenceGraphs;
}

//calculateRepeatBases((char **)seqs->list, (FLOAT_32 **)wvAlignment->list, seqNo, wvAlignment->length);

char *calculateRepeatBases(char **leafSeqs, float **columnAlignment, int64_t seqNo, int64_t alignmentLength) {
    int64_t i;
    int64_t j;
    int64_t k;
    int64_t nodeNo = seqNo * 2 - 1;
    int64_t repeatCount;
    int64_t totalCount;
    int64_t *indices;
    char *repeatBases;
    float *column;

    i = sizeof(int64_t)*seqNo;
    indices = st_malloc(i);
    memset(indices, 0, i);
    i = sizeof(char)*nodeNo*alignmentLength;
    repeatBases = st_malloc(i);
    memset(repeatBases, 0, i);

    k=0;
    for(i=0; i<alignmentLength; i++) {
        column = columnAlignment[i];
        repeatCount = 0;
        totalCount = 0;
        for(j=0; j<nodeNo; j+=2) {
            if(column[j*ALPHABET_SIZE] >= 0.0) {
                totalCount++;
                //is not a gap
                if(islower(leafSeqs[j/2][indices[j/2]++])) {
                    repeatCount++;
                    repeatBases[k + j] = 1;
                }
            }
        }
        if(repeatCount >= totalCount/2) {
            for(j=1; j<nodeNo; j+=2) {
                repeatBases[k + j] = 1;
            }
        }
        k += nodeNo;
    }
    //free memory
    free(indices);
    return repeatBases;
}

/*INT_32 expectedLength(struct BinaryTree *binaryTree, INT_32 *seqLengths) {
    INT_32 i, j;

    if(binaryTree->internal) {
        i = expectedLength(binaryTree->left, seqLengths);
        j = expectedLength(binaryTree->right, seqLengths);
        return i > j ? i : j;
    }
    return seqLengths[binaryTree->traversalID->leafNo];
}*/

int64_t lastInternalNodeState_TreeNode(struct TreeNode *treeNode, int64_t i, int64_t *state) {
    int64_t k;

    k = treeNode->traversalID->mid;
    if(k == i && treeNode->type != TREE_NODE_DELETE)  {
        //found it
        *state = treeNode->transitionID;
        return TRUE;
    }
    if(treeNode->treeNodeX != NULL) {
        if(lastInternalNodeState_TreeNode(treeNode->treeNodeX, i, state)) {
            return TRUE;
        }
    }
    if(treeNode->treeNodeY != NULL) {
        if(lastInternalNodeState_TreeNode(treeNode->treeNodeY, i, state)) {
            return TRUE;
        }
    }
    return FALSE;
}

int64_t lastInternalNodeState(struct List *alignment, int64_t i) {
    int64_t j;
    int64_t k;
    struct TreeNodeHolder *treeNodeHolder;
    struct TreeNode *treeNode;

    k = INT_MAX;
    for(j=alignment->length-1; j>= 0; j--) {
        treeNodeHolder = alignment->list[j];
        treeNode = treeNodeHolder->treeNode;
        if(lastInternalNodeState_TreeNode(treeNode, i, &k)) {
            break;
        }
    }
    return k;
}

void calculateTreeStates(struct List *alignment, int64_t nodeNo, int64_t *states) {
    //returns an array containing the end states of the different INT_32ernal node
    //alignments
    int64_t i;
    int64_t j;

    for(i=1; i<nodeNo; i+=2) {
        //i runs in mid order
        j = lastInternalNodeState(alignment, i);
        if(j != INT_MAX) {
            states[i] = j;
        }
    }
}

void shrinkAlignment(struct List *linearAlignment, int64_t columnGap) {
    struct TreeNodeHolder *treeNodeHolder;
    struct TreeNode *treeNode;

    while(linearAlignment->length > 0 && columnGap > 0) {
        treeNodeHolder = linearAlignment->list[--linearAlignment->length];
        treeNode = treeNodeHolder->treeNode;
        treeNode = getFirstNonSilentVertex(treeNode);
        if(treeNode != NULL) {
            columnGap--;
        }
    }
}

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//main script
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {
    int64_t i;
    int64_t j;
    //INT_32 k;
    //FLOAT_32 l;
    //float floatParser;
    char *newickTreeString = NULL;
    struct BinaryTree *newickTree;
    struct List *newickTreeLeafStrings;
    char **newickTreeStrings;
    struct BinaryTree **binaryTreeNodes;
    struct BinaryTree *binaryTreeNode;
    struct SubModel **subModels;
    struct CombinedTransitionModel **combinedTransitionModels;
    int64_t *seqLengths;
    struct List *seqs;
    int64_t nodeNumber;
    int64_t seqNo;
    int64_t seqInc = 1;
    int64_t startSeq = 0;
    char **inputSequences = NULL;
    int64_t inputSequenceNumber = INT_MAX;
    FILE *outputAlignment = NULL;

    char *constraintAlignmentFile = NULL;
    struct SequenceGraph **inputSequenceGraphs;
    struct Constraints ***primeConstraints;
    struct SequenceGraph *sequenceGraph;
    struct ParameterStruct *parameterStruct;

    float totalScore;
    float *ancestorProbs = NULL;
    struct List *alignment;
    char *mod;
    struct List *linearAlignment;
    struct List *wvAlignment;
    char *repeatBases;
    int64_t *treeStates;
    FILE *treeStatesFile;
    FILE *scoreOutputFile;

    int64_t seedTime = 0;

    if(argc == 1) {
        fprintf(stderr, "Ortheus [MODIFIER_ARGUMENTS]\n");
        fprintf(stderr, "A program for the inferral of ancestor sequences\n");
        fprintf(stderr, "Arguments:\n");
        fprintf(stderr, "\t-a [FILE]xN input sequence files (argument is essential)\n");
        fprintf(stderr, "\t-b [STRING] newick tree string (argument is essential, parser is pretty tolerant, but maybe dangerous)\n");
        fprintf(stderr, "\t-c [FILE] constraining alignment (if not present, assumes not constraints)\n");
        fprintf(stderr, "\t-d [FILE] output alignment file (if not present, writes to standard out)\n");
        fprintf(stderr, "\t-e set logging level to DEBUG (lowest level stuff, default is OFF)\n");
        fprintf(stderr, "\t-f set logging level to INFO (default is OFF)\n");
        fprintf(stderr, "\t-g don't write MFA alignment\n");
        fprintf(stderr, "\t-h output only leaf sequence alignment \n");
        fprintf(stderr, "\t-i [INTEGER] number of alignment samples\n");
        fprintf(stderr, "\t-j [INTEGER] 2 * anti-diagonal relax value\n");
        fprintf(stderr, "\t-k [INTEGER] use given seed for random number generator \n");
        //fprintf(stderr, "\t-l use Jukes-Cantor instead of Tamura-Nei substitution model\n");
        //fprintf(stderr, "\t-m [FLOAT] transition : tranversion ratio \n");
        //fprintf(stderr, "\t-n [FLOAT]x[ACGT] use following nucleotide frequencies (otherwise expected frequencies have stationary GC of 40 percent) \n");
        //fprintf(stderr, "\t-o [FLOAT] insertion open prob (scaled linearly with branch length) \n");
        //fprintf(stderr, "\t-p [FLOAT] insertion continue prob \n");
        //fprintf(stderr, "\t-q [FLOAT] deletion open prob (scaled linearly with branch length) \n");
        //fprintf(stderr, "\t-r [FLOAT] deletion continue prob \n");
        fprintf(stderr, "\t-s [INTEGER] leave gap at end of the alignment equal to at least given number of columns \n");
        fprintf(stderr, "\t-t [FILE] tree states files for reading in \n");
        fprintf(stderr, "\t-u [FILE] tree states files for writing out to \n");
        fprintf(stderr, "\t-v use partial rather than total ordering on constraints imposed by input alignment \n");
        fprintf(stderr, "\t-w seed with current system time \n");
        fprintf(stderr, "\t-x output file for scores (otherwise write to command-line)\n");
        //Stuff for the params of the model
        printParamStruct();
        exit(0);
    }
    //parse arguments
    for(i=1; i<argc; i++) {
        mod = argv[i];
        if(mod[0] != '-') {
            continue;
        }
        assert(mod[0] == '-');
        switch(mod[1]) {
            case 'a':
                inputSequences = argv + ++i;
                inputSequenceNumber = 0;
                while(i < argc && argv[i][0] != '-') {
                    inputSequenceNumber++;
                    i++;
                }
                i--;
                break;
            case 'b':
                newickTreeString = argv[++i];
                break;
            case 'c':
                constraintAlignmentFile = argv[++i];
                break;
            case 'd':
                OUTPUT_ALIGNMENT_FILE = argv[++i];
                break;
            case 'e':
                st_setLogLevel(debug);
                break;
            case 'f':
                st_setLogLevel(info);
                break;
            case 'g':
                WRITE_MFA_ALIGNMENT = FALSE;
                break;
            case 'h':
                seqInc = 2;
                break;
            case 'i':
                sscanf(argv[++i], INT_STRING, &NUMBER_OF_SAMPLES);
                break;
            case 'j':
                sscanf(argv[++i], INT_STRING, &CONSTRAINT_HALF_ANTI_DIAGONAL_LOOSEN_SIZE);
                break;
            case 'k':
                sscanf(argv[++i], INT_STRING, &j);
                srand(j);
                break;
            case 's':
                sscanf(argv[++i], INT_STRING, &VITERBI_COLUMN_GAP);
                break;
            case 't':
                TREE_STATES_INPUT_FILE = argv[++i];
                break;
            case 'u':
                TREE_STATES_OUTPUT_FILE = argv[++i];
                break;
            case 'v':
                TOTAL_CONSTRAINTS = FALSE;
                break;
            case 'w':
                seedTime = (unsigned)time(NULL);
                break;
            case 'x':
                SCORE_OUTPUT_FILE = argv[++i];
                break;
            default:
                break;
        }
    }
    parameterStruct = constructParamStruct(argc, argv);

    for(i=0; i<argc; i++) {
        st_logInfo("Argument recieved, no : " INT_STRING ", value : %s\n", i, argv[i]);
    }
    st_logInfo("Program seeded with time %" PRIi64 " ", seedTime);
    srand(seedTime);

    st_logInfo("Anti diagonal relax value " INT_STRING " \n", CONSTRAINT_HALF_ANTI_DIAGONAL_LOOSEN_SIZE);
    st_logInfo("Number of samples " INT_STRING " \n", NUMBER_OF_SAMPLES);

    //tree stuff
    st_logInfo("Newick-Tree : %s\n", newickTreeString);
    newickTree = newickTreeParser(newickTreeString, DEFAULT_BINARY_TREE_DISTANCE, 0);
    newickTreeLeafStrings = binaryTree_getOrderedLeafStrings(newickTree);
    binaryTree_depthFirstNumbers(newickTree);
    st_logDebug("Parsed tree\n");
    nodeNumber = newickTree->traversalID->midEnd;
    newickTreeStrings = st_malloc(sizeof(void *)*nodeNumber);
    treeNodeNames(newickTree, newickTreeStrings, (char **)newickTreeLeafStrings->list);
    binaryTreeNodes = st_malloc(sizeof(void *)*nodeNumber);
    subModels = st_malloc(sizeof(void *)*nodeNumber);
    getBinaryTreeNodesInMidOrder(newickTree, binaryTreeNodes);
    if(st_getLogLevel() == debug) {
        printBinaryTree(stderr, newickTree);
    }
    st_logInfo("Newick-Tree seems okay\n");
    //parse sequences
    seqs = parseSequences(inputSequences, inputSequenceNumber, &seqNo, &seqLengths);
    combinedTransitionModels = st_malloc(sizeof(void *)*nodeNumber);
    for(i=0; i<nodeNumber; i++) {
        st_logInfo("Newick-Tree node names : %s\n", newickTreeStrings[i]);
        binaryTreeNode = binaryTreeNodes[i];
        if(binaryTreeNode->internal) {
            combinedTransitionModels[i] = constructCombinedTransitionModel(
            binaryTreeNode->left->distance, binaryTreeNode->right->distance, binaryTreeNode == newickTree, parameterStruct);
            ancestorProbs = combinedTransitionModels[i]->ancestorProbs;
            subModels[binaryTreeNode->left->traversalID->mid] = combinedTransitionModels[i]->subModelX;
            subModels[binaryTreeNode->right->traversalID->mid] = combinedTransitionModels[i]->subModelY;
        }
        else {
            combinedTransitionModels[i] = NULL;
        }
    }
    //get leaf sequence graphs
    inputSequenceGraphs = buildLeafSequenceGraphs(seqs, seqNo, seqLengths, nodeNumber, binaryTreeNodes);

    assert(1 + nodeNumber/2  == seqNo);
    st_logInfo("Sequence number tallies with tree : " INT_STRING "\n", seqNo);
    for(i=0; i<seqNo; i++) {
        st_logInfo("Sequence lengths : " INT_STRING " : %s \n", seqLengths[i], newickTreeLeafStrings->list[i]);
    }
    //get prime constraints
    primeConstraints = buildPrimeConstraints(constraintAlignmentFile, seqNo, seqLengths, newickTree);
    //do input tree states
    if(TREE_STATES_INPUT_FILE != NULL) {
        treeStates = st_malloc(sizeof(int64_t)*nodeNumber);
        treeStatesFile = fopen(TREE_STATES_INPUT_FILE, "r");
        readIntegers(treeStatesFile, nodeNumber, treeStates);
        fclose(treeStatesFile);
        for(i=0; i<nodeNumber; i++) {
            st_logInfo("Tree states read : " INT_STRING " \n", treeStates[i]);
        }
    }
    else {
        treeStates = NULL;
    }
    //now do the final alignment
    i=0;
    sequenceGraph = \
    align_BottomUpScript(newickTree, inputSequenceGraphs,
                         primeConstraints,
                         newickTreeStrings, //subModels,
                         combinedTransitionModels,
                         NUMBER_OF_SAMPLES, 0, &i, treeStates);
    st_logInfo("Alignment graph computed\n");
    alignment = viterbi(sequenceGraph, &totalScore);
    st_logInfo("Computed viterbi alignment\n");
    linearAlignment = lineariseAlignment(alignment);
    st_logInfo("Mapped alignment data-structure to columns for output \n");
    shrinkAlignment(linearAlignment, VITERBI_COLUMN_GAP);
    wvAlignment = convertAlignmentToColumns(linearAlignment, nodeNumber, subModels, ancestorProbs);
    st_logInfo("Created bases for alignment\n");
    repeatBases = calculateRepeatBases((char **)seqs->list, (float **)wvAlignment->list, seqNo, wvAlignment->length);
    if(WRITE_MFA_ALIGNMENT) {
        if(OUTPUT_ALIGNMENT_FILE != NULL) {
            outputAlignment = fopen(OUTPUT_ALIGNMENT_FILE, "w");
            writeFastaAlignment(outputAlignment, (float **)wvAlignment->list, wvAlignment->length, newickTreeStrings, startSeq, seqInc, nodeNumber, dNAMap_WVToACTGFn, repeatBases);
            fclose(outputAlignment);
        }
        else {
            writeFastaAlignment(stdout, (float **)wvAlignment->list, wvAlignment->length, newickTreeStrings, startSeq, seqInc, nodeNumber, dNAMap_WVToACTGFn, repeatBases);
        }
    }
    if(TREE_STATES_OUTPUT_FILE != NULL) {
        if(treeStates == NULL) {
            treeStates = st_malloc(sizeof(int64_t)*nodeNumber);
            memset(treeStates, 0, sizeof(int64_t)*nodeNumber); //assumes start state is 0 (but then you won't get much using this with zero length alignments!)
        }
        calculateTreeStates(linearAlignment, nodeNumber, treeStates);
        for(i=1; i<nodeNumber; i+=2) {
            if(VITERBI_COLUMN_GAP == 0) { //to avoid counting INT_32o end state
                convertTransitionIDToStates(stateNo(), treeStates[i], &treeStates[i], &j); //j is dummy
            }
            else {
                convertTransitionIDToStates(stateNo(), treeStates[i], &j, &treeStates[i]); //j is dummy
            }
        }
        for(i=0; i<nodeNumber; i++) {
            st_logInfo("Tree states to write : " INT_STRING " \n", treeStates[i]);
        }
        treeStatesFile = fopen(TREE_STATES_OUTPUT_FILE, "w");
        writeIntegers(treeStatesFile, nodeNumber, treeStates);
        fclose(treeStatesFile);
    }
    //memory clean up
    destructParamStruct(parameterStruct);
    //get rid of seqs now, as no longe needed
    destructList(seqs);
    if(treeStates != NULL) {
        free(treeStates);
    }
    //struct BinaryTree *newickTree;
    destructBinaryTree(newickTree);
    //struct List *newickTreeLeafStrings;
    newickTreeLeafStrings->destructElement = NULL; //strings are freed below
    destructList(newickTreeLeafStrings);
    //char **newickTreeStrings;
    for(i=0; i<nodeNumber; i++) {
        free(newickTreeStrings[i]);
    }
    free(newickTreeStrings);
    //struct BinaryTree **binaryTreeNodes;
    free(binaryTreeNodes); //elements already free from binary tree destruction
    free(subModels);
    for(i=0; i<nodeNumber; i++) {
        if(combinedTransitionModels[i] != NULL) {
            destructCombinedTransitionModel(combinedTransitionModels[i]);
        }
    }
    free(combinedTransitionModels);
    //INT_32 *seqLengths;
    free(seqLengths);

    //struct SequenceGraph **inputSequenceGraphs;
    free(inputSequenceGraphs); //have already been cleared
    //struct Constraints ***primeConstraints;
    for(i=0; i<seqNo; i++) {
        for(j=0; j<seqNo; j++) {
            destructConstraints(primeConstraints[i][j]);
        }
    }

    //struct SequenceGraph *sequenceGraph;
    destructSequenceGraph(sequenceGraph, TRUE);
    //struct List *alignment;
    destructList(alignment);
    destructList(linearAlignment);
    //char *mod;
    //struct List *wvAlignment;
    destructList(wvAlignment);
    free(repeatBases);
    //end memory clean up
    //extern INT_32 heapCounter;
    st_logInfo("Finished, total time taken : %f (seconds) \n", ((clock() + 0.0f)/CLOCKS_PER_SEC + 0.0f));

    if(SCORE_OUTPUT_FILE != NULL) {
        scoreOutputFile = fopen(SCORE_OUTPUT_FILE, "w");
        fprintf(scoreOutputFile, "%f\n", totalScore);
        fclose(scoreOutputFile);
    }
    else {
        printf("Total-score: %f \n", totalScore);
    }
    //while(TRUE)
    //    ;
    return 0;
}
