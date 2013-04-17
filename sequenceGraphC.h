/*
 * Copyright (C) 2008-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#ifndef SEQUENCEGRAPHC_H_
#define SEQUENCEGRAPHC_H_

#include <stdlib.h>
#include <limits.h>

#include "fastCMaths.h"
#include "commonC.h"
#include "substitutionC.h"

#include "constraintsC.h"
#include "xyzModelC.h"

#define ALPHABET_SIZE 4

//tree nodes
struct TreeNode {
    struct TreeNode *left;
    int64_t refCount;

    int64_t type;
    int64_t transitionID;
    struct TraversalID *traversalID;
    struct TreeNode *treeNodeX;
    struct TreeNode *treeNodeY;
    float *wV;
};

struct TreeNode *copyConstructTreeNode(int64_t type, int64_t transitionID, struct TraversalID *traversalID,
                                   struct TreeNode *treeNodeX, struct TreeNode *treeNodeY, float *wV);

void destructTreeNode(struct TreeNode *treeNode);

#define TREE_NODE_EFFECTIVELY_SILENT 64

#define TREE_NODE_INTERNAL 1
#define TREE_NODE_INSERT 66 //2 | 64;
#define TREE_NODE_DELETE 68 ////4 | 64;
#define TREE_NODE_LEAF 8
#define TREE_NODE_SILENT 80 //16 | 64;
#define TREE_NODE_PREVIOUSLY_SILENT 96 //32 | 64;

//edge

struct AlignmentDataStructures;

struct Edge {
    //edge in sequence graph

    int64_t from;
    int64_t to;
    float edgeScore;
    float subScore;
    float *wV;
    int64_t silent;
    struct TreeNode *treeNode;
    int64_t iD;
};

struct Edge *copyConstructEdge(int64_t from, int64_t to, float edgeScore, float subScore,
                               float *wV, int64_t silent, void *treeNode, int64_t iD);

void destructEdge(struct Edge *edge);

struct TraceBackEdge {
    int64_t from;
    int64_t to;
    float edgeScore;
    struct Edge *edgeX;
    struct Edge *edgeY;
    char silent;
    //void *getTreeNode;
    struct TreeNode *(*getTreeNode)(struct AlignmentDataStructures *, struct TraceBackEdge *, int64_t);
};


//struct TraceBackEdge *constructTraceBackEdge(int64_t from, int64_t to, float edgeScore, struct Edge *edgeX, struct Edge *edgeY, char silent, void *getTreeNode);

struct TraceBackEdge *constructTraceBackEdge(int64_t from, int64_t to, float edgeScore, struct Edge *edgeX, struct Edge *edgeY, char silent,
        struct TreeNode *(*getTreeNode)(struct AlignmentDataStructures *, struct TraceBackEdge *, int64_t));

void destructTraceBackEdge(struct TraceBackEdge *edge);

//sequence graphs

struct SequenceGraph {
    int64_t vertexNo;
    struct List *edges;
    struct List **edgesArrangedByToVertex;
    struct List **edgesArrangedByFromVertex;
};

struct SequenceGraph *constructSequenceGraph(struct List *edges, int64_t vertexNo);

void destructSequenceGraph(struct SequenceGraph *sequenceGraph, int64_t freeEdgeList);

//graph member holders

struct GraphMemberHolder {
    void *graphMember;
    int64_t *sequenceConstraints;
    void (*destructGraphMember)(void *);
};

struct GraphMemberHolder *constructGraphMember(void *graphMember, int64_t *sequenceConstraints, void (*destructGraphMember)(void *));

void destructGraphMember(struct GraphMemberHolder *graphMemberHolder);

int64_t  isSilent(int64_t state);
int64_t  isXInsert(int64_t state);
int64_t  isXDelete(int64_t state);
int64_t  isYInsert(int64_t state);
int64_t  isYDelete(int64_t state);
int64_t  isMatch(int64_t state);
int64_t  isXYDelete(int64_t state);

struct CombinedTransitionModel *constructCombinedTransitionModel(float distanceX, float distanceY, int64_t includeRoot, struct ParameterStruct *pS);

void destructCombinedTransitionModel(struct CombinedTransitionModel *model);

//end model
//over all data structure

struct AlignmentDataStructures {
    struct SequenceGraph *sequenceGraphX;
    struct SequenceGraph *sequenceGraphY;
    //void *allConstraints;
    struct CombinedTransitionModel *model;
    float *startStates;
    float *endStates;
    //void *traversalIDX;
    //void *traversalIDY;
    struct TraversalID *traversalID;
    int64_t leftMostSeqNoX;
    int64_t leftMostSeqNoY;
    int64_t leafSeqNoX;
    int64_t leafSeqNoY;
    //struct SubModel *subModel;
    int64_t numberOfSamples;

    struct Chunks *matrixChunks;
    float *fromCell;
    float *toCell;

    int64_t **vertexXSequenceCoordinates;
    //INT_32 **vertexYSequenceCoordinates;
    int64_t *sequenceGraphXSilentVertices_To;
    int64_t *sequenceGraphYSilentVertices_To;

    //these datastructures are used to compute the constraint envelope of the alignments during matrix computation
    void **mergedStartConstraints_Vertices;
    void **mergedStartConstraints_Edges;

    void **mergedEndConstraints_Vertices;
    void **mergedEndConstraints_Edges;

    void **mergedEndSequenceCoordinates_Vertices;
    void **mergedEndSequenceCoordinates_Edges;

    //added to avoid repeated access
    int64_t vertexXNo;
    int64_t vertexYNo;

    //used in scanning
    int64_t *newVertices;
    int64_t newVertices_Size;
    int64_t noOfNewVertices;
    int64_t *changes;
    int64_t noOfChanges;
    struct List *sequenceCoordinatesCollection;
    void **mergedEndConstraints_GraphMembers;
    int64_t (*getIDFromGraphMember)(void *graphMember);

    //transition starts and ends, used in compute matrix and tracebacl
    int64_t to;
    int64_t from;

    //following used by the traceback methods
    void **potentialEdges;
    float *potentialEdgeCosts;
    int64_t potentialEdges_Size;
    int64_t potentialEdges_Index;
    //branches for composing tree, used in traceback
    struct TreeNode *deleteNodeX;
    struct TreeNode *deleteNodeY;
    //struct SubModel *subModelX;
    //struct SubModel *subModelY;
    int64_t toCombined;
    int64_t state;
    struct TreeNode *(*getTreeNode)(struct AlignmentDataStructures *, struct TraceBackEdge*, int64_t);
    //void *getTreeNode;
    char silent;
    struct Edge *edgeX;
    struct Edge *edgeY;
    int64_t *treeStates;
};

int64_t stateNo();

float *startStates(struct CombinedTransitionModel *model);

float *endStates(struct CombinedTransitionModel *model);

//void turnOnDeleteXYLoopCorrection(struct CombinedTransitionModel *model);

//void turnOffDeleteXYLoopCorrection(struct CombinedTransitionModel *model);

void silentFn(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, float *cell,
                      void (*assignFn)(struct AlignmentDataStructures *aDS, int64_t, int64_t, float));

//void silentFn(struct AlignmentDataStructures *, struct CombinedTransitionModel *model, void (*assignFn)(struct AlignmentDataStructures *, INT_32, INT_32, FLOAT_32));

//void deleteFn(struct AlignmentDataStructures *, struct CombinedTransitionModel *model, void (*assignFn)(struct AlignmentDataStructures *, INT_32, INT_32, FLOAT_32));

//void deleteDeleteFn(struct AlignmentDataStructures *, struct CombinedTransitionModel *model, void (*assignFn)(struct AlignmentDataStructures *, INT_32, INT_32, FLOAT_32));

void silentFn_TraceBack(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, void (*assignFn)(struct AlignmentDataStructures *aDS, int64_t, int64_t, float, float));

void deleteFn_TraceBack(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, void (*assignFn)(struct AlignmentDataStructures *aDS, int64_t, int64_t, float, float));

void insertXFn(struct AlignmentDataStructures *, struct CombinedTransitionModel *model, struct Edge *edge,
                       void (*assignFn)(struct AlignmentDataStructures *, int64_t, int64_t, float));

void insertXFn_TraceBack(struct AlignmentDataStructures *, struct CombinedTransitionModel *model, struct Edge *edge,
                       void (*assignFn)(struct AlignmentDataStructures *, int64_t, int64_t, float, float));

void insertYFn(struct AlignmentDataStructures *, struct CombinedTransitionModel *model, struct Edge *edge,
                      void (*assignFn)(struct AlignmentDataStructures *, int64_t, int64_t, float));

void insertYFn_TraceBack(struct AlignmentDataStructures *, struct CombinedTransitionModel *model, struct Edge *edge,
                      void (*assignFn)(struct AlignmentDataStructures *, int64_t, int64_t, float, float));

void deleteXFn(struct AlignmentDataStructures *, struct CombinedTransitionModel *model, struct Edge *edge,
                      void (*assignFn)(struct AlignmentDataStructures *, int64_t, int64_t, float));

void deleteXFn_TraceBack(struct AlignmentDataStructures *, struct CombinedTransitionModel *model, struct Edge *edge,
                      void (*assignFn)(struct AlignmentDataStructures *, int64_t, int64_t, float, float));

void deleteYFn(struct AlignmentDataStructures *, struct CombinedTransitionModel *model, struct Edge *edge,
                      void (*assignFn)(struct AlignmentDataStructures *, int64_t, int64_t, float));

void deleteYFn_TraceBack(struct AlignmentDataStructures *, struct CombinedTransitionModel *model, struct Edge *edge,
                      void (*assignFn)(struct AlignmentDataStructures *, int64_t, int64_t, float, float));

void matchFn(struct AlignmentDataStructures *, struct CombinedTransitionModel *model,
            struct Edge *edgeX, struct Edge *edgeY, void (*assignFn)(struct AlignmentDataStructures *, int64_t, int64_t, float));

void matchFn_TraceBack(struct AlignmentDataStructures *, struct CombinedTransitionModel *model,
            struct Edge *edgeX, struct Edge *edgeY, void (*assignFn)(struct AlignmentDataStructures *, int64_t, int64_t, float, float));

//Functions for input output of transducers.
void printParamStruct();

struct ParameterStruct *constructParamStruct(int argc, char *argv[]);

void destructParamStruct(struct ParameterStruct *pM);

int64_t  *randomChoices(float *probs, int64_t sizeA, int64_t pathWeight);

struct SequenceGraph *computeEdgeGraph(struct SequenceGraph *sequenceGraphX, struct SequenceGraph *sequenceGraphY,
                        struct CombinedTransitionModel *model,
                        struct Constraints ***allConstraints,
                          //struct SubModel **subModels,
                          struct TraversalID *traversalID,
                          struct TraversalID *traversalIDX,
                          struct TraversalID *traversalIDY,
                          int64_t leftMostSeqNo, int64_t leafSeqNoX, int64_t leafSeqNoY,
                          int64_t numberOfSamples, int64_t *treeStates);

struct SequenceGraph *align_BottomUpScript(struct BinaryTree *binaryTree, struct SequenceGraph **sequenceGraphs, struct Constraints ***constraints,
                                           char **nodeNames, //struct SubModel **subModels,
                                           struct CombinedTransitionModel **combinedTransitionModels,
                                           int64_t numberOfSamples, int64_t leftMostSeqNo, int64_t *leafNo, int64_t *treeStates);

//FLOAT_32 *calculateAncestorSurvivalMatrix(char **alignment, INT_32 alignmentLength, struct BinaryTree *binaryTree, INT_32 seqNo);

//FLOAT_32 *calculateEdgeInsertionSurvivalBias(FLOAT_32 *ancestorSurvivalMatrix, INT_32 *sequenceLengths,
//                                             INT_32 internalNode, INT_32 alignmentLength, INT_32 seqNo, INT_32 leftMostSeqNo,
//                                             INT_32 **alignmentCoordinates, INT_32 **vertexSequenceCoordinates, struct SequenceGraph *sequenceGraph);

int64_t **calculateAlignmentCoordinates(char **alignment, int64_t alignmentLength, int64_t *seqLengths, int64_t seqNo);

struct List *viterbi(struct SequenceGraph *sequenceGraph, float *finalScore);

void convertTransitionIDToStates(int64_t stateNo, int64_t z, int64_t *fromState, int64_t *toState);

//felsensteins classic
float *felsensteins(struct TreeNode *treeNode, struct SubModel **subModels, int64_t nodeNumber, float *ancestorProbs);

#endif /*SEQUENCEGRAPHC_H_*/
