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
    int32_t refCount;

    int32_t type;
    int32_t transitionID;
    struct TraversalID *traversalID;
    struct TreeNode *treeNodeX;
    struct TreeNode *treeNodeY;
    float *wV;
};

struct TreeNode *copyConstructTreeNode(int32_t type, int32_t transitionID, struct TraversalID *traversalID,
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

struct Edge {
    //edge in sequence graph

    int32_t from;
    int32_t to;
    float edgeScore;
    float subScore;
    float *wV;
    int32_t silent;
    struct TreeNode *treeNode;
    int32_t iD;
};

struct Edge *copyConstructEdge(int32_t from, int32_t to, float edgeScore, float subScore,
                               float *wV, int32_t silent, void *treeNode, int32_t iD);

void destructEdge(struct Edge *edge);

struct TraceBackEdge {
    int64_t from;
    int64_t to;
    float edgeScore;
    struct Edge *edgeX;
    struct Edge *edgeY;
    char silent;
    void *getTreeNode;
};

struct TraceBackEdge *constructTraceBackEdge(int64_t from, int64_t to, float edgeScore, struct Edge *edgeX, struct Edge *edgeY, char silent, void *getTreeNode);

void destructTraceBackEdge(struct TraceBackEdge *edge);

//sequence graphs

struct SequenceGraph {
    int32_t vertexNo;
    struct List *edges;
    struct List **edgesArrangedByToVertex;
    struct List **edgesArrangedByFromVertex;
};

struct SequenceGraph *constructSequenceGraph(struct List *edges, int32_t vertexNo);

void destructSequenceGraph(struct SequenceGraph *sequenceGraph, int32_t freeEdgeList);

//graph member holders

struct GraphMemberHolder {
    void *graphMember;
    int32_t *sequenceConstraints;
    void (*destructGraphMember)(void *);
};

struct GraphMemberHolder *constructGraphMember(void *graphMember, int32_t *sequenceConstraints, void (*destructGraphMember)(void *));

void destructGraphMember(struct GraphMemberHolder *graphMemberHolder);

int32_t  isSilent(int32_t state);
int32_t  isXInsert(int32_t state);
int32_t  isXDelete(int32_t state);
int32_t  isYInsert(int32_t state);
int32_t  isYDelete(int32_t state);
int32_t  isMatch(int32_t state);
int32_t  isXYDelete(int32_t state);

struct CombinedTransitionModel *constructCombinedTransitionModel(float distanceX, float distanceY, int32_t includeRoot, struct ParameterStruct *pS);

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
    int32_t leftMostSeqNoX;
    int32_t leftMostSeqNoY;
    int32_t leafSeqNoX;
    int32_t leafSeqNoY;
    //struct SubModel *subModel;
    int32_t numberOfSamples;

    struct Chunks *matrixChunks;
    float *fromCell;
    float *toCell;

    int32_t **vertexXSequenceCoordinates;
    //INT_32 **vertexYSequenceCoordinates;
    int32_t *sequenceGraphXSilentVertices_To;
    int32_t *sequenceGraphYSilentVertices_To;

    //these datastructures are used to compute the constraint envelope of the alignments during matrix computation
    void **mergedStartConstraints_Vertices;
    void **mergedStartConstraints_Edges;

    void **mergedEndConstraints_Vertices;
    void **mergedEndConstraints_Edges;

    void **mergedEndSequenceCoordinates_Vertices;
    void **mergedEndSequenceCoordinates_Edges;

    //added to avoid repeated access
    int32_t vertexXNo;
    int32_t vertexYNo;

    //used in scanning
    int32_t *newVertices;
    int32_t newVertices_Size;
    int32_t noOfNewVertices;
    int32_t *changes;
    int32_t noOfChanges;
    struct List *sequenceCoordinatesCollection;
    void **mergedEndConstraints_GraphMembers;
    int32_t (*getIDFromGraphMember)(void *graphMember);

    //transition starts and ends, used in compute matrix and tracebacl
    int64_t to;
    int64_t from;

    //following used by the traceback methods
    void **potentialEdges;
    float *potentialEdgeCosts;
    int32_t potentialEdges_Size;
    int32_t potentialEdges_Index;
    //branches for composing tree, used in traceback
    struct TreeNode *deleteNodeX;
    struct TreeNode *deleteNodeY;
    //struct SubModel *subModelX;
    //struct SubModel *subModelY;
    int64_t toCombined;
    int32_t state;
    void *getTreeNode;
    char silent;
    struct Edge *edgeX;
    struct Edge *edgeY;
    int32_t *treeStates;
};

inline int32_t stateNo();

float *startStates(struct CombinedTransitionModel *model);

float *endStates(struct CombinedTransitionModel *model);

//void turnOnDeleteXYLoopCorrection(struct CombinedTransitionModel *model);

//void turnOffDeleteXYLoopCorrection(struct CombinedTransitionModel *model);

inline void silentFn(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, float *cell,
					  void (*assignFn)(struct AlignmentDataStructures *aDS, int32_t, int32_t, float));

//inline void silentFn(struct AlignmentDataStructures *, struct CombinedTransitionModel *model, void (*assignFn)(struct AlignmentDataStructures *, INT_32, INT_32, FLOAT_32));

//inline void deleteFn(struct AlignmentDataStructures *, struct CombinedTransitionModel *model, void (*assignFn)(struct AlignmentDataStructures *, INT_32, INT_32, FLOAT_32));

//inline void deleteDeleteFn(struct AlignmentDataStructures *, struct CombinedTransitionModel *model, void (*assignFn)(struct AlignmentDataStructures *, INT_32, INT_32, FLOAT_32));

inline void silentFn_TraceBack(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, void (*assignFn)(struct AlignmentDataStructures *aDS, int32_t, int32_t, float, float));

inline void deleteFn_TraceBack(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, void (*assignFn)(struct AlignmentDataStructures *aDS, int32_t, int32_t, float, float));

inline void insertXFn(struct AlignmentDataStructures *, struct CombinedTransitionModel *model, struct Edge *edge,
                       void (*assignFn)(struct AlignmentDataStructures *, int32_t, int32_t, float));

inline void insertXFn_TraceBack(struct AlignmentDataStructures *, struct CombinedTransitionModel *model, struct Edge *edge,
                       void (*assignFn)(struct AlignmentDataStructures *, int32_t, int32_t, float, float));

inline void insertYFn(struct AlignmentDataStructures *, struct CombinedTransitionModel *model, struct Edge *edge,
                      void (*assignFn)(struct AlignmentDataStructures *, int32_t, int32_t, float));

inline void insertYFn_TraceBack(struct AlignmentDataStructures *, struct CombinedTransitionModel *model, struct Edge *edge,
                      void (*assignFn)(struct AlignmentDataStructures *, int32_t, int32_t, float, float));

inline void deleteXFn(struct AlignmentDataStructures *, struct CombinedTransitionModel *model, struct Edge *edge,
                      void (*assignFn)(struct AlignmentDataStructures *, int32_t, int32_t, float));

inline void deleteXFn_TraceBack(struct AlignmentDataStructures *, struct CombinedTransitionModel *model, struct Edge *edge,
                      void (*assignFn)(struct AlignmentDataStructures *, int32_t, int32_t, float, float));

inline void deleteYFn(struct AlignmentDataStructures *, struct CombinedTransitionModel *model, struct Edge *edge,
                      void (*assignFn)(struct AlignmentDataStructures *, int32_t, int32_t, float));

inline void deleteYFn_TraceBack(struct AlignmentDataStructures *, struct CombinedTransitionModel *model, struct Edge *edge,
                      void (*assignFn)(struct AlignmentDataStructures *, int32_t, int32_t, float, float));

inline void matchFn(struct AlignmentDataStructures *, struct CombinedTransitionModel *model,
            struct Edge *edgeX, struct Edge *edgeY, void (*assignFn)(struct AlignmentDataStructures *, int32_t, int32_t, float));

inline void matchFn_TraceBack(struct AlignmentDataStructures *, struct CombinedTransitionModel *model,
            struct Edge *edgeX, struct Edge *edgeY, void (*assignFn)(struct AlignmentDataStructures *, int32_t, int32_t, float, float));

//Functions for input output of transducers.
inline void printParamStruct();

struct ParameterStruct *constructParamStruct(int argc, char *argv[]);

void destructParamStruct(struct ParameterStruct *pM);

int32_t  *randomChoices(float *probs, int32_t sizeA, int32_t pathWeight);

struct SequenceGraph *computeEdgeGraph(struct SequenceGraph *sequenceGraphX, struct SequenceGraph *sequenceGraphY,
                        struct CombinedTransitionModel *model,
                        struct Constraints ***allConstraints,
                          //struct SubModel **subModels,
                          struct TraversalID *traversalID,
                          struct TraversalID *traversalIDX,
                          struct TraversalID *traversalIDY,
                          int32_t leftMostSeqNo, int32_t leafSeqNoX, int32_t leafSeqNoY,
                          int32_t numberOfSamples, int32_t *treeStates);

struct SequenceGraph *align_BottomUpScript(struct BinaryTree *binaryTree, struct SequenceGraph **sequenceGraphs, struct Constraints ***constraints,
                                           char **nodeNames, //struct SubModel **subModels,
                                           struct CombinedTransitionModel **combinedTransitionModels,
                                           int32_t numberOfSamples, int32_t leftMostSeqNo, int32_t *leafNo, int32_t *treeStates);

//FLOAT_32 *calculateAncestorSurvivalMatrix(char **alignment, INT_32 alignmentLength, struct BinaryTree *binaryTree, INT_32 seqNo);

//FLOAT_32 *calculateEdgeInsertionSurvivalBias(FLOAT_32 *ancestorSurvivalMatrix, INT_32 *sequenceLengths,
//                                             INT_32 internalNode, INT_32 alignmentLength, INT_32 seqNo, INT_32 leftMostSeqNo,
//                                             INT_32 **alignmentCoordinates, INT_32 **vertexSequenceCoordinates, struct SequenceGraph *sequenceGraph);

int32_t **calculateAlignmentCoordinates(char **alignment, int32_t alignmentLength, int32_t *seqLengths, int32_t seqNo);

struct List *viterbi(struct SequenceGraph *sequenceGraph, float *finalScore);

void convertTransitionIDToStates(int32_t stateNo, int32_t z, int32_t *fromState, int32_t *toState);

//felsensteins classic
float *felsensteins(struct TreeNode *treeNode, struct SubModel **subModels, int32_t nodeNumber, float *ancestorProbs);

#endif /*SEQUENCEGRAPHC_H_*/
