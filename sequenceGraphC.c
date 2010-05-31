#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>

#include "fastCMaths.h"
#include "hashTableC.h"
#include "heapC.h"
#include "commonC.h"
#include "substitutionC.h"

#include "xyzModelC.h"
#include "sequenceGraphC.h"
#include "constraintsC.h"

#define CELL_TABLE_INITIAL_SIZE_PER_EDGE 10

void calculateConstraints_GreaterThan(int32_t *leafPositionsX, int32_t leftMostLeafNoX, int32_t leftMostLeafNoY, int32_t leafSeqNoX, int32_t leafSeqNoY,
                                     struct Constraints ***allConstraints, int32_t *constraints) {
    //Like calculateConstraints, except for 'x-is-greater-than-y' constraints.
    int32_t seqX;
    int32_t seqY;
    int32_t x, x2;
    int32_t dummy; //these values are dummies

    for(seqY=0; seqY<leafSeqNoY; seqY++) {
        //x = constraints[seqY];
        getYConstraint(allConstraints[leftMostLeafNoY + seqY][leftMostLeafNoX + 0], leafPositionsX[0], &x, &dummy, &dummy);
        for(seqX=1; seqX<leafSeqNoX; seqX++) {
            getYConstraint(allConstraints[leftMostLeafNoY + seqY][leftMostLeafNoX + seqX], leafPositionsX[seqX], &x2, &dummy, &dummy);
            if(x2 > x) {
                x = x2;
            }
        }
        constraints[seqY] = x;
    }
}

void calculateConstraints_LessThan(int32_t *leafPositionsX, int32_t leftMostLeafNoX, int32_t leftMostLeafNoY, int32_t leafSeqNoX, int32_t leafSeqNoY,
                                   struct Constraints ***allConstraints, int32_t *constraints) {
    //Finds the constraints for a given set of leaf positions (the x coordinates) on the y sequences.
    //Default function is for 'x-is-less-than-y' constraints
    //ignores constraINT_32 type
    int32_t seqX;
    int32_t seqY;
    int32_t y, y2;
    int32_t dummy; //these values are dummies

    for(seqY=0; seqY<leafSeqNoY; seqY++) {
        //y = constraints[seqY];
        getXConstraint(allConstraints[leftMostLeafNoX + 0][leftMostLeafNoY + seqY], leafPositionsX[0], &dummy, &y, &dummy);
        for(seqX=1; seqX<leafSeqNoX; seqX++) {
            getXConstraint(allConstraints[leftMostLeafNoX + seqX][leftMostLeafNoY + seqY], leafPositionsX[seqX], &dummy, &y2, &dummy);
            if(y2 < y) {
                y = y2;
            }
        }
        constraints[seqY] = y;
    }
}

void makeNonRedundantCollection(struct List *collection, int32_t (*legalCollectionFn)(int32_t *i, int32_t *j, int32_t k), int32_t length) {
    //Removes redundancy from a set. Note legalFunction is not necessarily ymmetric.
    int32_t i;
    int32_t j;
    void *k;
    void *l;

    static struct List list;

    listResize(&list, collection->length);
    list.length = 0;

    //for i in collection:
    for(i=0; i<collection->length; i++) {
        k = collection->list[i];
        //for j in collection:
        for(j=0; j<i; j++) {
            l = collection->list[j];
            if (!legalCollectionFn(k, l, length)) {
                goto outer;
            }
        }
        for(j=i+1; j<collection->length; j++) {
            l = collection->list[j];
            if ((!legalCollectionFn(k, l, length)) && memcmp(k, l, length*sizeof(int32_t))) {
                goto outer;
            }
        }
        list.list[list.length++] = k;
        outer:
        continue;
    }
    collection->length = list.length;
    memcpy(collection->list, list.list, list.length*sizeof(void *));
}

void filterCollection(struct List *collection, struct List *collection2, struct List *collection3, int32_t (*legalCollectionFn)(int32_t *i, int32_t *j, int32_t k), int32_t length) {
    //Removes redundancy from a set. Note legalFunction is not necessarily ymmetric.
    int32_t i;
    int32_t j;
    void *k;

    listResize(collection3, collection->length);
    collection3->length = 0;

    //for i in collection:
    for(i=0; i<collection->length; i++) {
        k = collection->list[i];
        //for j in collection:
        for(j=0; j<collection2->length; j++) {
            if (!legalCollectionFn(k, collection2->list[j], length)) {
                goto outer;
            }
        }
        collection3->list[collection3->length++] = k;
        outer:
        continue;
    }
}

void mergeNonRedundantCollections(struct List *collection1, struct List *collection2, int32_t (*legalCollectionFn)(int32_t *i, int32_t *j, int32_t k), int32_t length) {
    //Merges two sets of sequenceCoordinates, by finding minimal number of explainatory sets.
    static struct List collection3;
    int32_t i;

    filterCollection(collection1, collection2, &collection3, legalCollectionFn, length);
    filterCollection(collection2, &collection3, collection1, legalCollectionFn, length);
    listCopyResize(collection1, collection1->length + collection3.length);
    for(i=0; i<collection3.length; i++) {
        collection1->list[collection1->length++] = collection3.list[i];
    }
}

int32_t legalCollection_GreaterThan(int32_t *a, int32_t *b, int32_t length) {
    int32_t i;

    for(i=0; i<length; i++) {
        if (a[i] > b[i]) {
            return TRUE;
        }
    }
    return FALSE;
}

int32_t legalCollection_LessThan(int32_t *a, int32_t *b, int32_t length) {
    int32_t i;

    for(i=0; i<length; i++) {
        if (a[i] < b[i]) {
            return TRUE;
        }
    }
    return FALSE;
}

int32_t lessThan(int32_t i, int32_t j) {
    return i < j;
}

int32_t greaterThan(int32_t i, int32_t j) {
    return i > j;
}

int32_t edgeFrom(struct Edge *edge) {
    return edge->from;
}

int32_t edgeTo(struct Edge *edge) {
    return edge->to;
}

struct List *cloneCollection(struct List *oldList, int32_t length, struct Chunks *chunks) {
    int32_t i;
    struct List *newList;

    newList = constructEmptyList(oldList->length, NULL);
    for(i=0; i<oldList->length; i++) {
        newList->list[i] = memcpy(mallocChunk(chunks), oldList->list[i], sizeof(int32_t)*length);
    }
    return newList;
}

void getActiveLeaves(int32_t *activeLeaves, int32_t leftMostSeqNo, struct TreeNode *treeNode) {
    if(treeNode->type == TREE_NODE_INTERNAL) {
        getActiveLeaves(activeLeaves, leftMostSeqNo, treeNode->treeNodeX);
        getActiveLeaves(activeLeaves, leftMostSeqNo, treeNode->treeNodeY);
    }
    else if(treeNode->type == TREE_NODE_LEAF) {
        activeLeaves[treeNode->traversalID->leafNo - leftMostSeqNo] = TRUE;
    }
}

void calculateMergedSequenceCoordinates(int32_t **vertexSequenceCoordinates, struct SequenceGraph *sequenceGraph,
                                        int32_t leftMostSeqNo, int32_t leafSeqNo, void ***mergedVertices, void ***mergedEdges,
                                        int32_t vertexStart, int32_t vertexEnd, int32_t vertexIncrement, int32_t endIncrement,
                                        struct List **connections, int32_t (*legalCollectionFn)(int32_t *, int32_t *, int32_t),
                                        int32_t (*edgeTo)(struct Edge *), int32_t (*lessThan)(int32_t i, int32_t j), struct Chunks *chunks) {
    //Calculates the set of correct sequenceCoordinates for each vertex and non-silent edge, taking into
    //account the effect of skipping out sequences in the graph contained in insert
    //edges. Multiple sequence coordinate sets for each vertex are possible, due to multiple possible
    //paths to the right of a given vertex. Function will try to find minimum number possible.
    void **mergedSequenceCoordinates_Vertices;
    void **mergedSequenceCoordinates_Edges;
    int32_t vertex;
    int32_t i;
    int32_t j;
    int32_t seq;
    struct List *edges;
    struct Edge *edge;
    struct List *tempList;
    int32_t *tempInt;
    int32_t *activeLeaves;
    struct List *sequenceCoordinatesCollection;
    int32_t *sequenceCoordinatesTo;
    int32_t *sequenceCoordinatesFrom;
    int32_t *sequenceCoordinates;

    static struct List mergeList;

    mergedSequenceCoordinates_Vertices = st_malloc(sizeof(void *)*sequenceGraph->vertexNo);
    mergedSequenceCoordinates_Edges = st_malloc(sizeof(void *)*sequenceGraph->edges->length);
    activeLeaves = st_malloc(sizeof(int32_t)*leafSeqNo);

    //for vertex in vertexOrder(sequenceGraph):
    for(vertex=vertexStart; vertex!=vertexEnd; vertex+=vertexIncrement) {
        edges = connections[vertex];
        if (edges->length == 0) {
            tempList = constructEmptyList(1, NULL);
            mergedSequenceCoordinates_Vertices[vertex] = tempList;
            tempInt = mallocChunk(chunks);
            tempList->list[0] = tempInt;
            for(i=0; i<leafSeqNo; i++) {
                tempInt[i] = vertexSequenceCoordinates[vertex][i]+endIncrement;
            }
        }
        else {
            //for edge in edges:
            mergeList.length = 0;
            for(i=0; i<edges->length; i++) {
                edge = edges->list[i];
                if (!edge->silent) {
                    //modifies sequenceCoordinates from to vertices by transitions along edge
                    sequenceCoordinatesCollection = cloneCollection(mergedSequenceCoordinates_Vertices[edgeTo(edge)], leafSeqNo, chunks);
                    sequenceCoordinatesTo = vertexSequenceCoordinates[edge->to];
                    sequenceCoordinatesFrom = vertexSequenceCoordinates[edge->from];
                    //for seq in xrange(0, leafSeqNo):
                    memset(activeLeaves, FALSE, sizeof(int32_t)*leafSeqNo);
                    getActiveLeaves(activeLeaves, leftMostSeqNo, edge->treeNode);
                    for(seq=0; seq<leafSeqNo; seq++) {
                        if (sequenceCoordinatesTo[seq] - sequenceCoordinatesFrom[seq] > 0) {
                            //for sequenceCoordinates in sequenceCoordinatesCollection:
                            for(j=0; j<sequenceCoordinatesCollection->length; j++) {
                                sequenceCoordinates = sequenceCoordinatesCollection->list[j];
                                if (lessThan(sequenceCoordinatesTo[seq], sequenceCoordinates[seq]) && activeLeaves[seq]) { //this is correct, even though you probably think it's not
                                    sequenceCoordinates[seq] = sequenceCoordinatesTo[seq]; //sequenceCoordinatesTo[seq];
                                }
                            }
                        }
                    }
                    makeNonRedundantCollection(sequenceCoordinatesCollection, legalCollectionFn, leafSeqNo);
                    mergedSequenceCoordinates_Edges[edge->iD] = sequenceCoordinatesCollection;
                    mergeNonRedundantCollections(&mergeList, sequenceCoordinatesCollection, legalCollectionFn, leafSeqNo);
                }
                else {
                    mergedSequenceCoordinates_Edges[edge->iD] = mergedSequenceCoordinates_Vertices[edgeTo(edge)];
                    mergeNonRedundantCollections(&mergeList, mergedSequenceCoordinates_Vertices[edgeTo(edge)], legalCollectionFn, leafSeqNo);
                }
            }
            mergedSequenceCoordinates_Vertices[vertex] = copyConstructList(mergeList.list, mergeList.length, NULL);
        }
    }
    free(activeLeaves);
    *mergedVertices =  mergedSequenceCoordinates_Vertices;
    *mergedEdges =  mergedSequenceCoordinates_Edges;
}

void calculateMergedSequenceCoordinates_LeftToRight(int32_t **vertexSequenceCoordinates, struct SequenceGraph *sequenceGraph, int32_t leftMostSeqNo, int32_t leafSeqNo,
                                                   void ***mergedVertices, void ***mergedEdges, struct Chunks *chunks) {
    //As with calculateMergedSequenceCoordinates, but going from left to right
    //instead of right to left
    calculateMergedSequenceCoordinates(vertexSequenceCoordinates, sequenceGraph, leftMostSeqNo, leafSeqNo, mergedVertices, mergedEdges,
                                       0, sequenceGraph->vertexNo, 1, -1, sequenceGraph->edgesArrangedByToVertex,
                                       legalCollection_LessThan, edgeFrom, greaterThan, chunks);
}


void calculateMergedSequenceCoordinates_RightToLeft(int32_t **vertexSequenceCoordinates, struct SequenceGraph *sequenceGraph, int32_t leftMostSeqNo, int32_t leafSeqNo,
                                                   void ***mergedVertices, void ***mergedEdges, struct Chunks *chunks) {
    //As with calculateMergedSequenceCoordinates, but going from left to right
    //instead of right to left
    calculateMergedSequenceCoordinates(vertexSequenceCoordinates, sequenceGraph, leftMostSeqNo, leafSeqNo, mergedVertices, mergedEdges,
                                       sequenceGraph->vertexNo-1, -1, -1, 1, sequenceGraph->edgesArrangedByFromVertex,
                                       legalCollection_GreaterThan, edgeTo, lessThan, chunks);
}


void **convertSequenceCoordinatesToConstraints(void **mergedSequenceCoordinates, int32_t mergedSequenceCoordinateLength,
                                             int32_t leftMostLeafNoX, int32_t leftMostLeafNoY, int32_t leafSeqNoX, int32_t leafSeqNoY,
                                             struct Constraints ***allConstraints,
                                             void (*constraintsFn)(int32_t *, int32_t , int32_t , int32_t , int32_t , struct Constraints ***, int32_t *),
                                             int32_t (*legalCollectionFn)(int32_t *i, int32_t *j, int32_t length), struct Chunks *chunks) {
    //Converts sequence coordinates to constraints on other sequence. Used to convert merged
    //sequence coordinates into constraints
    int32_t i;
    int32_t j;
    struct List *collection;
    struct List *collection2;
    int32_t *temp;
    void **mergedSequenceConstraints;

    mergedSequenceConstraints = st_malloc(sizeof(void *)*mergedSequenceCoordinateLength);

    //for key in mergedSequenceCoordinates.keys():
    for(i=0; i<mergedSequenceCoordinateLength; i++) {

        collection = mergedSequenceCoordinates[i];
        collection2 = constructEmptyList(collection->length, NULL);
        mergedSequenceConstraints[i] = collection2;

        //for sequenceCoordinates in sequenceCoordinatesCollection:
        for(j=0; j<collection->length; j++) {
            temp = mallocChunk(chunks);//mallocLocal(sizeof(INT_32)*leafSeqNoY);
            constraintsFn(collection->list[j], leftMostLeafNoX, leftMostLeafNoY, leafSeqNoX, leafSeqNoY, allConstraints, temp);
            collection2->list[j] = temp;
        }
        makeNonRedundantCollection(collection2, legalCollectionFn, leafSeqNoY);
        //(struct List *collection, INT_32 (*legalCollectionFn)(INT_32 *i, INT_32 *j, INT_32 k), INT_32 length)
    }
    return mergedSequenceConstraints;
}

void **convertSequenceCoordinatesToConstraints_LeftToRight(void **mergedSequenceCoordinates, int32_t mergedSequenceCoordinateLength,
                                                         int32_t leftMostLeafNoX, int32_t leftMostLeafNoY, int32_t leafSeqNoX, int32_t leafSeqNoY,
                                                         struct Constraints ***allConstraints, struct Chunks *chunks) {
    return convertSequenceCoordinatesToConstraints(mergedSequenceCoordinates, mergedSequenceCoordinateLength,
                                                   leftMostLeafNoX, leftMostLeafNoY, leafSeqNoX, leafSeqNoY, allConstraints,
                                                   calculateConstraints_GreaterThan, legalCollection_LessThan, chunks);
}


void **convertSequenceCoordinatesToConstraints_RightToLeft(void **mergedSequenceCoordinates, int32_t mergedSequenceCoordinateLength,
                                                         int32_t leftMostLeafNoX, int32_t leftMostLeafNoY, int32_t leafSeqNoX, int32_t leafSeqNoY,
                                                         struct Constraints ***allConstraints, struct Chunks *chunks) {
    return convertSequenceCoordinatesToConstraints(mergedSequenceCoordinates, mergedSequenceCoordinateLength,
                                                   leftMostLeafNoX, leftMostLeafNoY, leafSeqNoX, leafSeqNoY, allConstraints,
                                                   calculateConstraints_LessThan, legalCollection_GreaterThan, chunks);
}

void cleanUpSequenceCoordinates(struct SequenceGraph *sequenceGraphX,
                                void **mergedSequenceCoordinates_Vertices,
                                void **mergedSequenceCoordinates_Edges) {
    int32_t i;
    //memory clean up
    for(i=0; i<sequenceGraphX->vertexNo; i++) {
        destructList(mergedSequenceCoordinates_Vertices[i]);
    }
    for(i=0; i<sequenceGraphX->edges->length; i++) {
        if(!((struct Edge *)sequenceGraphX->edges->list[i])->silent) {
            destructList(mergedSequenceCoordinates_Edges[i]);
        }
    }
    free(mergedSequenceCoordinates_Vertices);
    free(mergedSequenceCoordinates_Edges);
}

void calculateMergedConstraints_RightToLeft(int32_t **vertexXSequenceCoordinates,
                                            int32_t leftMostLeafNoX, int32_t leftMostLeafNoY, int32_t leafSeqNoX, int32_t leafSeqNoY,
                                            struct Constraints ***allConstraints, struct SequenceGraph *sequenceGraphX,
                                            void ***mergedVertices, void ***mergedEdges, struct Chunks *constraintChunks) {
    void **mergedSequenceCoordinates_Vertices;
    void **mergedSequenceCoordinates_Edges;
    struct Chunks *chunks;

    chunks = constructChunks(MEDIUM_CHUNK_SIZE, sizeof(int32_t)*leafSeqNoX);
    //Convenience function to wrap together two functions for calculating the constraints for each vertex
    //and edge in a graph in terms of a given set of sequences present on an opposite graph
    calculateMergedSequenceCoordinates_RightToLeft(vertexXSequenceCoordinates, sequenceGraphX, leftMostLeafNoX, leafSeqNoX,
                                                   &mergedSequenceCoordinates_Vertices, &mergedSequenceCoordinates_Edges, chunks);
    *mergedVertices = convertSequenceCoordinatesToConstraints_RightToLeft(mergedSequenceCoordinates_Vertices, sequenceGraphX->vertexNo,
                                                                          leftMostLeafNoX, leftMostLeafNoY, leafSeqNoX, leafSeqNoY,
                                                                          allConstraints, constraintChunks);
    *mergedEdges = convertSequenceCoordinatesToConstraints_RightToLeft(mergedSequenceCoordinates_Edges, sequenceGraphX->edges->length,
                                                                       leftMostLeafNoX, leftMostLeafNoY, leafSeqNoX, leafSeqNoY,
                                                                       allConstraints, constraintChunks);
    cleanUpSequenceCoordinates(sequenceGraphX, mergedSequenceCoordinates_Vertices, mergedSequenceCoordinates_Edges);
    destructChunks(chunks);
}

void calculateMergedConstraints_LeftToRight(int32_t **vertexXSequenceCoordinates,
                                            int32_t leftMostLeafNoX, int32_t leftMostLeafNoY, int32_t leafSeqNoX, int32_t leafSeqNoY,
                                            struct Constraints ***allConstraints, struct SequenceGraph *sequenceGraphX,
                                            void ***mergedVertices, void ***mergedEdges, struct Chunks *constraintChunks) {
    void **mergedSequenceCoordinates_Vertices;
    void **mergedSequenceCoordinates_Edges;
    struct Chunks *chunks;

    chunks = constructChunks(MEDIUM_CHUNK_SIZE, sizeof(int32_t)*leafSeqNoX);
    //Convenience function to wrap together two functions for calculating the constraints for each vertex
    //and edge in a graph in terms of a given set of sequences present on an opposite graph
    calculateMergedSequenceCoordinates_LeftToRight(vertexXSequenceCoordinates, sequenceGraphX, leftMostLeafNoX, leafSeqNoX,
                                                   &mergedSequenceCoordinates_Vertices, &mergedSequenceCoordinates_Edges, chunks);
    *mergedVertices = convertSequenceCoordinatesToConstraints_LeftToRight(mergedSequenceCoordinates_Vertices, sequenceGraphX->vertexNo, leftMostLeafNoX, leftMostLeafNoY, leafSeqNoX, leafSeqNoY,
                                                                          allConstraints, constraintChunks);
    *mergedEdges = convertSequenceCoordinatesToConstraints_LeftToRight(mergedSequenceCoordinates_Edges, sequenceGraphX->edges->length, leftMostLeafNoX, leftMostLeafNoY, leafSeqNoX, leafSeqNoY,
                                                                       allConstraints, constraintChunks);
    cleanUpSequenceCoordinates(sequenceGraphX, mergedSequenceCoordinates_Vertices, mergedSequenceCoordinates_Edges);
    destructChunks(chunks);
}

int32_t *calculateSilentVertices(struct SequenceGraph *sequenceGraph) {
    //labels vertices silent if has (only) silent incoming edges
    int32_t vertex;
    int32_t *silentVertices;
    int32_t i;
    struct Edge *edge;

    silentVertices = st_calloc(sequenceGraph->vertexNo, sizeof(int32_t));
    //for vertex in xrange(0, sequenceGraph.vertexNo):
    for(vertex=0; vertex<sequenceGraph->vertexNo; vertex++) {
        //for edge in sequenceGraph.edgesArrangedByToVertex[vertex]:
        for(i=0; i<sequenceGraph->edgesArrangedByToVertex[vertex]->length; i++) {
            edge = sequenceGraph->edgesArrangedByToVertex[vertex]->list[i];
            if (edge->silent) {
                silentVertices[vertex] = TRUE;
                break;
            }
        }
    }
    return silentVertices;
}

//treeNodes
struct TreeNode *copyConstructTreeNode(int32_t type, int32_t transitionID, struct TraversalID *traversalID,
                                   struct TreeNode *treeNodeX, struct TreeNode *treeNodeY, float *wV) {
    struct TreeNode *treeNode = st_malloc(sizeof(struct TreeNode));
    treeNode->left = NULL;
    treeNode->refCount = 0;

    treeNode->type = type;
    treeNode->transitionID = transitionID;
    treeNode->traversalID = traversalID;
    treeNode->treeNodeX = treeNodeX;
    treeNode->treeNodeY = treeNodeY;
    treeNode->wV = NULL;
    if (treeNodeX != NULL) {
        treeNodeX->refCount++;
    }
    if (treeNodeY != NULL) {
        treeNodeY->refCount++;
    }
    treeNode->wV = wV ? memcpy(st_malloc(ALPHABET_SIZE*sizeof(float)), wV, ALPHABET_SIZE*sizeof(float)) : NULL;
    return treeNode;
}

void destructTreeNode(struct TreeNode *treeNode) {
    assert(treeNode->refCount >= 0);
    if(treeNode->refCount == 0) {
        if (treeNode->left != NULL) {
            treeNode->left->refCount--;
            destructTreeNode(treeNode->left);
        }
        if (treeNode->treeNodeX != NULL) {
            treeNode->treeNodeX->refCount--;
            destructTreeNode(treeNode->treeNodeX);
        }
        if (treeNode->treeNodeY != NULL) {
            treeNode->treeNodeY->refCount--;
            destructTreeNode(treeNode->treeNodeY);
        }
        if(treeNode->wV != NULL) {
            free(treeNode->wV);
        }
        free(treeNode);
    }
}

//edges
struct Edge *copyConstructEdge(int32_t from, int32_t to, float edgeScore, float subScore,
                               float *wV, int32_t silent, void *treeNode, int32_t iD) {
    struct Edge *temp;
    //INT_32 i;
    temp = st_malloc(sizeof(struct Edge));
    temp->from = from;
    temp->to = to;
    temp->edgeScore = edgeScore;
    temp->subScore = subScore;
    //temp->wV = wV;
    temp->wV = wV ? memcpy(st_malloc(ALPHABET_SIZE*sizeof(float)), wV, ALPHABET_SIZE*sizeof(float)) : NULL;
    temp->silent = silent;
    temp->treeNode = treeNode;
    temp->iD = iD;

    return temp;
}

void destructEdge(struct Edge *edge) {
    if(edge->wV)
        free(edge->wV);
    //free(edge->treeNode);
    destructTreeNode(edge->treeNode);
    free(edge);
}

//struct TraceBackEdge *constructTraceBackEdge(int64_t from, int64_t to, float edgeScore, struct Edge *edgeX, struct Edge *edgeY, char silent, void *getTreeNode) {
struct TraceBackEdge *constructTraceBackEdge(int64_t from, int64_t to, float edgeScore, struct Edge *edgeX, struct Edge *edgeY, char silent,
		struct TreeNode *(*getTreeNode)(struct AlignmentDataStructures *, struct TraceBackEdge *, int32_t)) {
//struct TraceBackEdge *constructTraceBackEdge(LONG_64 from, LONG_64 to, FLOAT_32 edgeScore, struct Edge *edgeX, struct Edge *edgeY) {
    struct TraceBackEdge *temp;

    temp = st_malloc(sizeof(struct TraceBackEdge));
    temp->from = from;
    temp->to = to;
    temp->edgeScore = edgeScore;
    temp->edgeX = edgeX;
    temp->edgeY = edgeY;
    temp->silent = silent;
    temp->getTreeNode = getTreeNode;
    return temp;
}

void destructTraceBackEdge(struct TraceBackEdge *edge) {
    free(edge);
}

int32_t edgeComparator(struct Edge *edge1, struct Edge *edge2) {
    int32_t i;
    int32_t j;
    float k;
    float l;

    i = edge1->to;
    j = edge2->to;
    if (i < j) {
        return -1;
    }
    if (i > j) {
        return 1;
    }
    i = edge1->from;
    j = edge2->from;
    if (i < j) {
        return -1;
    }
    if (i > j) {
        return 1;
    }
    if (edge1 == edge2) {
        return 0;
    }
    k = edge1->edgeScore;
    l = edge2->edgeScore;
    if (k < l) {
        return -1;
    }
    if (k > l) {
        return 1;
    }
    i = edge1->silent;
    j = edge2->silent;
    if (i < j) {
        return -1;
    }
    if (i > j) {
        return 1;
    }
    return edge1->treeNode < edge2->treeNode ? -1 : edge1->treeNode > edge2->treeNode ? 1 : 0;
    //return 0; //a possible but bizarre occurrence!
}

int32_t edgeComparatorStub(struct AlignmentDataStructures *aDS, struct Edge *edge1, struct Edge *edge2) {
	assert(aDS != NULL);
    return edgeComparator(edge1, edge2);
}

int edgePointerComparator(struct Edge **edge1, struct Edge **edge2) {
    return edgeComparator(*edge1, *edge2);
}

void sortEdges(void **list, int32_t length) {
    //edgeComparator(struct Edge *edge1, struct Edge *edge2)
    qsort(list, length, sizeof(void *), (int (*)(const void *, const void*))edgePointerComparator);
}

//sequence graphs
struct SequenceGraph *constructSequenceGraph(struct List *edges, int32_t vertexNo) {
    int32_t i;
    struct Edge *edge;
    struct SequenceGraph *sequenceGraph;
    int32_t *counts;

    sequenceGraph = st_malloc(sizeof(struct SequenceGraph));
    sequenceGraph->edges = edges;
    sequenceGraph->vertexNo = vertexNo;
    sequenceGraph->edgesArrangedByToVertex = st_malloc(sizeof(struct List *)*vertexNo);
    sequenceGraph->edgesArrangedByFromVertex = st_malloc(sizeof(struct List *)*vertexNo);
    counts = (int32_t *)st_malloc(sizeof(int32_t)*vertexNo);
    //tos
    for(i=0; i<vertexNo; i++) {
        counts[i] = 0;
    }
    //for edge in edges:
    for(i=0; i<edges->length; i++) {
        edge = edges->list[i];
        counts[edge->to]++;
    }
    //for edge in edges:
    for(i=0; i<vertexNo; i++) {
        sequenceGraph->edgesArrangedByToVertex[i] = constructEmptyList(counts[i], NULL); //assumes that will be destructed in edges list
        counts[i] = 0;
    }
    for(i=0; i<edges->length; i++) {
        edge = edges->list[i];
        sequenceGraph->edgesArrangedByToVertex[edge->to]->list[counts[edge->to]++] = edge;
    }
    //froms
    for(i=0; i<vertexNo; i++) {
        counts[i] = 0;
    }
    //for edge in edges:
    for(i=0; i<edges->length; i++) {
        edge = edges->list[i];
        counts[edge->from]++;
    }
    //for edge in edges:
    for(i=0; i<vertexNo; i++) {
        sequenceGraph->edgesArrangedByFromVertex[i] = constructEmptyList(counts[i], NULL);
        counts[i] = 0;
    }
    for(i=0; i<edges->length; i++) {
        edge = edges->list[i];
        sequenceGraph->edgesArrangedByFromVertex[edge->from]->list[counts[edge->from]++] = edge;
    }
    //for(i=0; i<vertexNo; i++) { not necessary
    //    sortEdges(sequenceGraph->edgesArrangedByFromVertex[i]->list, sequenceGraph->edgesArrangedByFromVertex[i]->length);
    //    sortEdges(sequenceGraph->edgesArrangedByToVertex[i]->list, sequenceGraph->edgesArrangedByToVertex[i]->length);
    //}
    free(counts);
    return sequenceGraph;
}

void destructSequenceGraph(struct SequenceGraph *sequenceGraph, int32_t freeIncomingEdges) {
    int32_t i;

    if(freeIncomingEdges) {
        destructList(sequenceGraph->edges);
    }
    for (i = 0; i < sequenceGraph->vertexNo; ++i) {
        destructList(sequenceGraph->edgesArrangedByFromVertex[i]);
        destructList(sequenceGraph->edgesArrangedByToVertex[i]);
    }
    free(sequenceGraph->edgesArrangedByFromVertex);
    free(sequenceGraph->edgesArrangedByToVertex);
    free(sequenceGraph);
}

int32_t *randomChoices(float *probs, int32_t sizeA, int32_t pathWeight) {
    //method for choosing a number of items according to their given probs
    static int32_t *choices;
    static int32_t choicesSize;
    float randomChoice;
    float totalProb;
    float cumulativeProbs[sizeA];
    int32_t i, j;

    choices = arrayResize(choices, &choicesSize, sizeA, sizeof(int32_t));
    totalProb = probs[0];
    cumulativeProbs[0] = totalProb;
    choices[0] = 0;
    for(i=1; i< sizeA; i++) {
    	choices[i] = 0;
        LOG_PLUS_EQUALS(&totalProb, probs[i]);
        cumulativeProbs[i] = totalProb;
    }

    assert(totalProb != LOG_ZERO);

    //choiceProbs = [ LOG(RANDOM()) + totalProb for i in xrange(0, pathWeight) ]
    for(i=0; i< pathWeight; ++i) {
    	randomChoice = RANDOM_LOG() + totalProb;
    	for (j = 0; j < sizeA; ++j) {
			if (cumulativeProbs[j] >= randomChoice) {
				choices[j]++;
				break;
			}
		}
    }

    return choices;
}

//graph member holders

struct GraphMemberHolder *constructGraphMember(void *graphMember, int32_t *sequenceConstraints, void (*destructGraphMember)(void *)) {
	struct GraphMemberHolder *graphMemberHolder;
    graphMemberHolder = st_malloc(sizeof(struct GraphMemberHolder));
	graphMemberHolder->graphMember = graphMember;
	graphMemberHolder->sequenceConstraints = sequenceConstraints;
    graphMemberHolder->destructGraphMember = destructGraphMember;
	return graphMemberHolder;
}

void destructGraphMember(struct GraphMemberHolder *graphMemberHolder) {
    //don't destruct sequence constraints, as already handled
    if(graphMemberHolder->destructGraphMember != NULL) {
        graphMemberHolder->destructGraphMember(graphMemberHolder->graphMember);
    }
    free(graphMemberHolder);
}

int32_t getIDFromGraphMemberEdge(void *graphMember) {
    return ((struct Edge *)graphMember)->iD;
}

int32_t getIDfromGraphMemberVertex(void *graphMember) {
    return *((int32_t *)graphMember);
}

void computeUnionScratch(struct AlignmentDataStructures *aDS, struct List *list1, struct List *list2,
                        struct List *scratchList,
                        int32_t (*compare)(struct AlignmentDataStructures *, void *, void *)) {
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t l = 0;
    void *graphMember1;
    void *graphMember2;

    int32_t length1 = list1->length;
    int32_t length2 = list2->length;

    void **array1 = list1->list;
    void **array2 = list2->list;
    void **list3;

    listResize(scratchList, length1 + length2);
    list3 = scratchList->list;

    assert(list1->destructElement == list2->destructElement);
    while (i < length1) {
        if (j < length2) {
            graphMember1 = array1[i];
            graphMember2 = array2[j];
            k = (*compare)(aDS, graphMember1, graphMember2);
            if (k < 0) {
                list3[l] = graphMember1;
                i += 1;
            }
            else if (k > 0) {
                list3[l] = graphMember2;
                j += 1;
            }
            else {
                list3[l] = graphMember1;
                i += 1;
                j += 1;
            }
            l += 1;
        }
        else {
            k = (length1-i);
            memcpy(list3 + l, array1 + i, k*sizeof(void *));
            scratchList->length = l + k;
            return;
        }
    }
    k = (length2-j);
    memcpy(list3 + l, array2 + j, k*sizeof(void *));
    scratchList->length = l + k;
}

struct List *computeUnion(struct AlignmentDataStructures *aDS, struct List *list1, struct List *list2,
                   int32_t (*compare)(struct AlignmentDataStructures *, void *, void *)) {
    static struct List scratchList;

    computeUnionScratch(aDS, list1, list2, &scratchList, compare);
    return copyConstructList(scratchList.list, scratchList.length, list1->destructElement);
}

void diffScratch(struct AlignmentDataStructures *aDS, struct List *list1, struct List *list2,
                struct List *scratchList,
                int32_t (*compare)(struct AlignmentDataStructures *, void *, void *)) {
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t l = 0;
    void *graphMember1;
    void *graphMember2;

    static void **list3;

    listResize(scratchList, list1->length);
    list3 = scratchList->list;

    assert(list1->destructElement == list2->destructElement);
    while (i < list1->length) {
        graphMember1 = list1->list[i];
        if (j < list2->length) {
            graphMember2 = list2->list[j];
            k = (*compare)(aDS, graphMember1, graphMember2);
            if (k < 0) {
                list3[l] = graphMember1;
                l += 1;
                i += 1;
            }
            else if (k > 0) {
                j += 1;
            }
            else {
                i += 1;
                while (i < list1->length) {
                    graphMember1 = list1->list[i];
                    k = (*compare)(aDS, graphMember1, graphMember2);
                    if (k > 0) {
                        break;
                    }
                    assert(k == 0);
                    i += 1;
                }
                j += 1;
            }
        }
        else {
            k = (list1->length - i);
            memcpy(list3 + l, list1->list + i, k*sizeof(void *));
            scratchList->length = l + k;
            return;
        }
    }
    scratchList->length = l;
}

void filterScratch(struct AlignmentDataStructures *aDS, struct List *list,
                  struct List *scratchList,
                  int32_t (*filterFn)(struct AlignmentDataStructures *, void *)) {
    int32_t i = 0;
    int32_t j = 0;

    static void **list2;

    void **array = list->list;
    int32_t length = list->length;

    listResize(scratchList, length);
    list2 = scratchList->list;

    while (i < length) {
        if ((*filterFn)(aDS, array[i])) {
            list2[j] = array[i];
            j += 1;
        }
        i += 1;
    }
    scratchList->length = j;
}

void divideScratch(struct AlignmentDataStructures *aDS, struct List *list,
                   struct List *scratchList1, struct List *scratchList2,
                   int32_t (*filterFn)(struct AlignmentDataStructures *, void *)) {
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;

    static void **list2;
    static void **list3;

    listResize(scratchList1, list->length);
    list2 = scratchList1->list;

    listResize(scratchList2, list->length);
    list3 = scratchList2->list;

    while (i < list->length) {
        if ((*filterFn)(aDS, list->list[i])) {
            list2[j++] = list->list[i];
        }
        else {
            list3[k++] = list->list[i];
        }
        i += 1;
    }
    scratchList1->length = j;
    scratchList2->length = k;
}

int32_t graphMember_EdgeComparator(struct AlignmentDataStructures *aDS, struct GraphMemberHolder *graphMemberHolder1,
                               struct GraphMemberHolder *graphMemberHolder2) {
    int32_t i;

    i = edgeComparatorStub(aDS, graphMemberHolder1->graphMember, graphMemberHolder2->graphMember);
    if (i == 0) {
        return intsComparator(graphMemberHolder1->sequenceConstraints,
                              graphMemberHolder2->sequenceConstraints, aDS->leafSeqNoX);
    }
    return i;
}

int32_t graphMember_VertexComparator(struct AlignmentDataStructures *aDS, struct GraphMemberHolder *graphMemberHolder1,
                                 struct GraphMemberHolder *graphMemberHolder2) {
    int32_t i;

    i = intComparator(graphMemberHolder1->graphMember, graphMemberHolder2->graphMember);
    if (i == 0) {
        return intsComparator(graphMemberHolder1->sequenceConstraints, graphMemberHolder2->sequenceConstraints, aDS->leafSeqNoX);
    }
    return i;
}

int32_t edge_graphMember_EdgeComparator(struct AlignmentDataStructures *aDS, struct Edge *edge,
                                    struct GraphMemberHolder *graphMemberHolder) {
    return edgeComparatorStub(aDS, edge, graphMemberHolder->graphMember);
}

int32_t edgeVertexFilter(struct AlignmentDataStructures *aDS, struct Edge *edge) {
    int32_t *i;

    i = bsearch(&(edge->from), aDS->newVertices, aDS->noOfNewVertices, sizeof(int32_t),
                (int (*)(const void *, const void *))intComparator_Int);
    return i != NULL;
}

int32_t greaterThanConditionFilter(struct AlignmentDataStructures *aDS,
                               struct GraphMemberHolder *graphMemberHolder) {
    int32_t i;
    int32_t seq;
    int32_t seqCoordinate;

    //for seq, seqCoordinate in aDS.changes
    for (i=0; i<aDS->noOfChanges; i+=2) {
        seq = aDS->changes[i];
        seqCoordinate = aDS->changes[i+1];
        if (seqCoordinate >= graphMemberHolder->sequenceConstraints[seq]) {
            return FALSE;
        }
    }
    return TRUE;
}

int32_t lessThanConditionFilter(struct AlignmentDataStructures *aDS,
                            struct GraphMemberHolder *graphMemberHolder) {
    int32_t i;
    int32_t j;
    int32_t k;
    struct List *sequenceConstraintsCollection;
    int32_t *sequenceConstraints;
    struct List *sequenceCoordinatesCollection;
    int32_t *sequenceCoordinates;

    sequenceConstraintsCollection = aDS->mergedEndConstraints_GraphMembers[aDS->getIDFromGraphMember(graphMemberHolder->graphMember)];
    sequenceCoordinatesCollection = aDS->sequenceCoordinatesCollection;
    //for sequenceConstraints in aDS.mergedEndConstraints_GraphMembers[graphMember.graphMember]:
    for (i=0; i<sequenceConstraintsCollection->length; i++) {
        sequenceConstraints = sequenceConstraintsCollection->list[i];
        //for sequenceCoordinates in aDS.sequenceCoordinatesCollection
        j=0;
        outer:
        while (j<sequenceCoordinatesCollection->length) {
            sequenceCoordinates = sequenceCoordinatesCollection->list[j];
            //for seq in xrange(0, aDS.leafSeqNoX):
            for (k=0; k<aDS->leafSeqNoX; k++) {
                if (sequenceCoordinates[k] <= sequenceConstraints[k]) {
                    j++;
                    goto outer;
                }
            }
            return TRUE;
        }
    }
    return FALSE;
}

void deleteEndsScratch(struct AlignmentDataStructures *aDS, struct List *previousVertices, struct List *previousEdges,
                struct List *previousIllegalEdges, struct List *sequenceCoordinatesCollection,
                struct List *previousSequenceCoordinatesCollection,
                struct List *newVertices, struct List *newEdges, struct List *newIllegalEdges) {
    int32_t i;
    static struct List scratchList1;
    //Deletes edges that were previously legal, but are not longer.
    //if (previousSequenceCoordinatesCollection == sequenceCoordinatesCollection)
    //    return previousVertices, previousEdges, previousIllegalEdges
    if(sequenceCoordinatesCollection->length == previousSequenceCoordinatesCollection->length) {
        for (i=0; i<sequenceCoordinatesCollection->length; i++) {
            if(memcmp(sequenceCoordinatesCollection->list[i], previousSequenceCoordinatesCollection->list[i], aDS->leafSeqNoX*sizeof(int32_t))) {
                goto outer;
            }
        }
        copyList(previousVertices, newVertices);
        copyList(previousEdges, newEdges);
        copyList(previousIllegalEdges, newIllegalEdges);
        return;
    }
    outer:
    aDS->sequenceCoordinatesCollection = sequenceCoordinatesCollection;
    aDS->mergedEndConstraints_GraphMembers = aDS->mergedEndConstraints_Edges;
    aDS->getIDFromGraphMember = getIDFromGraphMemberEdge;

    divideScratch(aDS, previousEdges, newEdges, newIllegalEdges, (int32_t (*)(struct AlignmentDataStructures *, void *))lessThanConditionFilter);
    if(newIllegalEdges->length == 0) {
        copyList(previousVertices, newVertices);
        copyList(previousIllegalEdges, newIllegalEdges);
        return;
    }
    aDS->mergedEndConstraints_GraphMembers = aDS->mergedEndConstraints_Vertices;
    aDS->getIDFromGraphMember = getIDfromGraphMemberVertex;

    filterScratch(aDS, previousVertices, newVertices, (int32_t (*)(struct AlignmentDataStructures *, void *))lessThanConditionFilter);
    //dealing with illegal edges
    //newIllegalEdges = [ i.graphMember for i in newIllegalEdges ]
    for (i = 0; i < newIllegalEdges->length; ++i) {
        newIllegalEdges->list[i] = ((struct GraphMemberHolder *)newIllegalEdges->list[i])->graphMember;
    }
    computeUnionScratch(aDS, newIllegalEdges, previousIllegalEdges, &scratchList1,
                       (int32_t (*)(struct AlignmentDataStructures *, void *, void *))edgeComparatorStub);
    //aDS.newVertices = [ i.graphMember for i in newVertices ]
    aDS->newVertices = arrayResize(aDS->newVertices, &aDS->newVertices_Size, newVertices->length, sizeof(int32_t));
    for (i = 0; i < newVertices->length; ++i) {
        aDS->newVertices[i] = *((int32_t *)((struct GraphMemberHolder *)newVertices->list[i])->graphMember);
    }
    aDS->noOfNewVertices = newVertices->length;
    filterScratch(aDS, &scratchList1, newIllegalEdges, (int32_t (*)(struct AlignmentDataStructures *, void *))edgeVertexFilter);
}

struct Edge *getMinEdge(void** edges, int32_t *size) { //not the best way of doing this, but still
    int32_t i;
    int32_t j;
    int32_t k;
    struct Edge *edge1;
    struct Edge *edge2;

    edge1 = edges[0];
    k = 0;
    for(i=1; i<(*size); i++) {
        edge2 = edges[i];
        j = edgeComparator(edge1, edge2);
        if(j > 0) {
            k = i;
            edge1 = edge2;
        }
    }
    edges[k] = edges[--(*size)];
    return edge1;
}

int32_t extendEnds_isLegalEnds(struct AlignmentDataStructures *aDS, int32_t *sequenceCoordinates, int32_t *sequenceConstraints) {
    int32_t seq;
    for (seq=0; seq<aDS->leafSeqNoX; seq++) {
        if (sequenceCoordinates[seq] <= sequenceConstraints[seq]) {
            return FALSE;
        }
    }
    return TRUE;
}

void extendEnds(struct AlignmentDataStructures *aDS, struct List *vertices, struct List *edges, struct List *illegalEdges,
                struct List *sequenceCoordinatesCollection_Ends,
                struct List **newVertices, struct List **newEdges, struct List **newIllegalEdges) {
    //Extends the legal vertices along the y-sequence.
    int32_t i;
    int32_t j;
    int32_t k;
    int32_t l;
    int32_t hashDummy = 0;

    struct List temp;

    struct hashtable *newVertices_Set;

    static void **newVertices_List;
    static int32_t newVertices_ListSize;
    int32_t newVertices_List_Index = 0;

    static void **newEdges_List;
    static int32_t newEdges_ListSize;
    int32_t newEdges_List_Index = 0;

    static void **newIllegalEdges_List;
    static int32_t newIllegalEdges_ListSize;
    int32_t newIllegalEdges_List_Index = 0;

    static void **illegalEdges_List;
    static int32_t illegalEdges_ListSize;
    int32_t illegalEdges_List_Index = 0;

    void *previous;
    struct Edge *edge;

    struct List *sequenceConstraintsCollection;
    int32_t *sequenceConstraints;

    struct List *startSequenceConstraintsCollection;
    //set up vertex hash
    newVertices_Set = create_hashtable(vertices->length + SMALL_CHUNK_SIZE, hashtable_intHashKey, hashtable_intEqualKey, NULL, NULL);
    for (i = 0; i < vertices->length; ++i) {
        hashtable_insert(newVertices_Set, ((struct GraphMemberHolder *)vertices->list[i])->graphMember, &hashDummy); //value is dummy
    }

    //set up illegal edges
    illegalEdges_List = arrayResize(illegalEdges_List, &illegalEdges_ListSize, illegalEdges->length, sizeof(void *));
    for (i = 0; i < illegalEdges->length; ++i) {
        illegalEdges_List[i] = illegalEdges->list[i];
    }
    illegalEdges_List_Index = illegalEdges->length;
    previous = NULL;
    outer:
    while (illegalEdges_List_Index > 0) {
        edge = getMinEdge(illegalEdges_List, &illegalEdges_List_Index);
        if (edge == previous || hashtable_search(newVertices_Set, &edge->from) == NULL) { //the
            continue;
        }
        previous = edge;
        //for edgeConstraints in aDS.mergedEndConstraints_Edges[edge]:
        sequenceConstraintsCollection = aDS->mergedEndConstraints_Edges[edge->iD];
        for (i=0; i<sequenceConstraintsCollection->length; i++) {
            sequenceConstraints = sequenceConstraintsCollection->list[i];
            //for sequenceCoordinates in sequenceCoordinatesCollection_Ends:
            for (j=0; j<sequenceCoordinatesCollection_Ends->length; j++) {
                if (extendEnds_isLegalEnds(aDS, sequenceCoordinatesCollection_Ends->list[j], sequenceConstraints)) {
                    //is legal
                    //for startSequenceConstraints in aDS.mergedStartConstraints_Edges[edge]:
                    startSequenceConstraintsCollection = aDS->mergedStartConstraints_Edges[edge->iD];
                    newEdges_List = arrayPrepareAppend(newEdges_List, &newEdges_ListSize, newEdges_List_Index + startSequenceConstraintsCollection->length, sizeof(void *));
                    l = newEdges_List_Index;
                    for (k=0; k<startSequenceConstraintsCollection->length; k++) {
                        newEdges_List[newEdges_List_Index++] = startSequenceConstraintsCollection->list[k];
                    }
                    if(newEdges_List_Index > l) {
                        if (hashtable_search(newVertices_Set, &edge->to) == NULL) {
                            hashtable_insert(newVertices_Set, &edge->to, &hashDummy);
                            //for startSequenceConstraints in aDS.mergedStartConstraints_Vertices[edge.to]:
                            startSequenceConstraintsCollection = aDS->mergedStartConstraints_Vertices[edge->to];
                            newVertices_List = arrayPrepareAppend(newVertices_List, &newVertices_ListSize, newVertices_List_Index + startSequenceConstraintsCollection->length, sizeof(void *));
                            for (k=0; k<startSequenceConstraintsCollection->length; k++) {
                                newVertices_List[newVertices_List_Index++] = startSequenceConstraintsCollection->list[k];
                            }
                            //for edge2 in aDS.sequenceGraphY.edgesArrangedByFromVertex[edge.to]:
                            illegalEdges_List = arrayPrepareAppend(illegalEdges_List, &illegalEdges_ListSize, illegalEdges_List_Index + aDS->sequenceGraphY->edgesArrangedByFromVertex[edge->to]->length, sizeof(void *));
                            for (k=0; k<aDS->sequenceGraphY->edgesArrangedByFromVertex[edge->to]->length; k++) {
                                illegalEdges_List[illegalEdges_List_Index++] = aDS->sequenceGraphY->edgesArrangedByFromVertex[edge->to]->list[k];
                            }
                        }
                    }
                    goto outer;
                }
            }
        }
        newIllegalEdges_List = arrayPrepareAppend(newIllegalEdges_List, &newIllegalEdges_ListSize, newIllegalEdges_List_Index+1, sizeof(void *));
        newIllegalEdges_List[newIllegalEdges_List_Index++] = edge;
    }
    temp.list = newVertices_List; temp.length = newVertices_List_Index; temp.destructElement = NULL;
    (*newVertices) = computeUnion(aDS, &temp, vertices, (int32_t (*)(struct AlignmentDataStructures *, void *, void *))graphMember_VertexComparator);

    temp.list = newEdges_List; temp.length = newEdges_List_Index;
    (*newEdges)  = computeUnion(aDS, &temp, edges, (int32_t (*)(struct AlignmentDataStructures *, void *, void *))graphMember_EdgeComparator);

    (*newIllegalEdges) = copyConstructList(newIllegalEdges_List, newIllegalEdges_List_Index, NULL); //(void (*)(void *))destructEdge);
    //final memory clean up
    hashtable_destroy(newVertices_Set, FALSE, FALSE);
}

int64_t v(struct AlignmentDataStructures *aDS, int64_t x, int64_t y) {
    return (y + x*aDS->vertexYNo)*stateNo();
}


#define MAX_Z 10000

int64_t vZ(struct AlignmentDataStructures *aDS, int64_t x, int64_t y, int64_t z) {
    return (y + x*aDS->vertexYNo)*(stateNo()*MAX_Z) + (MAX_Z - z - 1)*stateNo();
}

void rVZ(struct AlignmentDataStructures *aDS, int64_t c, int32_t *x, int32_t *y, int32_t *z, int32_t *state) { //reverse vertex number function
    int32_t vertexNo;

    vertexNo = aDS->vertexYNo;
    *state = c % stateNo();
    *z = MAX_Z - (c % (MAX_Z*stateNo()))/stateNo() - 1;
    c /= MAX_Z*stateNo();
    *y = c % vertexNo;
    *x = c / vertexNo;
}

int32_t rVZ_State(struct AlignmentDataStructures *aDS, int64_t c) {
	assert(aDS != NULL);
    return c % stateNo();
}

int64_t rVZ_StatelessVertex(struct AlignmentDataStructures *aDS, int64_t c) {
    return c - rVZ_State(aDS, c);
}

int64_t rVZ_Zless(struct AlignmentDataStructures *aDS, int64_t c) {
	int32_t x, y, z, state;
	rVZ(aDS, c, &x, &y, &z, &state);
	return v(aDS, x, y);
}

float *getTempCell(struct AlignmentDataStructures *aDS, int64_t vertex, struct Chunks *matrixChunks, struct hashtable *matrix) {
	assert(aDS != NULL);
    int64_t *chunk;
    float *cell;
    int32_t i;

    if ((cell = hashtable_search(matrix, &vertex)) == NULL) {
        chunk = mallocChunk(matrixChunks);
        *chunk = vertex;
        cell = (float *)(chunk + 1);
        for(i=0; i<stateNo(); i++) {
            cell[i] = LOG_ZERO;
        }
        hashtable_insert(matrix, chunk, cell);
    }
    return cell;
}

float *getCell(struct AlignmentDataStructures *aDS, int64_t vertex) {
    int64_t *chunk;
    struct List *chunkList;
    int32_t i;

    //search in the binary list
    chunkList = aDS->matrixChunks->chunkList;
    chunk = chunkList->list[chunkList->length-1];
    if(vertex >= *chunk) {
        chunk = bsearch(&vertex, chunk, aDS->matrixChunks->chunkSize - aDS->matrixChunks->remaining, aDS->matrixChunks->elementSize, (int (*)(const void *, const void *))longComparator_Int);
    }
    else {
        for(i=chunkList->length-2; i>=0; i--) {
            chunk = chunkList->list[i];
            if(vertex >= *chunk) {
                break;
            }
        }
        chunk = bsearch(&vertex, chunk, aDS->matrixChunks->chunkSize, aDS->matrixChunks->elementSize, (int (*)(const void *, const void *))longComparator_Int);
    }
    if(chunk == NULL) {
        return NULL;
    }
    return (float *)(chunk + 1);
}

void addCell(struct AlignmentDataStructures *aDS, int64_t vertex, float *cell, int64_t *pVertex) {
    assert(*pVertex < vertex);
    *pVertex = vertex;
    int64_t *chunk;
    float *cell2;

    chunk = mallocChunk(aDS->matrixChunks);
    *chunk = vertex;
    cell2 = (float *)(chunk + 1);
    memcpy(cell2, cell, stateNo()*sizeof(float));
}

void assignFromTo_Sum(struct AlignmentDataStructures *aDS, int32_t fromState, int32_t toState, float edgeScore) {
    LOG_PLUS_EQUALS(&aDS->toCell[toState], aDS->fromCell[fromState] + edgeScore);
}

int32_t getChanges(int32_t **vertexSequenceCoordinates, struct Edge *edge, int32_t *changes, int32_t leafSeqNo, int32_t leftMostSeqNo) {
    //function to calculate leaf sequence changes along edge
    int32_t i;
    int32_t *sequenceCoordinatesFrom;
    int32_t *sequenceCoordinatesTo;
    int32_t seq;
    static int32_t *activeLeaves;
    static int32_t activeLeavesSize;

    sequenceCoordinatesFrom = vertexSequenceCoordinates[edge->from];
    sequenceCoordinatesTo = vertexSequenceCoordinates[edge->to];

    activeLeaves = arrayResize(activeLeaves, &activeLeavesSize, leafSeqNo, sizeof(int32_t));
    memset(activeLeaves, FALSE, sizeof(int32_t)*leafSeqNo);
    getActiveLeaves(activeLeaves, leftMostSeqNo, edge->treeNode);
    //for seq in xrange(0, aDS.leafSeqNoX) {
    i = 0;
    for (seq=0; seq < leafSeqNo; seq++) {
        if ((sequenceCoordinatesTo[seq] - sequenceCoordinatesFrom[seq] > 0) && activeLeaves[seq]) {
            changes[i++] = seq;
            changes[i++] = sequenceCoordinatesTo[seq];
        }
    }
    //assert(i > 0); -- this is not true, as delete-delete node will affect this
    return i;
}

int32_t filterDuplicatesFunction(struct List *graphMembers, void **list) {
    int32_t i;
    int32_t j;
    void *previous;
    void *current;

    previous = NULL;
    j=0;
    //for i in xrange(0, len(graphMembers)) {
    for(i=0; i < graphMembers->length; i++) {
        current = ((struct GraphMemberHolder *)graphMembers->list[i])->graphMember;
        if (current != previous) {
            list[j++] = current;
            previous = current;
        }
    }
    return j;
}

void printTreeNode(struct TreeNode *treeNode) {
    if(treeNode == NULL) {
        return;
    }
    if(!(treeNode->type & TREE_NODE_EFFECTIVELY_SILENT)) {
        printTreeNode(treeNode->treeNodeX);
        printTreeNode(treeNode->treeNodeY);
    }
    if(treeNode->type == TREE_NODE_LEAF) {
        printf(" active leaf " INT_STRING " \n", treeNode->traversalID->leafNo);
    }
}

void debugSets(struct AlignmentDataStructures *aDS, struct List *newEdgesPointer, struct List *newVerticesPointer, struct List *illegalEdgesPointer) {
	if(DEBUG) { //debug code
        int32_t xx; //debug variables
        int32_t i;
        //printf(" the vertex is " INT_STRING ", of " INT_STRING " " INT_STRING " " INT_STRING " \n", toX, aDS->sequenceGraphX->vertexNo, aDS->sequenceGraphY->vertexNo, aDS->sequenceGraphXSilentVertices_To[toX]);
        //for(xx=0; xx<aDS->leafSeqNoX; xx++) { printf(" coordinate " INT_STRING " \n", ((INT_32*)aDS->vertexXSequenceCoordinates[toX])[xx]); }
        //for(xx=0; xx<newVerticesPointer->length; xx++) { INT_32 *yy = ((struct GraphMemberHolder *)newVerticesPointer->list[xx])->graphMember; printf(" 1351vert " INT_STRING " " , *yy); } printf("\n");
        //for(xx=0; xx<newEdgesPointer->length; xx++) { struct Edge *yy = ((struct GraphMemberHolder *)newEdgesPointer->list[xx])->graphMember; printf(" 1352edge " INT_STRING " " INT_STRING " " INT_STRING " " , yy->from, yy->to, yy->silent); } printf("\n");
        //for(xx=0; xx<illegalEdgesPointer->length; xx++) { struct Edge *yy = illegalEdgesPointer->list[xx]; printf(" 1353ill " INT_STRING " " INT_STRING " " INT_STRING " " , yy->from, yy->to, yy->silent); } printf("\n");  printf("\n");
        int32_t xxx = -100000;
        int32_t yyy = -100000;
        int32_t *prev = NULL;
        int32_t *current = NULL;
        struct Edge *pEdge = NULL;
        for(xx=0; xx<newEdgesPointer->length; xx++) {
            struct Edge *yy = (struct Edge *)((struct GraphMemberHolder *)newEdgesPointer->list[xx])->graphMember;
            assert(xxx <= yy->to);
            if (yy->to == xxx) {
                assert(yyy <= yy->from);
                if(edgeComparator(pEdge, yy) == 0) {
                    if(memcmp(prev, ((struct GraphMemberHolder *)newEdgesPointer->list[xx])->sequenceConstraints, aDS->leafSeqNoX*sizeof(int32_t)) == 0) {
                        fprintf(stderr, "New edges the same " INT_STRING " " INT_STRING " " INT_STRING " \n", yy->from, yy->to, yy->silent);
                        assert(FALSE);
                    }
                }
            }
            prev = ((struct GraphMemberHolder *)newEdgesPointer->list[xx])->sequenceConstraints;
            pEdge = yy;
            xxx = yy->to;
            yyy = yy->from;
        }
        xxx = -100000;
        for(xx=0; xx<newVerticesPointer->length; xx++) {
            int32_t yy = *((int32_t *)((struct GraphMemberHolder *)newVerticesPointer->list[xx])->graphMember);
            assert(xxx <= yy);
            if(xxx == yy) {
                if(memcmp(prev, ((struct GraphMemberHolder *)newVerticesPointer->list[xx])->sequenceConstraints, aDS->leafSeqNoX*sizeof(int32_t)) == 0) {
                    current = ((struct GraphMemberHolder *)newVerticesPointer->list[xx])->sequenceConstraints;
                    for(i=0; i<aDS->leafSeqNoX; i++) {
                        fprintf(stderr, " Failed " INT_STRING " " INT_STRING " \n ", prev[i], current[i]);
                    }
                    fprintf(stderr, " Failed, vertices the same : " INT_STRING "  \n", xxx);
                    assert(FALSE);
                }
            }
            prev = ((struct GraphMemberHolder *)newVerticesPointer->list[xx])->sequenceConstraints;
            xxx = yy;
        }
        //assert(newVerticesPointer->length > 0); (this will not be true when
        xxx = -100000;
        yyy = -100000;
        for(xx=0; xx<illegalEdgesPointer->length; xx++) {
            struct Edge *yy = (struct Edge *)illegalEdgesPointer->list[xx];
            assert(xxx <= yy->to);
            if (yy->to == xxx) {
                assert(yyy <= yy->from);
                assert(edgeComparator(pEdge, yy));
            }
            pEdge = yy;
            xxx = yy->to;
            yyy = yy->from;
        }
    }  //#end debug code
}

void computeMatrix(struct AlignmentDataStructures *aDS) {
    //memory that needs cleaning up
    void **previousVertices;
    void **previousEdges;
    void **previousIllegalEdges;
    int32_t *rightMostVertices;
    //end
    //core structures computed during each loop
    static struct List newVertices;
    static struct List illegalEdges;
    static struct List newEdges;

    static struct List filteredVertices;
    static struct List filteredIllegalEdges;
    static struct List filteredEdges;

    static struct List scratchList1;
    static struct List scratchList2;

    struct List *newVerticesPointer;
    struct List *newEdgesPointer;
    struct List *illegalEdgesPointer;

    static void **scratch;
    static int32_t scratchSize;
    static void **scratch2;
    static int32_t scratchSize2;

    struct Edge *edgeX;
    struct Edge *edgeY;

    struct hashtable *columnCellsHash;
    struct Chunks *columnCellsChunks;
    float* cell;

    int32_t i;
    int32_t j;
    int32_t k;
    int32_t l;
    int32_t toX;
    int32_t fromX;
    int32_t toY;
    int32_t vertexY;
    int32_t state;
    int32_t rightmostVertex;
    int64_t pVertex;

    scratch = arrayResize(scratch, &scratchSize, aDS->sequenceGraphY->edges->length + 1, sizeof(void *)); //can't have duplicates loaded on, so okay
    scratch2 = arrayResize(scratch2, &scratchSize2, aDS->sequenceGraphY->edges->length + 1, sizeof(void *));

    previousVertices = st_calloc(aDS->sequenceGraphX->vertexNo, sizeof(void *));
    previousEdges = st_calloc(aDS->sequenceGraphX->vertexNo, sizeof(void *));
    previousIllegalEdges = st_calloc(aDS->sequenceGraphX->vertexNo, sizeof(void *));
    rightMostVertices = st_calloc(aDS->sequenceGraphX->vertexNo, sizeof(int32_t));

    //core structures computed during each loop
    copyList(aDS->mergedStartConstraints_Vertices[0], &newVertices);
    copyList(aDS->sequenceGraphY->edgesArrangedByFromVertex[0], &illegalEdges);
    sortEdges(illegalEdges.list, illegalEdges.length);
    assert(newEdges.length == 0);

    pVertex = -1;
    columnCellsHash = create_hashtable(MEDIUM_CHUNK_SIZE, hashtable_longHashKey, hashtable_longEqualKey, NULL, NULL);
    columnCellsChunks = constructChunks(SMALL_CHUNK_SIZE, sizeof(int64_t) + sizeof(float)*stateNo());
    cell = getTempCell(aDS, 0, columnCellsChunks, columnCellsHash);
    if(aDS->treeStates == NULL || TRUE) {
        for (i=0; i<stateNo(); i++) {
            cell[i] = aDS->startStates[i];
        }
    }
    else {
        for (i=0; i<stateNo(); i++) {
            cell[i] = LOG_ZERO; //defensive (should be unnecessary)
        }
        //cell[0] = LOG_ONE;
        cell[aDS->treeStates[aDS->traversalID->mid]] = LOG_ONE;
    }
    //turnOnDeleteXYLoopCorrection(aDS->model);
    //for toX in xrange(0, aDS.sequenceGraphX.vertexNo):
    for (toX=0; toX < aDS->sequenceGraphX->vertexNo; toX++) {
        if (aDS->sequenceGraphXSilentVertices_To[toX]) {
            //horizontal incoming transitions
            //for edgeX in aDS.sequenceGraphX.edgesArrangedByToVertex[toX]: {
            for (i=0; i< aDS->sequenceGraphX->edgesArrangedByToVertex[toX]->length; i++) {
                edgeX = aDS->sequenceGraphX->edgesArrangedByToVertex[toX]->list[i];
                fromX = edgeX->from;
                deleteEndsScratch(aDS, previousVertices[fromX], previousEdges[fromX], previousIllegalEdges[fromX],
                                  aDS->mergedEndSequenceCoordinates_Edges[edgeX->iD],
                                  aDS->mergedEndSequenceCoordinates_Vertices[edgeX->from],
                                  &filteredVertices, &filteredEdges, &filteredIllegalEdges);
                k = filterDuplicatesFunction(&filteredVertices, scratch);
                for (j=0; j<k; j++) {
                    vertexY = *((int32_t *)scratch[j]);
                    aDS->fromCell = getCell(aDS, v(aDS, fromX, vertexY));
                    aDS->toCell = getTempCell(aDS, v(aDS, toX, vertexY), columnCellsChunks, columnCellsHash);
                    //for state in xrange(0, aDS.stateNo):
                    for (state=0; state<stateNo(); state++) {
                        assignFromTo_Sum(aDS, state, state, edgeX->edgeScore);
                    }
                }
                computeUnionScratch(aDS, &filteredVertices, &newVertices, &scratchList1, (int32_t (*)(struct AlignmentDataStructures *, void *, void *))graphMember_VertexComparator);
                swapListFields(&newVertices, &scratchList1);

                computeUnionScratch(aDS, &filteredEdges, &newEdges, &scratchList1, (int32_t (*)(struct AlignmentDataStructures *, void *, void *))graphMember_EdgeComparator);
                swapListFields(&newEdges, &scratchList1);

                computeUnionScratch(aDS, &filteredIllegalEdges, &illegalEdges, &scratchList1, (int32_t (*)(struct AlignmentDataStructures *, void *, void *))edgeComparatorStub);
                swapListFields(&illegalEdges, &scratchList1);
            }
            diffScratch(aDS, &illegalEdges, &newEdges, &scratchList1, (int32_t (*)(struct AlignmentDataStructures *, void *, void *))edge_graphMember_EdgeComparator);
            swapListFields(&illegalEdges, &scratchList1);

            newVerticesPointer = copyConstructList(newVertices.list, newVertices.length, NULL);
            newEdgesPointer = copyConstructList(newEdges.list, newEdges.length, NULL);
            illegalEdgesPointer = copyConstructList(illegalEdges.list, illegalEdges.length, NULL);
        }
        else {
            //for edgeX in aDS.sequenceGraphX.edgesArrangedByToVertex[toX] {
            for (i=0; i< aDS->sequenceGraphX->edgesArrangedByToVertex[toX]->length; i++) {
                edgeX = aDS->sequenceGraphX->edgesArrangedByToVertex[toX]->list[i];
                fromX = edgeX->from;
                //horizontal incoming transitions
                aDS->noOfChanges = getChanges(aDS->vertexXSequenceCoordinates, edgeX, aDS->changes, aDS->leafSeqNoX, aDS->leftMostSeqNoX);
                filterScratch(aDS, previousVertices[fromX], &scratchList1, (int32_t (*)(struct AlignmentDataStructures *, void *))greaterThanConditionFilter);
                filterScratch(aDS, previousEdges[fromX], &scratchList2, (int32_t (*)(struct AlignmentDataStructures *, void *))greaterThanConditionFilter);

                deleteEndsScratch(aDS, &scratchList1, &scratchList2, previousIllegalEdges[fromX],
                           aDS->mergedEndSequenceCoordinates_Edges[edgeX->iD],
                           aDS->mergedEndSequenceCoordinates_Vertices[edgeX->from],
                           &filteredVertices, &filteredEdges, &filteredIllegalEdges);
                //for vertexY in filterDuplicatesFunction(filteredVertices) {
                k = filterDuplicatesFunction(&filteredVertices, scratch);
                for (j=0; j<k; j++) {
                    vertexY = *((int32_t *)scratch[j]);
                    if (!aDS->sequenceGraphYSilentVertices_To[vertexY]) {
                        //x-gap
                        aDS->fromCell = getCell(aDS, v(aDS, fromX, vertexY));
                        aDS->toCell = getTempCell(aDS, v(aDS, toX, vertexY), columnCellsChunks, columnCellsHash);
                        insertXFn(aDS, aDS->model, edgeX, assignFromTo_Sum);
                        deleteYFn(aDS, aDS->model, edgeX, assignFromTo_Sum);
                    }
                }
                //matches
                //for edgeY in filterDuplicatesFunction(filteredEdges) {
                k = filterDuplicatesFunction(&filteredEdges, scratch);
                for (j=0; j<k; j++) {
                    edgeY = scratch[j];
                    if (!edgeY->silent) {
                        aDS->fromCell = getCell(aDS, v(aDS, fromX, edgeY->from));
                        aDS->toCell = getTempCell(aDS, v(aDS, toX, edgeY->to), columnCellsChunks, columnCellsHash);
                        matchFn(aDS, aDS->model, edgeX, edgeY, assignFromTo_Sum);
                    }
                }
                computeUnionScratch(aDS, &filteredVertices, &newVertices, &scratchList1, (int32_t (*)(struct AlignmentDataStructures *, void *, void *))graphMember_VertexComparator);
                swapListFields(&newVertices, &scratchList1);

                computeUnionScratch(aDS, &filteredEdges, &newEdges, &scratchList1, (int32_t (*)(struct AlignmentDataStructures *, void *, void *))graphMember_EdgeComparator);
                swapListFields(&newEdges, &scratchList1);

                computeUnionScratch(aDS, &filteredIllegalEdges, &illegalEdges, &scratchList1, (int32_t (*)(struct AlignmentDataStructures *, void *, void *))edgeComparatorStub);
                swapListFields(&illegalEdges, &scratchList1);
            }
            diffScratch(aDS, &illegalEdges, &newEdges, &scratchList1, (int32_t (*)(struct AlignmentDataStructures *, void *, void *))edge_graphMember_EdgeComparator);
            swapListFields(&illegalEdges, &scratchList1);

            extendEnds(aDS, &newVertices, &newEdges, &illegalEdges,
                       //aDS.mergedEndSequenceCoordinates_Vertices[toX],
                       aDS->mergedEndSequenceCoordinates_Vertices[toX],
                       &newVerticesPointer, &newEdgesPointer, &illegalEdgesPointer);
            i = filterDuplicatesFunction(newEdgesPointer, scratch);
            j = 0;
            //for toY in filterDuplicatesFunction(newVertices) {
            k = filterDuplicatesFunction(newVerticesPointer, scratch2);
            for (l=0; l<k; l++) {
                toY = *((int32_t *)scratch2[l]);
                aDS->toCell = getTempCell(aDS, v(aDS, toX, toY), columnCellsChunks, columnCellsHash);
                if (aDS->sequenceGraphYSilentVertices_To[toY]) {
                    while (j < i) {
                        edgeY = scratch[j];
                        assert(edgeY->to >= toY);
                        if (edgeY->to == toY) {
                            aDS->fromCell = getTempCell(aDS, v(aDS, toX, edgeY->from), columnCellsChunks, columnCellsHash);
                            //for state in xrange(0, aDS.stateNo):
                            for (state=0; state < stateNo(); state++) {
                                assignFromTo_Sum(aDS, state, state, edgeY->edgeScore);
                            }
                            j += 1;
                        } else {
                            break;
                        }
                    }
                }
                else {
                    while (j < i) {
                        edgeY = scratch[j];
                        assert(edgeY->to >= toY);
                        if (edgeY->to == toY) {
                            aDS->fromCell = getTempCell(aDS, v(aDS, toX, edgeY->from), columnCellsChunks, columnCellsHash);
                            //for state in xrange(0, aDS.stateNo):
                            insertYFn(aDS, aDS->model, edgeY, assignFromTo_Sum);
                            deleteXFn(aDS, aDS->model, edgeY, assignFromTo_Sum);
                            j += 1;
                        } else {
                            break;
                        }
                    }
                    aDS->fromCell = aDS->toCell;
                    silentFn(aDS, aDS->model, aDS->fromCell, assignFromTo_Sum);
                    //deleteFn(aDS, aDS->model, assignFromTo_Sum);
                }
            }
            assert(j == i);
        }
        k = filterDuplicatesFunction(newVerticesPointer, scratch2);
        for (l=0; l<k; l++) {
            toY = *((int32_t *)scratch2[l]);
            addCell(aDS, v(aDS, toX, toY), getTempCell(aDS, v(aDS, toX, toY), columnCellsChunks, columnCellsHash), &pVertex);
        }
        //clean up loop memory of previous start vertices if possible
        rightmostVertex = 0;//INT_MIN;
        //for edgeX in aDS.sequenceGraphX.edgesArrangedByFromVertex[toX] {
        for (i=0; i < aDS->sequenceGraphX->edgesArrangedByFromVertex[toX]->length; i++) {
            edgeX = aDS->sequenceGraphX->edgesArrangedByFromVertex[toX]->list[i];
            if (edgeX->to > rightmostVertex) {
                rightmostVertex = edgeX->to;
            }
        }
        //for fromVertex in Set([edgeX.from for edgeX in aDS.sequenceGraphX.edgesArrangedByToVertex[toX]]) {
        j = 0;
        sortEdges(aDS->sequenceGraphX->edgesArrangedByToVertex[toX]->list, aDS->sequenceGraphX->edgesArrangedByToVertex[toX]->length);
        for (i=0; i<aDS->sequenceGraphX->edgesArrangedByToVertex[toX]->length; i++) {
            edgeX = aDS->sequenceGraphX->edgesArrangedByToVertex[toX]->list[i];
            assert(edgeX->from >= j);
            if (edgeX->from > j) {
                if (rightMostVertices[edgeX->from] == toX) {
                    destructList(previousVertices[edgeX->from]);
                    previousVertices[edgeX->from] = NULL;
                    destructList(previousEdges[edgeX->from]);
                    previousEdges[edgeX->from] = NULL;
                    destructList(previousIllegalEdges[edgeX->from]);
                    previousIllegalEdges[edgeX->from] = NULL;
                }
                j = edgeX->from;
            }
        }
        debugSets(aDS, newEdgesPointer, newVerticesPointer, illegalEdgesPointer);
        //end loop memory clean up
        //for the next recursion
        previousVertices[toX] = newVerticesPointer;
        previousEdges[toX] = newEdgesPointer;
        previousIllegalEdges[toX] = illegalEdgesPointer;

        newVertices.length = 0;
        newEdges.length = 0;
        illegalEdges.length = 0;

        filteredVertices.length = 0;
        filteredEdges.length = 0;
        filteredIllegalEdges.length = 0;

        scratchList1.length = 0;
        scratchList2.length = 0;

        newVerticesPointer = NULL;
        newEdgesPointer = NULL;
        illegalEdgesPointer = NULL;

        destructChunks(columnCellsChunks);
        hashtable_destroy(columnCellsHash, FALSE, FALSE);
        columnCellsHash = create_hashtable(MEDIUM_CHUNK_SIZE, hashtable_longHashKey, hashtable_longEqualKey, NULL, NULL);
        columnCellsChunks = constructChunks(SMALL_CHUNK_SIZE, sizeof(int64_t) + sizeof(float)*stateNo());//constructEmptyList(0, free);
    }
    //memory clean up
    for(i=0; i<aDS->sequenceGraphX->vertexNo; i++) {
        if(previousVertices[i] != NULL) {
            assert(((struct List *)previousVertices[i])->destructElement == NULL);
            destructList(previousVertices[i]);
        }
        if(previousEdges[i] != NULL) {
            assert(((struct List *)previousEdges[i])->destructElement == NULL);
            destructList(previousEdges[i]);
        }
        if(previousIllegalEdges[i] != NULL) {
            assert(((struct List *)previousIllegalEdges[i])->destructElement == NULL);
            destructList(previousIllegalEdges[i]);
        }
    }
    free(previousVertices);
    free(rightMostVertices);
    free(previousEdges);
    free(previousIllegalEdges);
    destructChunks(columnCellsChunks);
    hashtable_destroy(columnCellsHash, FALSE, FALSE);
    //end memory clean up
    //end forward recursion
    j = 1;
    cell = getCell(aDS, v(aDS, aDS->sequenceGraphX->vertexNo-1, aDS->sequenceGraphY->vertexNo-1));
    if(cell != NULL) {
        st_logDebug("End state probabilities ");
        for (i = 0; i < stateNo(); ++i) {
            st_logDebug(" %f ", cell[i]);
            if (cell[i] + aDS->endStates[i] > LOG_ZERO) {
                j = 0;
            }
        }
        st_logDebug("\n");
    }
    else {
        fprintf(stderr, " No path through alignment possible, so I have no choice but to exit, sorry! \n");
        //logDebug(" No path through alignment possible, so I have no choice but to exit, sorry! \n");
        exit(73);
    }
    if(j) {
    	fprintf(stderr, " Zero prob, so I have no choice but to exit, sorry! \n");
    	//logDebug(" No path through alignment possible, so I have no choice but to exit, sorry! \n");
    	exit(73);
    }
}

void assignSampling(struct AlignmentDataStructures *aDS,
		int32_t fromState, int32_t toState, float edgeScore, float correctedEdgeScore) {
    float *fromVertices;
    struct TraceBackEdge *edge;

    if(aDS->potentialEdges_Index >= aDS->potentialEdges_Size) {
        //resize potential edges
        aDS->potentialEdges = arrayCopyResize(aDS->potentialEdges, &aDS->potentialEdges_Size, aDS->potentialEdges_Size*2, sizeof(void *));
        aDS->potentialEdges_Size /= 2;
        aDS->potentialEdgeCosts = arrayCopyResize(aDS->potentialEdgeCosts, &aDS->potentialEdges_Size, aDS->potentialEdges_Size*2, sizeof(float));
    }
    fromVertices = aDS->fromCell;
    aDS->potentialEdgeCosts[aDS->potentialEdges_Index] = fromVertices[fromState] + edgeScore;
    edge = constructTraceBackEdge(aDS->from + fromState, aDS->to + toState, correctedEdgeScore, aDS->edgeX, aDS->edgeY, aDS->silent, aDS->getTreeNode);
    aDS->potentialEdges[aDS->potentialEdges_Index++] = edge;
}

void assignSampling_CheckState(struct AlignmentDataStructures *aDS,
		int32_t fromState, int32_t toState, float edgeScore, float correctedEdgeScore) {
    //neccesary because model
    //may have multiple states for each type of transition
    if (aDS->state == toState) {
        assignSampling(aDS, fromState, toState, edgeScore, correctedEdgeScore);
    }
}

struct TreeNode *getTreeNode_insertX(struct AlignmentDataStructures *aDS, struct TraceBackEdge* traceBackEdge, int32_t transitionID) {
    return copyConstructTreeNode(TREE_NODE_INSERT, transitionID,
                                 aDS->traversalID, traceBackEdge->edgeX->treeNode, NULL, NULL);
}

struct TreeNode *getTreeNode_deleteX(struct AlignmentDataStructures *aDS, struct TraceBackEdge* traceBackEdge, int32_t transitionID, float *wV) {
	copyWV(traceBackEdge->edgeY->wV, wV, ALPHABET_SIZE);
	//multiplyWV(aDS->subModelX->deletionDistribution, traceBackEdge->edgeY->wV, wV);
    return copyConstructTreeNode(TREE_NODE_INTERNAL, transitionID,
                            aDS->traversalID, aDS->deleteNodeX, traceBackEdge->edgeY->treeNode, NULL);
}

struct TreeNode *getTreeNode_silentX(struct AlignmentDataStructures *aDS, struct TraceBackEdge* traceBackEdge, int32_t transitionID) {
    return copyConstructTreeNode(TREE_NODE_PREVIOUSLY_SILENT, transitionID,
                                 aDS->traversalID, traceBackEdge->edgeX->treeNode, NULL, NULL);
}

struct TreeNode *getTreeNode_insertY(struct AlignmentDataStructures *aDS, struct TraceBackEdge* traceBackEdge, int32_t transitionID) {
    return copyConstructTreeNode(TREE_NODE_INSERT, transitionID,
                                 aDS->traversalID, traceBackEdge->edgeY->treeNode, NULL, NULL);
}

struct TreeNode *getTreeNode_deleteY(struct AlignmentDataStructures *aDS, struct TraceBackEdge* traceBackEdge, int32_t transitionID, float *wV) {
    copyWV(traceBackEdge->edgeX->wV, wV, ALPHABET_SIZE);
	//multiplyWV(traceBackEdge->edgeX->wV, aDS->subModelY->deletionDistribution, wV);
    return copyConstructTreeNode(TREE_NODE_INTERNAL, transitionID,
                                 aDS->traversalID, traceBackEdge->edgeX->treeNode, aDS->deleteNodeY, NULL);
}

struct TreeNode *getTreeNode_silentY(struct AlignmentDataStructures *aDS, struct TraceBackEdge* traceBackEdge, int32_t transitionID) {
    return copyConstructTreeNode(TREE_NODE_PREVIOUSLY_SILENT, transitionID,
                              aDS->traversalID, traceBackEdge->edgeY->treeNode, NULL, NULL);
}

struct TreeNode *getTreeNode_matchXY(struct AlignmentDataStructures *aDS, struct TraceBackEdge* traceBackEdge, int32_t transitionID, float *wV) {
    multiplyWV(traceBackEdge->edgeX->wV, traceBackEdge->edgeY->wV, wV, ALPHABET_SIZE);
    normaliseWV(wV, wV, ALPHABET_SIZE);
    return copyConstructTreeNode(TREE_NODE_INTERNAL, transitionID,
                             aDS->traversalID, traceBackEdge->edgeX->treeNode, traceBackEdge->edgeY->treeNode, NULL);
}

struct TreeNode *getTreeNode_deleteXY(struct AlignmentDataStructures *aDS, struct TraceBackEdge* traceBackEdge, int32_t transitionID, float *wV) {
	assert(traceBackEdge != NULL);
	static float fA[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	copyWV(fA, wV, ALPHABET_SIZE);
	//multiplyWV(aDS->subModelX->deletionDistribution, aDS->subModelY->deletionDistribution, wV);
    return copyConstructTreeNode(TREE_NODE_INTERNAL, transitionID, aDS->traversalID, aDS->deleteNodeX, aDS->deleteNodeY, NULL);
}

struct TreeNode *getTreeNode_silentXY(struct AlignmentDataStructures *aDS, struct TraceBackEdge* traceBackEdge, int32_t transitionID) {
    assert(traceBackEdge != NULL);
	return copyConstructTreeNode(TREE_NODE_SILENT, transitionID, aDS->traversalID, NULL, NULL, NULL);
}

void updateIndices(struct TreeNode *treeNode, int32_t *currentIndices, int32_t leftMostLeafNo) {
    if (treeNode->left != NULL) { //go left
        updateIndices(treeNode->left, currentIndices, leftMostLeafNo);
    }
    if (treeNode->treeNodeX != NULL) {
        updateIndices(treeNode->treeNodeX, currentIndices, leftMostLeafNo);
    }
    if (treeNode->treeNodeY != NULL) {
        updateIndices(treeNode->treeNodeY, currentIndices, leftMostLeafNo);
    }
    if (treeNode->type == TREE_NODE_LEAF) {
        currentIndices[treeNode->traversalID->leafNo - leftMostLeafNo] += 1;
    }
}

int32_t **calculateVertexSequenceCoordinates(struct SequenceGraph *sequenceGraph, int32_t leftMostLeafNo, int32_t leafSeqNo, struct Chunks *chunks) {
    int32_t vertex;
    struct Edge *previousEdge;
    int32_t **sequenceCoordinates;
    int32_t *currentIndices;

    //calculates the coordinates for each vertex of the sequences in the graph
    //sequence coordinates are respresented in arrays offset from left-most leaf no
    sequenceCoordinates = st_malloc(sizeof(int32_t *)*sequenceGraph->vertexNo); //{ 0:[-1]*leafSeqNo }
    sequenceCoordinates[0] = mallocChunk(chunks); //callocLocal(leafSeqNo, sizeof(INT_32));
    memset(sequenceCoordinates[0], 0, sizeof(int32_t) * leafSeqNo);
    for(vertex=1; vertex<sequenceGraph->vertexNo; vertex++) { // in xrange(1, sequenceGraph.vertexNo):
        previousEdge = sequenceGraph->edgesArrangedByToVertex[vertex]->list[0];
        currentIndices = memcpy(mallocChunk(chunks), sequenceCoordinates[previousEdge->from], sizeof(int32_t)*leafSeqNo);
        updateIndices(previousEdge->treeNode, currentIndices, leftMostLeafNo);
        sequenceCoordinates[vertex] = currentIndices;
    }
    return sequenceCoordinates;
}

void compactVertex(struct List *toEdges, struct List *fromEdges, struct SequenceGraph *sequenceGraph) {
    struct Edge *edge;
    struct TreeNode *treeNode;
    struct Edge *edgeFrom;
    int32_t i;
    int32_t j;

    static struct List tempList;

    listResize(&tempList, toEdges->length * fromEdges->length);
    tempList.length = 0;

    for(i=0; i<toEdges->length; i++) {
        edge = toEdges->list[i];
        assert(edge->silent);
    }
    for(i=0; i<toEdges->length-1; i++) {
        edgeFrom = toEdges->list[i];
        for(j=0; j<fromEdges->length; j++) {
            edge = fromEdges->list[j];
            treeNode = edge->treeNode;
            treeNode = copyConstructTreeNode(treeNode->type, treeNode->transitionID, treeNode->traversalID, treeNode->treeNodeX, treeNode->treeNodeY, edge->treeNode->wV);
            treeNode->left = edgeFrom->treeNode;
            edgeFrom->treeNode->refCount++;
            listAppend(sequenceGraph->edgesArrangedByToVertex[edge->to], copyConstructEdge(edgeFrom->from, edge->to, edgeFrom->edgeScore + edge->edgeScore, edge->subScore,
                                                                                           edge->wV, edge->silent, treeNode, INT32_MAX));
        }
    }
    edgeFrom = toEdges->list[toEdges->length-1];
    for(j=0; j<fromEdges->length; j++) {
        edge = fromEdges->list[j];
        assert(edgeFrom->from < edge->to);
        edge->from = edgeFrom->from;
        edge->edgeScore += edgeFrom->edgeScore;
        assert(edge->treeNode->left == NULL);
        edge->treeNode->left = edgeFrom->treeNode;
        edgeFrom->treeNode->refCount++;
    }
    for(i=0; i<toEdges->length; i++) {
        destructEdge(toEdges->list[i]); //clean up incoming edges
    }
}

struct SequenceGraph *compactSilentVertices(struct SequenceGraph *sequenceGraph) {
    //simplifies silent-edge/vertice subgraphs by removing greedily vertices whose
    //combined incoming and outgoing edges would be better served by direct connections
    int32_t vertex;
    struct List *toEdges;
    struct List *fromEdges;
    struct List *newEdges;
    int32_t vertexShift = 0;
    int32_t *verticeShifts;
    struct Edge *edge;
    int32_t i;
    struct SequenceGraph *finalSequenceGraph;

    verticeShifts = st_calloc(sequenceGraph->vertexNo, sizeof(int32_t));
    newEdges = sequenceGraph->edges;
    newEdges->length = 0;
    for(vertex=0; vertex <sequenceGraph->vertexNo; vertex++) {
        toEdges = sequenceGraph->edgesArrangedByToVertex[vertex];
        fromEdges = sequenceGraph->edgesArrangedByFromVertex[vertex];
        if(fromEdges->length > 0 && toEdges->length > 0 && ((struct Edge *)toEdges->list[0])->silent &&
        //((struct Edge *)fromEdges->list[0])->silent &&
        fromEdges->length + toEdges->length + 1 >= fromEdges->length * toEdges->length) { //the incoming edges are silent, so we can skip them out if needed
            vertexShift++;
            verticeShifts[vertex] = INT32_MAX;
            compactVertex(toEdges, fromEdges, sequenceGraph);
        }
        else {
            listAppendArray(newEdges, toEdges->list, toEdges->length);
            verticeShifts[vertex] = vertexShift;
        }
    }
    for(i=0; i<newEdges->length; i++) {
        edge = newEdges->list[i];
        edge->iD = i;
        edge->from -= verticeShifts[edge->from];
        edge->to -= verticeShifts[edge->to];
    }
    finalSequenceGraph = constructSequenceGraph(newEdges, sequenceGraph->vertexNo-vertexShift);

    //memory clean up
    destructSequenceGraph(sequenceGraph, FALSE);
    free(verticeShifts);

    return finalSequenceGraph;
}

void convertTransitionIDToStates(int32_t stateNo, int32_t z, int32_t *fromState, int32_t *toState) {
    *fromState = z/stateNo;
    *toState = z%stateNo;
}

struct SequenceGraph* convertTraceBackEdgesToSequenceGraph(struct AlignmentDataStructures *aDS, void **newEdges, int32_t newEdgeNumber) {
    int32_t i;
    int32_t transitionID;
    struct Edge *edge;
    struct TraceBackEdge *traceBackEdge;
    struct TreeNode *treeNode;
    float wV[SMALL_CHUNK_SIZE];
    //
    struct SequenceGraph *finalSequenceGraph;
    int64_t vertex;
    int32_t compactVertex;
    struct hashtable *vertexMap;
    struct Chunks *intChunks;
    struct Chunks *longChunks;

    vertex=INT64_MIN;
    compactVertex=0;

    intChunks = constructChunks(newEdgeNumber + 1, sizeof(int32_t));
    longChunks = constructChunks(newEdgeNumber + 1, sizeof(int64_t));
    vertexMap = create_hashtable(newEdgeNumber*2 + 1, hashtable_longHashKey, hashtable_longEqualKey, NULL, NULL);
    hashtable_insert(vertexMap, constructChunkLong(0, longChunks), constructChunkInt(0, intChunks));

    for(i=0; i<newEdgeNumber; i++) {
        traceBackEdge = newEdges[i];
        transitionID = stateNo()*rVZ_State(aDS, traceBackEdge->from) + rVZ_State(aDS, traceBackEdge->to);
        if(traceBackEdge->to != vertex) {
            assert(!hashtable_search(vertexMap, &traceBackEdge->to));
            hashtable_insert(vertexMap, constructChunkLong((vertex = traceBackEdge->to), longChunks), constructChunkInt(++compactVertex, intChunks));
            traceBackEdge->to = compactVertex;
        }
        else {
            traceBackEdge->to = compactVertex;
        }
        traceBackEdge->from = *((int32_t *)hashtable_search(vertexMap, &traceBackEdge->from));
        if(traceBackEdge->silent) {
        	//treeNode = ((struct TreeNode *(*)(struct AlignmentDataStructures *, struct TraceBackEdge *, int32_t))traceBackEdge->getTreeNode)(aDS, traceBackEdge, transitionID);
            treeNode = traceBackEdge->getTreeNode(aDS, traceBackEdge, transitionID);
            edge = copyConstructEdge(traceBackEdge->from, traceBackEdge->to, traceBackEdge->edgeScore,
            LOG_ZERO, NULL, TRUE, treeNode, i);
        }
        else {
            treeNode = ((struct TreeNode *(*)(struct AlignmentDataStructures *, struct TraceBackEdge *, int32_t, float *))traceBackEdge->getTreeNode)(aDS, traceBackEdge, transitionID, wV);
            edge = copyConstructEdge(traceBackEdge->from, traceBackEdge->to, traceBackEdge->edgeScore,
            LOG_ZERO, //these values not yet set
            wV, FALSE, treeNode, i);
        }
        destructTraceBackEdge(traceBackEdge);
        newEdges[i] = edge;
    }
    i = hashtable_count(vertexMap);
    finalSequenceGraph = constructSequenceGraph(copyConstructList((void **)newEdges, newEdgeNumber, (void (*)(void *))destructEdge), i);
    st_logDebug("Total number of edges and vertices (without compaction) : " INT_STRING " , " INT_STRING " \n", finalSequenceGraph->edges->length, finalSequenceGraph->vertexNo);
    finalSequenceGraph = compactSilentVertices(finalSequenceGraph);
    st_logDebug("Total number of edges and vertices (after compaction) : " INT_STRING " , " INT_STRING " \n", finalSequenceGraph->edges->length, finalSequenceGraph->vertexNo);
    //memory clean up
    hashtable_destroy(vertexMap, FALSE, FALSE);
    destructChunks(intChunks);
    destructChunks(longChunks);
    return finalSequenceGraph;
}

struct SequenceGraph *traceBackMatrix(struct AlignmentDataStructures *aDS) {
    int32_t i;
    int32_t j;
    int64_t iL;
    int32_t state;
    int32_t pathWeight;
    int32_t pathWeight2;
    int32_t toX;
    int32_t toY;
    int32_t toZ;
    int32_t fromX;
    int32_t fromY;
    int32_t *temp;
    int32_t hashDummy;

    static void **newEdges; //do not clean up returned array
    static int32_t newEdgesSize;
    int32_t newEdgesIndex = 0;

    struct hashtable *pathWeightsHash;
    struct heap *vertexHeap;
    int64_t finalVertex;
    int64_t terminationVertex;
    float *finalVertices;
    float *endProbs;
    struct hashtable *startStateNos;
    int64_t previous;
    int32_t *choices;

    struct List *edgesX;
    struct List *edgesY;
    struct TraceBackEdge *traceBackEdge;
    struct Chunks *intChunks;
    struct Chunks *longChunks;

    newEdges = arrayResize(newEdges, &newEdgesSize, (aDS->sequenceGraphY->edges->length + aDS->sequenceGraphX->edges->length + stateNo() + 2)*10, sizeof(void *)); //estimated max

    finalVertex = vZ(aDS, aDS->sequenceGraphX->vertexNo-1, aDS->sequenceGraphY->vertexNo-1, 0); //starting vertices
    terminationVertex = finalVertex + 100*stateNo(); //vZ(aDS, aDS->sequenceGraphX->vertexNo-1, aDS->sequenceGraphY->vertexNo-1, 0); //starting vertices
    //special vertex used for creating extra vertices for double-deletes
    finalVertices = getCell(aDS, rVZ_Zless(aDS, finalVertex));
    endProbs = st_malloc(sizeof(float)*stateNo());
    for (state = 0; state < stateNo(); ++state) {
        endProbs[state] = finalVertices[state] + aDS->endStates[state];
    }
    //for state, pathWeight in randomChoices(endProbs, aDS.numberOfSamples) {
    vertexHeap = heap_create(aDS->numberOfSamples+MEDIUM_CHUNK_SIZE);
    pathWeightsHash = create_hashtable(MEDIUM_CHUNK_SIZE, hashtable_longHashKey, hashtable_longEqualKey, NULL, NULL);
    choices = randomChoices(endProbs, stateNo(), aDS->numberOfSamples);

    intChunks = constructChunks(MEDIUM_CHUNK_SIZE, sizeof(int32_t));
    longChunks = constructChunks(MEDIUM_CHUNK_SIZE, sizeof(int64_t));

    for (state = 0; state < stateNo(); state++) {
        pathWeight = choices[state];
        if (pathWeight > 0) {
            iL = finalVertex + state;
            heap_insert(vertexHeap, iL);  //heap.pushOnHeap(iL);
            hashtable_insert(pathWeightsHash, constructChunkLong(iL, longChunks), constructChunkInt(pathWeight, intChunks)); //pathWeightsHash[i] = pathWeight;
            newEdges[newEdgesIndex++] = constructTraceBackEdge(iL, terminationVertex, LOG_ONE, NULL, NULL, TRUE, getTreeNode_silentXY); //end vertex is 0 state
        }
    }

    startStateNos = create_hashtable(stateNo(), hashtable_intHashKey, hashtable_intEqualKey, NULL, NULL);
    //for state in xrange(0, aDS.stateNo) {
    if(aDS->treeStates == NULL || TRUE) {
        for (state=0; state < stateNo(); state++) {
            if (aDS->startStates[state] > LOG_ZERO) {
                hashtable_insert(startStateNos, constructChunkInt(state, intChunks), &hashDummy);
            }
        }
    }
    else {
        //hashtable_insert(startStateNos, constructChunkInt(0), &hashDummy);
        hashtable_insert(startStateNos, constructChunkInt(aDS->treeStates[aDS->traversalID->mid], intChunks), &hashDummy);
    }
    //main loop
    previous = INT64_MAX;
    while (!heap_empty(vertexHeap)) {
        aDS->toCombined = heap_extract(vertexHeap); // heap.popTopOfHeap()
        //while not heap.empty() and heap.peekTopOfHeap() == aDS.toCombined {
        while (!heap_empty(vertexHeap) && heap_peek(vertexHeap) == aDS->toCombined) {
            heap_extract(vertexHeap); //heap.popTopOfHeap()
        }
        //uglyf(" current, " LONG_INT_STRING "  previous " LONG_INT_STRING " \n", aDS->toCombined, previous);
        assert(aDS->toCombined <= previous);
        previous = aDS->toCombined;
        aDS->to = rVZ_StatelessVertex(aDS, aDS->toCombined);
        rVZ(aDS, aDS->toCombined, &toX, &toY, &toZ, &aDS->state);
        assert(getCell(aDS, rVZ_Zless(aDS, aDS->to))[aDS->state] != LOG_ZERO);
        temp = hashtable_remove(pathWeightsHash, &aDS->toCombined, FALSE);
        pathWeight = *temp;
        assert(pathWeight >= 1);
        if (rVZ_Zless(aDS, aDS->to) == 0 && hashtable_search(startStateNos, &aDS->state) != NULL) {
            continue;
        }
        //for safety
        aDS->silent = FALSE;
        aDS->getTreeNode = NULL;
        aDS->edgeX = NULL;
        aDS->edgeY = NULL;
        aDS->from = INT32_MAX;
        //end
        edgesX = aDS->sequenceGraphX->edgesArrangedByToVertex[toX];
        edgesY = aDS->sequenceGraphY->edgesArrangedByToVertex[toY];
        if (aDS->sequenceGraphXSilentVertices_To[toX]) {
            aDS->silent = TRUE;
            aDS->getTreeNode = getTreeNode_silentX;
            //for edgeX in edgesX {
            for (i = 0; i < edgesX->length; ++i) {
                aDS->edgeX = edgesX->list[i];
                fromX = aDS->edgeX->from;
                aDS->from = vZ(aDS, fromX, toY, 0); //toZ);
                aDS->fromCell = getCell(aDS, rVZ_Zless(aDS, aDS->from));
                if(aDS->fromCell != NULL) {
                    assignSampling(aDS, aDS->state, aDS->state,
                                   aDS->edgeX->edgeScore, aDS->edgeX->edgeScore);
                }
            }
        }
        else if (aDS->sequenceGraphYSilentVertices_To[toY]) {
            aDS->silent = TRUE;
            aDS->getTreeNode = getTreeNode_silentY;
            //for edgeY in edgesY {
            for (i = 0; i < edgesY->length; ++i) {
                aDS->edgeY = edgesY->list[i];
                fromY = aDS->edgeY->from;
                aDS->from = vZ(aDS, toX, fromY, 0); //toZ);
                aDS->fromCell = getCell(aDS, rVZ_Zless(aDS, aDS->from));
                if(aDS->fromCell != NULL) {
                    assignSampling(aDS, aDS->state, aDS->state,
                                   aDS->edgeY->edgeScore, aDS->edgeY->edgeScore);
                }
            }
        }
        else if (isSilent(aDS->state)) {
            aDS->silent = TRUE;
            aDS->getTreeNode = getTreeNode_silentXY;
            aDS->from = aDS->to;
            aDS->fromCell = getCell(aDS, rVZ_Zless(aDS, aDS->from));
            if(aDS->fromCell != NULL) {
                silentFn_TraceBack(aDS, aDS->model, assignSampling_CheckState);
            }
        }
        else if (isXInsert(aDS->state)) {
            aDS->silent = TRUE;
            aDS->getTreeNode = getTreeNode_insertX;
            //for edgeX in edgesX {
            for (i = 0; i < edgesX->length; ++i) {
                aDS->edgeX = edgesX->list[i];
                fromX = aDS->edgeX->from;
                aDS->from = vZ(aDS, fromX, toY, 0); //toZ);
                aDS->fromCell = getCell(aDS, rVZ_Zless(aDS, aDS->from));
                if(aDS->fromCell != NULL) {
                    insertXFn_TraceBack(aDS, aDS->model, aDS->edgeX, assignSampling_CheckState);
                }
            }
        }
        else if (isYDelete(aDS->state)) {
            aDS->silent = FALSE;
            //aDS->getTreeNode = getTreeNode_deleteY;
            aDS->getTreeNode = (struct TreeNode *(*)(struct AlignmentDataStructures *, struct TraceBackEdge*, int32_t))getTreeNode_deleteY;
            //for edgeX in edgesX {
            for (i = 0; i < edgesX->length; ++i) {
                aDS->edgeX = edgesX->list[i];
                fromX = aDS->edgeX->from;
                aDS->from = vZ(aDS, fromX, toY, 0); //toZ);
                aDS->fromCell = getCell(aDS, rVZ_Zless(aDS, aDS->from));
                if(aDS->fromCell != NULL) {
                    deleteYFn_TraceBack(aDS, aDS->model, aDS->edgeX, assignSampling_CheckState);
                }
            }
        }
        else if (isYInsert(aDS->state)) {
            aDS->silent = TRUE;
            aDS->getTreeNode = getTreeNode_insertY;
            //for edgeY in edgesY {
            for (i = 0; i < edgesY->length; ++i) {
                aDS->edgeY = edgesY->list[i];
                fromY = aDS->edgeY->from;
                aDS->from = vZ(aDS, toX, fromY, 0); //toZ);
                aDS->fromCell = getCell(aDS, rVZ_Zless(aDS, aDS->from));
                if(aDS->fromCell != NULL) {
                    insertYFn_TraceBack(aDS, aDS->model, aDS->edgeY, assignSampling_CheckState);
                }
            }
        }
        else if (isXDelete(aDS->state)) {
            aDS->silent = FALSE;
            //aDS->getTreeNode = getTreeNode_deleteX;
            aDS->getTreeNode = (struct TreeNode *(*)(struct AlignmentDataStructures *, struct TraceBackEdge*, int32_t))getTreeNode_deleteX;
            //for edgeY in edgesY {
            for (i = 0; i < edgesY->length; ++i) {
                aDS->edgeY = edgesY->list[i];
                fromY = aDS->edgeY->from;
                aDS->from = vZ(aDS, toX, fromY, 0); //toZ);
                aDS->fromCell = getCell(aDS, rVZ_Zless(aDS, aDS->from));
                if(aDS->fromCell != NULL) {
                    deleteXFn_TraceBack(aDS, aDS->model, aDS->edgeY, assignSampling_CheckState);
                }
            }
        }
        else if (isMatch(aDS->state)) {  //is match
            aDS->silent = FALSE;
            //aDS->getTreeNode = getTreeNode_matchXY;
            aDS->getTreeNode = (struct TreeNode *(*)(struct AlignmentDataStructures *, struct TraceBackEdge*, int32_t))getTreeNode_matchXY;
            //for edgeX in edgesX {
            for (i = 0; i < edgesX->length; ++i) {
                aDS->edgeX = edgesX->list[i];
                fromX = aDS->edgeX->from;
                //for edgeY in edgesY {
                for (j = 0; j < edgesY->length; ++j) {
                    aDS->edgeY = edgesY->list[j];
                    fromY = aDS->edgeY->from;
                    aDS->from = vZ(aDS, fromX, fromY, 0); //toZ);
                    aDS->fromCell = getCell(aDS, rVZ_Zless(aDS, aDS->from));
                    if(aDS->fromCell != NULL) {
                        matchFn_TraceBack(aDS, aDS->model, aDS->edgeX, aDS->edgeY, assignSampling_CheckState);
                    }
                }
            } //is match
        }
        else {
            //logDebug("delete-delete state " INT_STRING " \n", aDS->state);
            //exit(0);
            assert(isXYDelete(aDS->state));
            aDS->silent = FALSE;
            //aDS->getTreeNode = getTreeNode_deleteXY;
            aDS->getTreeNode = (struct TreeNode *(*)(struct AlignmentDataStructures *, struct TraceBackEdge*, int32_t))getTreeNode_deleteXY;
            aDS->from = vZ(aDS, toX, toY, toZ+1);
            if(toZ+1 >= MAX_Z) {
            	fprintf(stderr, "Reached max z-threshold\n");
            	exit(1);
            }
            aDS->fromCell = getCell(aDS, rVZ_Zless(aDS, aDS->from));
            deleteFn_TraceBack(aDS, aDS->model, assignSampling_CheckState);
        }
        choices = randomChoices(aDS->potentialEdgeCosts, aDS->potentialEdges_Index, pathWeight);
        for (i = 0; i < aDS->potentialEdges_Index; i++) {
            pathWeight2 = choices[i];
            traceBackEdge = aDS->potentialEdges[i];
            if (pathWeight2 > 0) {
                newEdges = arrayPrepareAppend(newEdges, &newEdgesSize, newEdgesIndex+1, sizeof(void *));
                //remove the bias introduced by the correction score
                newEdges[newEdgesIndex++] = traceBackEdge;
                heap_insert(vertexHeap, traceBackEdge->from); //heap.pushOnHeap(j);

                if ((temp = hashtable_search(pathWeightsHash, &traceBackEdge->from)) != NULL) {
                    (*temp) += pathWeight2;
                    //uglyf(" from to3 " LONG_INT_STRING " " LONG_INT_STRING " \n", aDS->from, aDS->to);
                    //pathWeightsHash[aDS.from] += pathWeight2
                }
                else {
                	 //uglyf(" from to2 " LONG_INT_STRING " " LONG_INT_STRING " \n", aDS->from, aDS->to);
                    hashtable_insert(pathWeightsHash, constructChunkLong(traceBackEdge->from, longChunks), constructChunkInt(pathWeight2, intChunks));
                    //pathWeightsHash[aDS.from] = pathWeight2
                }
            }
            else { //clean up edge
                aDS->potentialEdges[i] = NULL;
                destructTraceBackEdge(traceBackEdge);
            }
        }
        aDS->potentialEdges_Index = 0;
    }
    //correct start and end conditions
    //for edge in newEdges {
    for (i=0; i<newEdgesIndex; i++) {
        traceBackEdge = newEdges[i];
        j = rVZ_State(aDS, traceBackEdge->from);
        if (rVZ_Zless(aDS, traceBackEdge->from) == 0 && hashtable_search(startStateNos, &j) != NULL) {
            traceBackEdge->from = 0;
        }
    }
    //convert edges to real edges
    for(i=0; i<newEdgesIndex/2; i++) { //alignment.reverse()
        traceBackEdge = newEdges[i];
        newEdges[i] = newEdges[newEdgesIndex-1-i];
        newEdges[newEdgesIndex-1-i] = traceBackEdge;
    }
    //final memory clean up
    //hashtables
    //hashtable_destroy(newVerticesHash, FALSE, FALSE);
    hashtable_destroy(pathWeightsHash, FALSE, FALSE);
    hashtable_destroy(startStateNos, FALSE, FALSE);
    //heaps
    heap_destroy(vertexHeap);
    //pieces of array
    free(endProbs);
    //chunks arrays
    destructChunks(intChunks);
    destructChunks(longChunks);
    //now convert newEdges to proper edges -- not implemented yet in C
    return convertTraceBackEdgesToSequenceGraph(aDS, newEdges, newEdgesIndex);
}

int32_t wrapInGraphHolders_Length;
int wrapInGraphHolders_intComp(int32_t **iAA, int32_t **iAA2) {
    int32_t i;
    int32_t *iA = *iAA;
    int32_t *iA2 = *iAA2;

    for(i=0; i<wrapInGraphHolders_Length; i++) {
        if(iA[i] < iA2[i]) {
            return -1;
        }
        if(iA[i] > iA2[i]) {
            return 1;
        }
    }
    assert(FALSE);
}

void wrapInGraphHolders(struct AlignmentDataStructures *aDS, struct List *list, void *key, int32_t(*comp)(struct AlignmentDataStructures *, struct GraphMemberHolder *, struct GraphMemberHolder *)) {
    int32_t i;

    wrapInGraphHolders_Length = aDS->leafSeqNoX;
    qsort(list->list, list->length, sizeof(void *), (int (*)(const void *, const void *))wrapInGraphHolders_intComp);
    for(i=0; i<list->length; i++) {
        list->list[i] = constructGraphMember(key, list->list[i], NULL);
    }
    if(DEBUG) {
        for(i=1; i<list->length; i++) {
            assert(comp(aDS, list->list[i-1], list->list[i]) < 0);
        }
    }
}

void prepareGraphForAlignment(struct SequenceGraph *sequenceGraphX, struct CombinedTransitionModel *model, int32_t flipped) {
    int32_t i;
    struct Edge *edge;
    struct SubModel *subModelX;
    struct SubModel *subModelY;

    subModelX = model->subModelX;
    subModelY = model->subModelY;
    if(flipped) {
    	subModelX = model->subModelY;
    	subModelY = model->subModelX;
    }
    for(i=0; i<sequenceGraphX->edges->length; i++) {
        edge = sequenceGraphX->edges->list[i];
        if(!edge->silent) {
            transformWVByDistance(edge->wV, subModelX->backward, edge->wV, subModelX->alphabetSize);
            edge->subScore = LOG(combineWV(edge->wV, subModelX->stationaryDistribution, subModelX->alphabetSize));
            //edge->deleteBranchCost = combineWV(edge->wV, subModelY->deletionDistribution);
        }
    }
}

struct SequenceGraph *computeEdgeGraph(struct SequenceGraph *sequenceGraphX, struct SequenceGraph *sequenceGraphY,
                                      struct CombinedTransitionModel *model,
                                      struct Constraints ***allConstraints,
                                      //struct SubModel **subModels,
                                      struct TraversalID *traversalID,
                                      struct TraversalID *traversalIDX,
                                      struct TraversalID *traversalIDY,
                                      int32_t leftMostSeqNo, int32_t leafSeqNoX, int32_t leafSeqNoY,
                                      int32_t numberOfSamples, int32_t *treeStates) {
    int32_t i;
    struct Chunks *mergedChunksX;
    struct Chunks *mergedChunksY;
    struct Chunks *intChunks;
    int32_t **vertexYSequenceCoordinates;
    struct SequenceGraph *newSequenceGraph;
    struct AlignmentDataStructures *aDS;
    //computes a specified number of sample alignments using the forward algorithm
    st_logDebug("Number of edgesX: " INT_STRING " Number of edgesY: " INT_STRING " Number of verticesX: " INT_STRING " Number of verticesY: " INT_STRING "\n",
           sequenceGraphX->edges->length, sequenceGraphY->edges->length, sequenceGraphX->vertexNo, sequenceGraphY->vertexNo);

    aDS = st_malloc(sizeof(struct AlignmentDataStructures));

    aDS->sequenceGraphX = sequenceGraphX;
    aDS->sequenceGraphY = sequenceGraphY;
    prepareGraphForAlignment(sequenceGraphX, model, FALSE);
    prepareGraphForAlignment(sequenceGraphY, model, TRUE);

    //struct Edge *edge;
    //for(i=0; i<aDS->sequenceGraphX->edges->length; i++) {
    //    edge = aDS->sequenceGraphX->edges->list[i];
    //    fprintf(stderr, "edgeX " INT_STRING " " INT_STRING " " INT_STRING " " FLOAT_STRING " " INT_STRING " \n", edge->from, edge->to, edge->silent, edge->edgeScore, edge->treeNode);
    //    printTreeNode(edge->treeNode);
    //}
    //for(i=0; i<aDS->sequenceGraphY->edges->length; i++) {
    //    edge = aDS->sequenceGraphY->edges->list[i];
    //    fprintf(stderr, "edgeY " INT_STRING " " INT_STRING " " INT_STRING " " FLOAT_STRING " " INT_STRING " \n", edge->from, edge->to, edge->silent, edge->edgeScore, edge->treeNode);
    //    printTreeNode(edge->treeNode);
    //}
    //void *allConstraints;
    aDS->model = model;
    //aDS->subModel = subModels[traversalID->mid];

    aDS->startStates = startStates(model);
    aDS->endStates = endStates(model);

    //aDS->traversalIDX = traversalIDX;
    //aDS->traversalIDY = traversalIDY;
    aDS->traversalID = traversalID;
    aDS->leftMostSeqNoX = leftMostSeqNo;
    aDS->leftMostSeqNoY = leftMostSeqNo + leafSeqNoX;
    aDS->leafSeqNoX = leafSeqNoX;
    aDS->leafSeqNoY = leafSeqNoY;

    aDS->numberOfSamples = numberOfSamples;
    aDS->matrixChunks = constructChunks(LARGE_CHUNK_SIZE, sizeof(int64_t) + sizeof(float)*stateNo());//constructEmptyList(0, free);
    mergedChunksX = constructChunks(MEDIUM_CHUNK_SIZE, sizeof(int32_t)*aDS->leafSeqNoX);
    mergedChunksY = constructChunks(MEDIUM_CHUNK_SIZE, sizeof(int32_t)*aDS->leafSeqNoY);
    aDS->vertexXSequenceCoordinates = calculateVertexSequenceCoordinates(sequenceGraphX, aDS->leftMostSeqNoX, aDS->leafSeqNoX, mergedChunksX);
    vertexYSequenceCoordinates = calculateVertexSequenceCoordinates(sequenceGraphY, aDS->leftMostSeqNoY, aDS->leafSeqNoY, mergedChunksY);

    aDS->sequenceGraphXSilentVertices_To = calculateSilentVertices(sequenceGraphX);//sequenceGraphXSilentVertices_To;
    aDS->sequenceGraphYSilentVertices_To = calculateSilentVertices(sequenceGraphY);//sequenceGraphYSilentVertices_To;
    //these datastructures are used to compute the constraINT_32 envelope of the alignments during matrix computation

    //used to allocate memory
    calculateMergedConstraints_RightToLeft(vertexYSequenceCoordinates,
                                           aDS->leftMostSeqNoY, aDS->leftMostSeqNoX, aDS->leafSeqNoY, aDS->leafSeqNoX,
                                           allConstraints, sequenceGraphY,
                                           &aDS->mergedStartConstraints_Vertices, &aDS->mergedStartConstraints_Edges, mergedChunksX);
    intChunks = constructChunks(sequenceGraphY->vertexNo, sizeof(int32_t));
    for (i=0; i < sequenceGraphY->vertexNo; i++) {
        ((struct List *)aDS->mergedStartConstraints_Vertices[i])->destructElement = (void (*)(void *))destructGraphMember;
        wrapInGraphHolders(aDS, aDS->mergedStartConstraints_Vertices[i], constructChunkInt(i, intChunks), graphMember_VertexComparator); //leaf numbers must be initialised at this point
    }
    for (i=0; i < sequenceGraphY->edges->length; i++) {
        ((struct List *)aDS->mergedStartConstraints_Edges[i])->destructElement = (void (*)(void *))destructGraphMember;
        wrapInGraphHolders(aDS, aDS->mergedStartConstraints_Edges[i], sequenceGraphY->edges->list[i], graphMember_EdgeComparator);
    }
    calculateMergedConstraints_LeftToRight(vertexYSequenceCoordinates,
                                           aDS->leftMostSeqNoY, aDS->leftMostSeqNoX, aDS->leafSeqNoY, aDS->leafSeqNoX,
                                           allConstraints, sequenceGraphY,
                                           &aDS->mergedEndConstraints_Vertices, &aDS->mergedEndConstraints_Edges, mergedChunksX);
    calculateMergedSequenceCoordinates_RightToLeft(aDS->vertexXSequenceCoordinates, sequenceGraphX, aDS->leftMostSeqNoX, aDS->leafSeqNoX,
                                                   &aDS->mergedEndSequenceCoordinates_Vertices, &aDS->mergedEndSequenceCoordinates_Edges,
                                                   mergedChunksX);
    //added to avoid repeated access
    aDS->vertexXNo = sequenceGraphX->vertexNo;
    aDS->vertexYNo = sequenceGraphY->vertexNo;

    //used in scanning
    aDS->newVertices_Size = MEDIUM_CHUNK_SIZE;
    aDS->newVertices = st_malloc(sizeof(int32_t) * aDS->newVertices_Size);
    aDS->noOfNewVertices = 0;
    aDS->changes = st_malloc(sizeof(int32_t) * leafSeqNoX * 2);
    aDS->noOfChanges = 0;
    aDS->sequenceCoordinatesCollection = NULL;
    aDS->mergedEndConstraints_GraphMembers = NULL;
    aDS->getIDFromGraphMember = NULL;

    //transition starts and ends, used in compute matrix and traceback
    aDS->to = INT64_MAX;
    aDS->from = INT64_MAX;
    aDS->toCombined = INT64_MAX;

    //following used by the traceback methods
    //aDS->potentialEdges;
    //aDS->potentialEdgeCosts;
    aDS->potentialEdges_Size = MEDIUM_CHUNK_SIZE;
    aDS->potentialEdges = st_malloc(sizeof(void *)*aDS->potentialEdges_Size);
    aDS->potentialEdgeCosts = st_malloc(sizeof(float)*aDS->potentialEdges_Size);
    aDS->potentialEdges_Index = 0;
    //branches for composing tree, used in traceback
    aDS->deleteNodeX = copyConstructTreeNode(TREE_NODE_DELETE, INT32_MAX, traversalIDX, NULL, NULL, NULL);
    //aDS->subModelX = subModels[traversalIDX->mid];
    aDS->deleteNodeY = copyConstructTreeNode(TREE_NODE_DELETE, INT32_MAX, traversalIDY, NULL, NULL, NULL);
    //aDS->subModelY = subModels[traversalIDY->mid];
    aDS->state = INT32_MAX;
    aDS->edgeX = NULL;
    aDS->edgeY = NULL;
    aDS->treeStates = treeStates;

    //the important calls
    computeMatrix(aDS);
    newSequenceGraph = traceBackMatrix(aDS);
    //memory clean up
    destructChunks(aDS->matrixChunks);
    //destructList(aDS->matrixChunks);

    for (i=0; i < sequenceGraphY->vertexNo; i++) {
        destructList(aDS->mergedEndConstraints_Vertices[i]);
        destructList(aDS->mergedStartConstraints_Vertices[i]);
    }
    free(aDS->mergedStartConstraints_Vertices);
    free(aDS->mergedEndConstraints_Vertices);
    destructChunks(intChunks);

    for (i=0; i < sequenceGraphY->edges->length; i++) {
        destructList(aDS->mergedEndConstraints_Edges[i]);
        destructList(aDS->mergedStartConstraints_Edges[i]);
    }
    free(aDS->mergedStartConstraints_Edges);
    free(aDS->mergedEndConstraints_Edges);

    cleanUpSequenceCoordinates(aDS->sequenceGraphX, aDS->mergedEndSequenceCoordinates_Vertices,
                               aDS->mergedEndSequenceCoordinates_Edges);

    free(aDS->sequenceGraphXSilentVertices_To);
    free(aDS->sequenceGraphYSilentVertices_To);

    free(vertexYSequenceCoordinates);
    free(aDS->vertexXSequenceCoordinates);

    destructChunks(mergedChunksX);
    destructChunks(mergedChunksY);

    destructTreeNode(aDS->deleteNodeX);
    destructTreeNode(aDS->deleteNodeY);

    free(aDS->startStates);
    free(aDS->endStates);
    free(aDS->newVertices);
    free(aDS->changes);
    free(aDS->potentialEdges);
    free(aDS->potentialEdgeCosts);
    free(aDS);

    //end memory clean up
    return newSequenceGraph;
}

struct List *viterbi(struct SequenceGraph *sequenceGraph, float *finalScore) {
    //get the viterbi alignment from the sequence graph
    //assumes alignment starts from 0 vertex and must end in highest
    //numbered vertex
    //all in order
    float *viterbiMatrix;
    void **pointers;
    int32_t i;
    int32_t j;
    struct Edge *edge;
    float score;
    void **alignment;
    int32_t vertex;
    struct List *list;

    viterbiMatrix = st_malloc(sizeof(float)*sequenceGraph->vertexNo);
    for(i=0; i<sequenceGraph->vertexNo; i++) {
        viterbiMatrix[i] = LOG_ZERO;
    }
    pointers = st_calloc(sequenceGraph->vertexNo, sizeof(void *));
    alignment = st_malloc(sizeof(void *)*sequenceGraph->vertexNo);

    //zero is terminator
    pointers[0] = 0;
    viterbiMatrix[0] = LOG_ONE;
    //for edge in sequenceGraph.edges:
    for(i=0; i<sequenceGraph->edges->length; i++) {
        edge = sequenceGraph->edges->list[i];
        score = viterbiMatrix[edge->from] + edge->edgeScore;
        if(score > viterbiMatrix[edge->to]) {
            viterbiMatrix[edge->to] = score;
            pointers[edge->to] = edge;
        }
    }
    vertex = sequenceGraph->vertexNo-1;
    i=0;
    while(pointers[vertex] != 0) {
        edge = pointers[vertex];
        alignment[i++] = edge;
        vertex = edge->from;
    }
    *finalScore = viterbiMatrix[sequenceGraph->vertexNo-1];
    for(j=0; j<i/2; j++) { //alignment.reverse()
        edge = alignment[j];
        alignment[j] = alignment[i-1-j];
        alignment[i-1-j] = edge;
    }
    list = copyConstructList(alignment, i, NULL);
    //memory clean up
    free(viterbiMatrix);
    free(pointers);
    free(alignment);
    return list;
}

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//graph iterative/em thing
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

float *calcForwardMatrix(struct SequenceGraph *sequenceGraph) {
    /*calculates the forward values for the vertices and edges
    *in the graph
    */
    //all in order
    int32_t i;
    float *forwardMatrix;
    struct Edge *edge;

    forwardMatrix = st_malloc(sequenceGraph->vertexNo*sizeof(float));
    for(i=0; i<sequenceGraph->vertexNo; i++) {
        forwardMatrix[i] = LOG_ZERO;
    }
    //zero vertex is always terminator
    forwardMatrix[0] = LOG_ONE;

    for(i=0; i<sequenceGraph->edges->length; i++) {
        edge = sequenceGraph->edges->list[i];
        //for edge in sequenceGraph.edges:
        LOG_PLUS_EQUALS(&forwardMatrix[edge->to], forwardMatrix[edge->from] + edge->edgeScore);
    }
    return forwardMatrix;
}

float *calcBackwardMatrix(struct SequenceGraph *sequenceGraph) {
    /*calculates the forward values for the vertices and edges
    *in the graph
    */
    //all in order
    int32_t i;
    float *backwardMatrix;
    struct Edge *edge;

    backwardMatrix = st_malloc(sequenceGraph->vertexNo*sizeof(float));
    for(i=0; i<sequenceGraph->vertexNo; i++) {
        backwardMatrix[i] = LOG_ZERO;
    }
    //zero vertex is always terminator
    backwardMatrix[sequenceGraph->vertexNo-1] = LOG_ONE;

    for(i=sequenceGraph->edges->length-1; i>=0; i--) {
        edge = sequenceGraph->edges->list[i];
        //for edge in sequenceGraph.edges:
        LOG_PLUS_EQUALS(&backwardMatrix[edge->from], backwardMatrix[edge->to] + edge->edgeScore);
    }
    return backwardMatrix;
}

float *posteriorProbabilities(float *forwardMatrix, float *backwardMatrix, struct SequenceGraph *sequenceGraph) {
    /*
    *calculates the posterior probabilities of edges in the
    *sequence graph
    */
    int32_t i;
    float total;
    float *edgeProbs;
    struct Edge *edge;

    total = forwardMatrix[sequenceGraph->vertexNo-1];
    assert(total <= backwardMatrix[0] + 0.00001);
    assert(total >= backwardMatrix[0] - 0.00001);
    edgeProbs = st_malloc(sequenceGraph->edges->length*sizeof(float));
    for(i=0; i<sequenceGraph->edges->length; i++) {
    //for(edge in sequenceGraph.edges) {
        edge = sequenceGraph->edges->list[i];
        edgeProbs[edge->iD] = forwardMatrix[edge->from] + backwardMatrix[edge->to] + edge->edgeScore - total;
    }
    return edgeProbs;
}

int32_t tranCoord(int32_t from, int32_t to, int32_t node, int32_t stateNo, int32_t nodeNo) {
    return nodeNo * stateNo * stateNo * node + stateNo * from + to;
}


/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//felsensteins
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

float *upPass(struct TreeNode *treeNode, float *working, struct SubModel **subModels, float *totalP,  int32_t alphabetSize);

void branchUp(struct TreeNode *treeNode, float *working, float *currentResult, struct SubModel **subModels, float *totalP,  int32_t alphabetSize) {
    float *i;
    int32_t j;

    if (treeNode->type != TREE_NODE_DELETE) {
        //return subModels[treeNode->traversalID->mid].transformForwardInTime(upPass(treeNode)); --python error!
        i = upPass(treeNode, working, subModels, totalP, alphabetSize);
        transformWVByDistance(i, subModels[treeNode->traversalID->mid]->backward, currentResult, alphabetSize);
    }
    else {
        //i = subModels[treeNode->traversalID->mid]->deletionDistribution;
        //memcpy(currentResult, i, sizeof(FLOAT_32)*ALPHABET_SIZE);
        for(j=0; j<alphabetSize; j++) {
        	currentResult[j] = 1.0f;
        }
        totalP[treeNode->traversalID->mid] = 0.0f;
    }
    //return subModels[treeNode->traversalID->mid].deletionDistribution
}

//calculates the un-normalised probabilties of each non-gap residue position
float *upPass(struct TreeNode *treeNode, float *working, struct SubModel **subModels, float *totalP, int32_t alphabetSize) {
    float *i;
    float *j;
    float *k;

    if(treeNode->type == TREE_NODE_INTERNAL) { //is internal treeNode
        i = working + treeNode->traversalID->mid*3*alphabetSize;
        j = i + alphabetSize;
        k = j + alphabetSize;
        branchUp(treeNode->treeNodeX, working, i, subModels, totalP, alphabetSize);
        branchUp(treeNode->treeNodeY, working, j, subModels, totalP, alphabetSize);
        multiplyWV(i, j, k, alphabetSize);
        //uglyf("%i %f %f %f \n", treeNode->traversalID->mid, sumWV(k), totalP[treeNode->treeNodeX->traversalID->mid], totalP[treeNode->treeNodeY->traversalID->mid]);
        totalP[treeNode->traversalID->mid] = LOG(sumWV(k, alphabetSize)) + totalP[treeNode->treeNodeX->traversalID->mid] + totalP[treeNode->treeNodeY->traversalID->mid];
        normaliseWV(k, k, alphabetSize);
        return k;
    }
    assert(treeNode->type != TREE_NODE_DELETE);
    assert((treeNode->type & TREE_NODE_EFFECTIVELY_SILENT) == 0); //strong than above
    totalP[treeNode->traversalID->mid] = 0.0f;
    //l = sumWV(treeNode->wV);
    //assert(l < 1.01f);
    //assert(l > 0.99f);
    return treeNode->wV;
}

void downPass(struct TreeNode *treeNode, float *ancestorProb, float *working, float *results, struct SubModel **subModels, float *totalP, float totalAncP, float globalTotalP,  int32_t alphabetSize);

void branchDown(struct TreeNode *treeNode, float *ancestorProbs, float *working, float *results, struct SubModel **subModels, float *totalP, float totalAncP, float globalTotalP, int32_t alphabetSize) {
    transformWVByDistance(ancestorProbs, subModels[treeNode->traversalID->mid]->forward, ancestorProbs, alphabetSize);
    downPass(treeNode, ancestorProbs, working, results, subModels, totalP, totalAncP, globalTotalP, alphabetSize);
}

void downPass(struct TreeNode *treeNode, float *ancestorProbs, float *working, float *results,
			  struct SubModel **subModels, float *totalP, float totalAncP, float globalTotalP, int32_t alphabetSize) {
    float *i;
    float *j;
    float *k;
    float l;
    float *m;

    m = results + treeNode->traversalID->mid*alphabetSize;
    if (treeNode->type == TREE_NODE_INTERNAL) { //is internal treeNode
        i = working + treeNode->traversalID->mid*3*alphabetSize;
        j = i + alphabetSize;
        k = j + alphabetSize;
        multiplyWV(ancestorProbs, k, m, alphabetSize);
        l = totalP[treeNode->traversalID->mid] + totalAncP + LOG(sumWV(m, alphabetSize));
        normaliseWV(m, m, alphabetSize);
        //uglyf("%i %f %f \n", treeNode->traversalID->mid, l, globalTotalP);
        assert(l < globalTotalP + 0.01f);
        assert(l > globalTotalP - 0.01f);

        if(treeNode->treeNodeX->type != TREE_NODE_DELETE) {
            multiplyWV(ancestorProbs, j, j, alphabetSize);
            l = LOG(sumWV(j, alphabetSize)) + totalAncP + totalP[treeNode->treeNodeY->traversalID->mid];
            normaliseWV(j, j, alphabetSize);
            branchDown(treeNode->treeNodeX, j, working, results, subModels, totalP, l, globalTotalP, alphabetSize);
        }
        if(treeNode->treeNodeY->type != TREE_NODE_DELETE) {
            multiplyWV(ancestorProbs, i, i, alphabetSize);
            l = LOG(sumWV(i, alphabetSize)) + totalAncP + totalP[treeNode->treeNodeX->traversalID->mid];
            normaliseWV(i, i, alphabetSize);
            branchDown(treeNode->treeNodeY, i, working, results, subModels, totalP, l, globalTotalP, alphabetSize);
        }
    }
    else {
        assert(treeNode->type != TREE_NODE_DELETE);
        assert((treeNode->type & TREE_NODE_EFFECTIVELY_SILENT) == 0); //stronger than above
        memcpy(m, treeNode->wV, sizeof(float)*alphabetSize);
        l = sumWV(m, alphabetSize);
        //uglyf("%i %f %f \n", treeNode->traversalID->mid, l, globalTotalP);
        //assert(l < 1.01f);
        //assert(l > 0.99f);
    }
}

float *felsensteins(struct TreeNode *treeNode, struct SubModel **subModels, int32_t nodeNumber, float *ancestorProbs) {
    static float *working;
    static int32_t workSize;
    float *results;
    int32_t i;
    int32_t j;
    float *totalP;
    float totalAncP;
    float globalTotalP;
    int32_t alphabetSize;
    float fA[SMALL_CHUNK_SIZE];
    float *fA2;

    alphabetSize = 4; //subModels[treeNode->traversalID->mid]->alphabetSize;
    assert(alphabetSize < SMALL_CHUNK_SIZE);

    i = nodeNumber*alphabetSize;

    working = arrayResize(working, &workSize, i*3, sizeof(float));
    results = st_malloc(sizeof(float)*i);
    totalP = st_malloc(sizeof(float)*nodeNumber);
    for(j=0; j<i; j++) { //this is 'gap value'
        results[j] = -1.0f;
    }
    if(treeNode->type == TREE_NODE_INSERT) {
        treeNode = treeNode->treeNodeX;
    }
    totalAncP = 0.0f;
    assert((treeNode->type & TREE_NODE_EFFECTIVELY_SILENT) == 0);
    //ancestorProbs = ((struct SubModel *)subModels[treeNode->traversalID->mid])->insertionDistribution;
    fA2 = upPass(treeNode, working, subModels, totalP, alphabetSize);
    multiplyWV(ancestorProbs, fA2, fA, alphabetSize);
    globalTotalP = totalP[treeNode->traversalID->mid] + LOG(sumWV(fA, alphabetSize));
    assert(globalTotalP > -1000000);
    assert(globalTotalP < 1000000);
    downPass(treeNode, ancestorProbs, working, results, subModels, totalP, totalAncP, globalTotalP, alphabetSize);
    free(totalP);
    /*for(i=0; i<nodeNumber; i++) {
    	fA2 = results + i*ALPHABET_SIZE;
    	for(j=0; j<ALPHABET_SIZE; j++) {
    		if(fA2[j] >= -0.5) {
    			goto end;
    		}
    	}
    }
    assert(FALSE);
    end:*/
    return results;
}

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//control scripts
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

struct SequenceGraph *align_BottomUpScript(struct BinaryTree *binaryTree, struct SequenceGraph **sequenceGraphs, struct Constraints ***constraints,
                                           char **nodeNames, /*struct SubModel **subModels,*/ struct CombinedTransitionModel **combinedTransitionModels,
                                           int32_t numberOfSamples, int32_t leftMostSeqNo, int32_t *leafNo, int32_t *treeStates) {
    struct TraversalID *traversalID;
    char *nodeName;
    struct CombinedTransitionModel *combinedTransitionModel;

    struct BinaryTree *binaryTreeX;
    struct TraversalID *traversalIDX;
    //struct SubModel *subModelX;
    char *nodeNameX;
    struct SequenceGraph *sequenceGraphX;

    struct BinaryTree *binaryTreeY;
    struct TraversalID *traversalIDY;
    //struct SubModel *subModelY;
    char *nodeNameY;
    struct SequenceGraph *sequenceGraphY;

    int32_t leafSeqNoX;
    int32_t leafSeqNoY;

    struct SequenceGraph *sequenceGraph;

    traversalID = binaryTree->traversalID;
    nodeName = nodeNames[traversalID->mid];
    st_logInfo("Starting recursion to create node : %s \n", nodeName);
    if (binaryTree->internal) {
        binaryTreeX = binaryTree->left;
        traversalIDX = binaryTreeX->traversalID;
        //subModelX = subModels[traversalIDX->mid];
        nodeNameX = nodeNames[traversalIDX->mid];

        binaryTreeY = binaryTree->right;
        traversalIDY = binaryTreeY->traversalID;
        //subModelY = subModels[traversalIDY->mid];
        nodeNameY = nodeNames[traversalIDY->mid];

        leafSeqNoX = *leafNo;
        sequenceGraphX = align_BottomUpScript(binaryTreeX, sequenceGraphs, constraints, nodeNames, //subModels,
        combinedTransitionModels, numberOfSamples, leftMostSeqNo, leafNo, treeStates);
        leafSeqNoX = *leafNo - leafSeqNoX;

        leafSeqNoY = *leafNo;
        sequenceGraphY = align_BottomUpScript(binaryTreeY, sequenceGraphs, constraints, nodeNames, //subModels,
        combinedTransitionModels, numberOfSamples, leftMostSeqNo + leafSeqNoX, leafNo, treeStates);
        leafSeqNoY = *leafNo - leafSeqNoY;

        st_logInfo("Node %s is internal \n", nodeName);
        st_logInfo("Child x sequence : %s \n", nodeNames[traversalIDX->mid]);
        st_logInfo("X branch has length %f \n", binaryTreeX->distance);
        st_logInfo("Child y sequence : %s \n", nodeNames[traversalIDY->mid]);
        st_logInfo("Y branch has length %f \n", binaryTreeY->distance);
        combinedTransitionModel = combinedTransitionModels[traversalID->mid];
        sequenceGraph = computeEdgeGraph(sequenceGraphX, sequenceGraphY,
                                        combinedTransitionModel, constraints, //subModels,
                                        traversalID, traversalIDX, traversalIDY,
                                        leftMostSeqNo, leafSeqNoX, leafSeqNoY,
                                        numberOfSamples, treeStates);
        //memory clean up of input graphs
        destructSequenceGraph(sequenceGraphX, TRUE);
        destructSequenceGraph(sequenceGraphY, TRUE);
        //end clean up
        return sequenceGraph;
    }
    st_logInfo("Node %s is a leaf sequence \n", nodeName);
    //sequenceGraph = convertSeqToSeqGraph(seqsIt.next(), traversalIDs[binaryTree])
    return sequenceGraphs[(*leafNo)++];
}
