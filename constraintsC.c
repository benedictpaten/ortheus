/*
 * Copyright (C) 2008-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>

#include "hashTableC.h"
#include "commonC.h"
#include "fastCMaths.h"

#include "constraintsC.h"

#define CONSTRAINT_UNDECIDED 3 //reserved value
#define CONSTRAINT_BASE_SIZE 1000
#define CONSTRAINT_HASHTABLE_BASE_SIZE 1000
#define STACK_BASE_SIZE 1000

#define CONSTRAINT_MIN INT64_MIN
#define CONSTRAINT_MAX INT64_MAX


/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//alignment constraINT_32 data structures
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

struct Constraints *constructConstraintsBackward() {
    struct Constraints *constraints;
    constraints = st_malloc(sizeof(struct Constraints));
    assert(CONSTRAINT_BASE_SIZE >= 2);
    constraints->xList = st_malloc(sizeof(int64_t)*CONSTRAINT_BASE_SIZE);
    constraints->yList = st_malloc(sizeof(int64_t)*CONSTRAINT_BASE_SIZE);
    constraints->constraintsList = st_malloc(sizeof(int64_t)*CONSTRAINT_BASE_SIZE);
    constraints->maxLength = CONSTRAINT_BASE_SIZE;
    constraints->length = 2;

    constraints->xList[0] = CONSTRAINT_MAX;
    constraints->yList[0] = CONSTRAINT_MAX;

    constraints->xList[1] = CONSTRAINT_MIN;
    constraints->yList[1] = CONSTRAINT_MIN;

    constraints->constraintsList[0] = CONSTRAINT_LESS_THAN;
    constraints->constraintsList[1] = CONSTRAINT_LESS_THAN;

    return constraints;
}

void destructConstraints(struct Constraints *constraints) {
    //data-structure for storing set of pairwise constraints
    free(constraints->xList);
    free(constraints->yList);
    free(constraints->constraintsList);
    free(constraints);
}

void appendConstraintBackwards(struct Constraints *constraints, int64_t x, int64_t y, int64_t type) {
    //add prime constraINT_32 to end of list of constraints
    assert(constraints->xList[constraints->length-2] > x);
    assert(constraints->yList[constraints->length-2] > y || (constraints->yList[constraints->length-2] == y &&
                                  constraints->constraintsList[constraints->length-2] == CONSTRAINT_LESS_THAN_OR_EQUAL &&
                                  type == CONSTRAINT_LESS_THAN));

    if(constraints->length == constraints->maxLength) {
        void *new;
        int64_t newSize;

        newSize = 2*(constraints->maxLength);

        new = memcpy(st_malloc(sizeof(int64_t)*newSize), constraints->xList, sizeof(int64_t)*(constraints->maxLength));
        free(constraints->xList);
        constraints->xList = new;

        new = memcpy(st_malloc(sizeof(int64_t)*newSize), constraints->yList, sizeof(int64_t)*(constraints->maxLength));
        free(constraints->yList);
        constraints->yList = new;

        new = memcpy(st_malloc(sizeof(int64_t)*newSize), constraints->constraintsList, sizeof(int64_t)*(constraints->maxLength));
        free(constraints->constraintsList);
        constraints->constraintsList = new;

        constraints->maxLength = newSize;
    }

    constraints->xList[constraints->length-1] = x;
    constraints->yList[constraints->length-1] = y;
    constraints->constraintsList[constraints->length-1] = type;

    constraints->xList[constraints->length] = CONSTRAINT_MIN;
    constraints->yList[constraints->length] = CONSTRAINT_MIN;
    constraints->constraintsList[constraints->length++] = CONSTRAINT_LESS_THAN;

}

static int bSearchComparatorX(int64_t *i, int64_t *j) {
    if (*i > *j) {
        return -1;
    }
    if (*i == *j) {
        return 0;
    }
    if (*i > *(j+1)) {
        return 0;
    }
    return 1;
}

void getXConstraint(struct Constraints *constraints, int64_t x, int64_t *xConstraint, int64_t *yConstraint, int64_t *constraintType) {
    int64_t *i;
    int64_t j;
    j = constraints->length-1;
    if(x <= constraints->xList[j]) {
        *xConstraint = constraints->xList[j];
        *yConstraint = constraints->yList[j];
        *constraintType = constraints->constraintsList[j];
        return;
    }
    i = bsearch(&x, constraints->xList, j, sizeof(int64_t),
                (int (*)(const void *, const void *))bSearchComparatorX);
    //get prime constraINT_32 for poINT_32 x
    //i = bisect.bisect_left(self.xList, x)
    *xConstraint = *i;
    *yConstraint = constraints->yList[i - constraints->xList];
    *constraintType = constraints->constraintsList[i - constraints->xList];
    //return Constraint(constraints->xList[i], constraints->yList[i], constraints->constraintTypes[i])
}

static int bSearchComparatorY(int64_t *i, int64_t *j) {
    if (*i < *j) {
        return 1;
    }
    if (*i == *j) {
        return 0;
    }
    if (*i < *(j-1)) {
        return 0;
    }
    return -1;
}

void getYConstraint(struct Constraints *constraints, int64_t y, int64_t *xConstraint, int64_t *yConstraint, int64_t *constraintType) {
    //inverse of getXConstraINT_32
    int64_t *i;

    if(y >= constraints->yList[0]) {
        *xConstraint = constraints->xList[0];
        *yConstraint = constraints->yList[0];
        *constraintType = constraints->constraintsList[0];
        return;
    }
    i = bsearch(&y, constraints->yList+1, constraints->length-1, sizeof(int64_t),
                (int (*)(const void *, const void *))bSearchComparatorY);
    //get prime constraINT_32 for poINT_32 x
    //i = bisect.bisect_left(self.xList, x)
    *xConstraint = constraints->xList[i - constraints->yList];
    *yConstraint = *i;
    *constraintType = constraints->constraintsList[i - constraints->yList];
    //return Constraint(constraints->xList[i], constraints->yList[i], constraints->constraintTypes[i])
}

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//core constraINT_32 building functions
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

struct hashtable **getEmptyConstraints(const int64_t seqNo, const int64_t boundaries, struct Chunks *intChunks) {
    int64_t i;
    struct hashtable **tables;

    tables = st_malloc(sizeof(struct hashtable *)*seqNo*seqNo);
    for(i=0; i<seqNo*seqNo; i++) {
        tables[i] = create_hashtable(CONSTRAINT_HASHTABLE_BASE_SIZE, hashtable_intHashKey, hashtable_intEqualKey, (void (*)(void *))destructInt, (void (*)(void *))destructInt);
    }
    if(boundaries) {
        for(i=0; i<seqNo*seqNo; i++) {
            hashtable_insert(tables[i], constructChunkInt(0, intChunks), constructChunkInt(0, intChunks)); //constructInt(0), constructInt(1));
        }
    }
    //return [ [ { 0:1 } for j in xrange(0, seqNo) ] for i in xrange(0, seqNo) ]
    return tables;
}

void addConstraint(struct hashtable **constraintLists, int64_t *seqLengths, int64_t seqNo, int64_t seqX, int64_t seqY, int64_t x, int64_t y, struct Chunks *intChunks) {
    struct hashtable *constraints;
    int64_t *i;
    constraints = constraintLists[seqX * seqNo + seqY];
    if(y <= seqLengths[seqY]+1 && x >= 0) {
        //if not constraints.has_key(x):
        if((i = hashtable_search(constraints, &x)) == NULL) {
            //printf(" hello %" PRIi64 " %" PRIi64 " %" PRIi64 " %" PRIi64 " \n", seqX, seqY, x, y);
            hashtable_insert(constraints, constructChunkInt(x, intChunks), constructChunkInt(y, intChunks)); //constructInt(x), constructInt(y));
        }
        else {
            if(*i > y) {
            //if(*i < y) {
                //printf(" changing of the constraINT_32 %" PRIi64 " %" PRIi64 " %" PRIi64 " %" PRIi64 " % \n", seqX, seqY, x, y, i);
                *i = y;
            }
            //assert(*i <= y);
        }
    }
}

void addLessThanConstraints(int64_t *list, int64_t listLength, int64_t *indices, int64_t relaxValue,
                            struct hashtable **constraints, int64_t *seqLengths, int64_t seqNo, struct Chunks *intChunks) {
    int64_t i;
    int64_t j;
    int64_t seqX;
    int64_t seqY;
    int64_t posX;
    int64_t posY;
    for(i=0; i < listLength; i++) {
        //for j in xrange(i+1, len(list)):
        for(j=i+1; j<listLength; j++) {
            seqX = list[i];
            seqY = list[j];
            posX = indices[seqX];
            posY = indices[seqY];

            //posX = indices[seqX]-relaxValue;
            //posY = indices[seqY]+relaxValue;
            //addConstraint(constraints, seqLengths, seqNo, seqX, seqY, posX, posY+1, intChunks);
            //addConstraint(constraints, seqLengths, seqNo, seqX, seqY, posX-1, posY, intChunks);
            addConstraint(constraints, seqLengths, seqNo, seqX, seqY, posX, posY+1+relaxValue, intChunks);
            addConstraint(constraints, seqLengths, seqNo, seqX, seqY, posX-1-relaxValue, posY, intChunks);

            //posX = indices[seqX]+relaxValue;
            //posY = indices[seqY]-relaxValue;
            //addConstraint(constraints, seqLengths, seqNo, seqY, seqX, posY, posX+1, intChunks);
            //addConstraint(constraints, seqLengths, seqNo, seqY, seqX, posY-1, posX, intChunks);
            addConstraint(constraints, seqLengths, seqNo, seqY, seqX, posY, posX+1+relaxValue, intChunks);
            addConstraint(constraints, seqLengths, seqNo, seqY, seqX, posY-1-relaxValue, posX, intChunks);
        }
    }
}

void addMiddleConstraints(int64_t *nonGaps, int64_t nonGapsLength,
                          int64_t *gaps, int64_t gapsLength,
                          int64_t *indices, int64_t relaxValue,
                          struct hashtable **constraints, int64_t *seqLengths, int64_t seqNo, struct Chunks *intChunks) {
    int64_t i;
    int64_t j;
    int64_t seqX;
    int64_t seqY;
    int64_t posX;
    int64_t posY;

    for(i=0; i < nonGapsLength; i++) {
        //for j in xrange(i+1, len(nonGaps)):
        for(j=i+1; j<gapsLength; j++) {
            seqX = nonGaps[i];
            seqY = gaps[j];
            posX = indices[seqX];
            posY = indices[seqY];

            //posX = indices[seqX]-relaxValue;
            //posY = indices[seqY]+relaxValue;
            //addConstraint(constraints, seqLengths, seqNo, seqX, seqY, posX, posY+1, intChunks);
            addConstraint(constraints, seqLengths, seqNo, seqX, seqY, posX, posY+1+relaxValue, intChunks);

            //posX = indices[seqX]+relaxValue;
            //posY = indices[seqY]-relaxValue;
            //addConstraint(constraints, seqLengths, seqNo, seqY, seqX, posY, posX, intChunks);
            addConstraint(constraints, seqLengths, seqNo, seqY, seqX, posY-relaxValue, posX, intChunks);
        }
    }
}

struct hashtable **convertAlignmentToInputConstraints(char **alignment, int64_t alignmentLength, int64_t seqNo, int64_t *seqLengths,
                                                     int64_t relaxValue, char GAP, struct Chunks *intChunks) {
    //reads in a given alignment and converts it to a map of constraints
    struct hashtable **constraints;
    int64_t *indices;
    int64_t i;
    int64_t j;
    int64_t k;
    int64_t l;
    int64_t *list;
    int64_t *list2;

    constraints = getEmptyConstraints(seqNo, TRUE, intChunks);
    indices = st_malloc(sizeof(int64_t)*seqNo); //[0]*seqNo
    for(i=0; i<seqNo; i++) {
        indices[i] = 0;
    }
    list = st_malloc(sizeof(int64_t)*seqNo);
    list2 = st_malloc(sizeof(int64_t)*seqNo);
    //while((column = alignment()) != NULL) {
    for(i=0; i<alignmentLength; i++) {
        j=0;
        k=0;
        for(l=0; l<seqNo; l++) {
            //if column[columnIndex] != substitution.GAP:
            if(alignment[l][i] != GAP) {
                list[j++] = l;
                indices[l]++;
            }
            else {
                list2[k++] = l;
            }
        }
        addLessThanConstraints(list, j, indices, relaxValue, constraints, seqLengths, seqNo, intChunks);
        addMiddleConstraints(list, j, list2, k, indices, relaxValue, constraints, seqLengths, seqNo, intChunks);
    }
    for(i=0; i<seqNo; i++) {
        assert(indices[i] == seqLengths[i]);
    }
    //assert [ i for i in indices ] == seqLengths
    free(indices);
    free(list);
    free(list2);
    return constraints;
}

int64_t isAligned(int64_t *nonGaps, int64_t nonGapsLength, char **alignment, int64_t columnIndex, char GAP) {
    int64_t i;

    for(i=0; i<nonGapsLength; i++) {
        if(alignment[nonGaps[i]][columnIndex] != GAP) {
            return TRUE;
        }
    }
    return FALSE;
}

int64_t periodIsZero_Iter(int64_t *nonGaps, int64_t nonGapsLength, int64_t *gaps, int64_t gapsLength, int64_t from, int64_t to, int64_t change, char **alignment, char GAP, int64_t *l3) {
    int64_t i;
    int64_t j;
    int64_t lLength;
    int64_t l2Length;
    int64_t l3Length;

    static int64_t *l;
    static int64_t *l2;
    static int64_t lSize;
    static int64_t l2Size;

    l = arrayResize(l, &lSize, nonGapsLength + gapsLength, sizeof(int64_t));
    l2 = arrayResize(l2, &l2Size, nonGapsLength + gapsLength, sizeof(int64_t));
    lLength = nonGapsLength;
    l2Length = gapsLength;
    l3Length = 0;
    memcpy(l, nonGaps, nonGapsLength*sizeof(int64_t));
    memcpy(l2, gaps, gapsLength*sizeof(int64_t));

    for(i=from; i != to && l2Length != 0; i += change) {
        for(j=0; j<l2Length;) {
            if(alignment[l2[j]][i] != GAP) {
                if(isAligned(l, lLength, alignment, i, GAP)) {
                    l[lLength++] = l2[j];
                }
                else {
                    l3[l3Length++] = l2[j];
                }
                memmove(l2 + j, l2 + j + 1, (l2Size-(j+1))*sizeof(int64_t));
                l2Length--;
            }
            else {
                j++;
            }
        }
    }
    return l3Length;
}

struct BinaryTree *inSingleClade(struct BinaryTree *binaryTree, int64_t *gapSeqs, int64_t gapsLength) {
    struct BinaryTree *binarySubTree;
    if(binaryTree == NULL) {
        return NULL;
    }
    if(leftMostLeafNo(binaryTree->traversalID) <= gapSeqs[0] &&
       rightMostLeafNo(binaryTree->traversalID) >= gapSeqs[gapsLength-1]) {
        binarySubTree = inSingleClade(binaryTree->left, gapSeqs, gapsLength);
        if(binarySubTree != NULL) {
            return binarySubTree;
        }
        binarySubTree = inSingleClade(binaryTree->right, gapSeqs, gapsLength);
        if(binarySubTree != NULL) {
            return binarySubTree;
        }
        return binaryTree;
    }
    else {
        return NULL;
    }
}

int64_t merge(int64_t *l, int64_t *l2, int64_t lLength, int64_t l2Length) {
    int64_t i;
    int64_t j;
    int64_t k;

    for(i=0; i<l2Length;) {
        j = l2[i];
        for(k=0; k<lLength; k++) {
            if(j <= l[k]) {
                if(j < l[k]) {
                    memmove(l+k+1, l+k, (lLength-k)*sizeof(int64_t));
                    l[k] = j;
                    lLength++;
                }
                goto outer;
            }
        }
        l[k] = j;
        lLength++;
        outer:
        i++;
    }
    return lLength;
}

int64_t periodIsZero(int64_t *nonGaps, int64_t nonGapsLength, int64_t *gapSeqs, int64_t gapsLength, int64_t columnIndex, char **alignment, int64_t alignmentLength, char GAP, struct BinaryTree *binaryTree) {
    int64_t lLength;
    int64_t l2Length;

    static int64_t *l;
    static int64_t *l2;
    static int64_t lSize;
    static int64_t l2Size;

    l = arrayResize(l, &lSize, gapsLength, sizeof(int64_t));
    l2 = arrayResize(l2, &l2Size, gapsLength, sizeof(int64_t));

    lLength = periodIsZero_Iter(nonGaps, nonGapsLength, gapSeqs, gapsLength, columnIndex-1, -1, -1, alignment, GAP, l);
    l2Length = periodIsZero_Iter(nonGaps, nonGapsLength, gapSeqs, gapsLength, columnIndex+1, alignmentLength, 1, alignment, GAP, l2);

    merge(l, l, 0, lLength); //sorts l
    lLength = merge(l, l2, lLength, l2Length);
    if(lLength == 0) {
        return TRUE;
    }
    binaryTree = inSingleClade(binaryTree, l, lLength);
    return leafNoInSubtree(binaryTree->traversalID) == lLength;
}

void convertAlignmentToPhylogeneticInputConstraints_Recursion(char **alignment, int64_t alignmentLength, int64_t columnIndex, struct BinaryTree *binaryTree,
                                                              int64_t *indices, int64_t relaxValue, struct hashtable **constraints, int64_t *seqLengths, int64_t seqNo,
                                                              struct Chunks *intChunks, char GAP) {
    int64_t i;
    int64_t j;
    int64_t k;
    int64_t lLength;
    int64_t l2Length;

    static int64_t *l;
    static int64_t *l2;
    static int64_t lSize;
    static int64_t l2Size;

    if(binaryTree->internal) {
        i = leftMostLeafNo(binaryTree->traversalID);
        j = rightMostLeafNo(binaryTree->traversalID)+1;
        l = arrayResize(l, &lSize, j-i, sizeof(int64_t));
        l2 = arrayResize(l2, &l2Size, j-i, sizeof(int64_t));
        lLength = 0;
        l2Length = 0;

        for(k=i; k<j; k++) {
            if (alignment[k][columnIndex] == GAP) {
                l[lLength++] = k;
            }
            else {
                l2[l2Length++] = k;
            }
        }
        if(!periodIsZero(l2, l2Length, l, lLength, columnIndex, alignment, alignmentLength, GAP, binaryTree)) {
            convertAlignmentToPhylogeneticInputConstraints_Recursion(alignment, alignmentLength, columnIndex, binaryTree->left,
                                                             indices, relaxValue, constraints, seqLengths, seqNo, intChunks, GAP);
            convertAlignmentToPhylogeneticInputConstraints_Recursion(alignment, alignmentLength, columnIndex, binaryTree->right,
                                                             indices, relaxValue, constraints, seqLengths, seqNo, intChunks, GAP);
            return;
        }
        //success, can add constraints and return
        addLessThanConstraints(l2, l2Length, indices, relaxValue, constraints, seqLengths, seqNo, intChunks);
    }
}

struct hashtable **convertAlignmentToPhylogeneticInputConstraints(char **alignment, int64_t alignmentLength, int64_t seqNo, int64_t *seqLengths,
                                                   int64_t relaxValue, char GAP, struct Chunks *intChunks, struct BinaryTree *binaryTree) {
    struct hashtable **constraints;
    int64_t *indices;
    int64_t i;
    int64_t j;

    //reads in a given alignment and converts it to a map of constraints
    constraints =  getEmptyConstraints(seqNo, TRUE, intChunks);
    indices = st_calloc(seqNo, sizeof(int64_t)); //[0]*seqNo

    for(i=0; i<alignmentLength; i++) {
        for(j=0; j<seqNo; j++) {
            if(alignment[j][i] != GAP) {
                indices[j]++;
            }
        }
        convertAlignmentToPhylogeneticInputConstraints_Recursion(alignment, alignmentLength, i, binaryTree,
                                                                 indices, relaxValue, constraints, seqLengths, seqNo, intChunks, GAP);
    }
    for(i=0; i<seqNo; i++) {
        assert(indices[i] == seqLengths[i]);
    }
    free(indices);
    return constraints;
}

static int64_t selectedSeq;
static struct hashtable **lessThanConstraints;
static struct hashtable **lessThanOrEqualConstraints;
static int64_t seqNo;
static int64_t *seqLengths;

static struct Constraints **primeConstraints;
static int64_t *prime;
static int64_t *constraintType;
static int64_t *list;
static int64_t listIndex;
static int64_t **vertices;

static int64_t *stack;
static int64_t stackLength;
static int64_t stackMaxLength;

void searchFrom(int64_t vertexSeq, int64_t vertexPos, int64_t type, struct hashtable **constraintsLists) {
    int64_t seq;
    struct hashtable *constraints;
    int64_t *i;

    //for seq in xrange(seqNo-1, -1, -1):
    for(seq=seqNo-1; seq>=0; seq--) {
        if(seq != vertexSeq) {
            constraints = constraintsLists[vertexSeq*seqNo + seq];
            //if constraints.has_key(vertexPos):
            if ((i = hashtable_search(constraints, &vertexPos)) != NULL) {
                stack = arrayPrepareAppend(stack, &stackMaxLength, stackLength+3, sizeof(int64_t));
                stack[stackLength++] = type;
                stack[stackLength++] = *i;
                stack[stackLength++] = seq;
                //stack.append((seq, i, type));
            }
        }
    }
}

void search(int64_t vertexSeq, int64_t vertexPos, int64_t type) {
    int64_t vertexMark;

    while(TRUE) {
        vertexMark = vertices[vertexSeq][vertexPos];
        if (vertexMark == CONSTRAINT_UNDECIDED ||
        (vertexMark == CONSTRAINT_LESS_THAN_OR_EQUAL &&
         type == CONSTRAINT_LESS_THAN)) {
            vertices[vertexSeq][vertexPos] = type;
            if (prime[vertexSeq] >= vertexPos) {
                if (constraintType[vertexSeq] == CONSTRAINT_UNDECIDED && vertexSeq != selectedSeq) {
                    list[listIndex++] = vertexSeq;
                }
                prime[vertexSeq] = vertexPos;
                constraintType[vertexSeq] = type;
            }
            searchFrom(vertexSeq, vertexPos, type, lessThanOrEqualConstraints);
            searchFrom(vertexSeq, vertexPos, CONSTRAINT_LESS_THAN, lessThanConstraints);
            if(vertexPos+1 < seqLengths[vertexSeq]) {
                vertexPos++;
                type = CONSTRAINT_LESS_THAN;
                continue;
            }
        }
        if(stackLength == 0) {
            break;
        }
        vertexSeq = stack[--stackLength];
        vertexPos = stack[--stackLength];
        type = stack[--stackLength];
    }
}

void markup(int64_t vertexSeq, int64_t vertexPos) {
    struct hashtable *constraints;
    int64_t vertexMark;
    int64_t seq;
    int64_t *i;

    while(TRUE) {
        vertexMark = vertices[vertexSeq][vertexPos];
        if(vertexMark != CONSTRAINT_LESS_THAN) {
            vertices[vertexSeq][vertexPos] = CONSTRAINT_LESS_THAN;
            for(seq=seqNo-1; seq >= 0; seq--) {
                constraints = lessThanOrEqualConstraints[vertexSeq * seqNo + seq];
                if((i = hashtable_search(constraints, &vertexPos)) != NULL) { //constraints.has_key(vertexPos)) {
                    //i = constraints[vertexPos];
                    stack = arrayPrepareAppend(stack, &stackMaxLength, stackLength+2, sizeof(int64_t));
                    stack[stackLength++] = *i;
                    stack[stackLength++] = seq;
                    //stack.append((seq, i))
                }
            }
        }
        if(stackLength == 0) {
            break;
        }
        vertexSeq = stack[--stackLength];
        vertexPos = stack[--stackLength];
        //vertexSeq, vertexPos = stack.pop()
    }
}

struct Constraints **buildConstraints(struct hashtable **lessThanConstraintsA, struct hashtable **lessThanOrEqualConstraintsA, int64_t selectedSeqA, int64_t seqNoA, int64_t *seqLengthsA) {
    //computes set of prime constraints from set of input constraints
    //
    //exact copied, non-optimised implementation of Myers and Millers prime constraINT_32 algorithm
    //see page 14 of Progressive Multiple Alignment with Constraints, Myers et al.

    int64_t i;
    int64_t j;
    int64_t k;
    int64_t vertexPos;
    int64_t seq;
    int64_t *temp;

    selectedSeq = selectedSeqA;
    lessThanConstraints = lessThanConstraintsA;
    lessThanOrEqualConstraints = lessThanOrEqualConstraintsA;
    seqNo = seqNoA;
    seqLengths = seqLengthsA;

    primeConstraints = st_malloc(sizeof(void *)*seqNo); ///[ [] for i in xrange(0, seqNo) ]
    for(i=0; i<seqNo; i++) {
        primeConstraints[i] = constructConstraintsBackward();
    }

    prime = st_malloc(sizeof(int64_t)*seqNo); //[sys.maxint]*seqNo
    for(i=0; i<seqNo; i++) {
        prime[i] = INT64_MAX;
    }

    constraintType = st_malloc(sizeof(int64_t)*seqNo); //[CONSTRAINT_UNDECIDED]*seqNo
    for(i=0; i<seqNo; i++) {
        constraintType[i] = CONSTRAINT_UNDECIDED;
    }

    list = st_malloc(sizeof(int64_t)*seqNo); //[]

    vertices = st_malloc(sizeof(int64_t *)*seqNo);
    for(i=0; i<seqNo; i++) {
        k = seqLengths[i];
        temp = st_malloc(sizeof(int64_t)*k);
        vertices[i] = temp;
        for(j=0; j<k; j++) {
            temp[j] = CONSTRAINT_UNDECIDED;
        }
    }

    stackMaxLength = STACK_BASE_SIZE;
    stack = st_malloc(sizeof(int64_t)*stackMaxLength);
    stackLength = 0;
    //vertices = [ [ CONSTRAINT_UNDECIDED for pos in xrange(0, seqLengths[seq])]
    //            for seq in xrange(0, seqNo) ]
    //for vertexPos in xrange(seqLengths[selectedSeq]-1, -1, -1):
    for(vertexPos=seqLengths[selectedSeq]-1; vertexPos>=0; vertexPos--) {
        listIndex = 0;
        search(selectedSeq, vertexPos, CONSTRAINT_LESS_THAN_OR_EQUAL);
        assert(stackLength == 0);
        //for seq in list:
        for(i=0; i<listIndex; i++) {
            seq = list[i];
            appendConstraintBackwards(primeConstraints[seq], vertexPos, prime[seq], constraintType[seq]);
            constraintType[seq] = CONSTRAINT_UNDECIDED;
        }
        markup(selectedSeq, vertexPos);
        assert(stackLength == 0);
    }

    //memory clean up
    free(prime);

    free(constraintType);

    free(list);
    listIndex = INT64_MAX;

    for(i=0; i<seqNo; i++) {
        free(vertices[i]);
    }
    free(vertices);

    free(stack);
    stackLength = INT64_MAX;
    stackMaxLength = INT64_MAX;

    selectedSeq = INT64_MAX;
    lessThanConstraints = NULL;
    lessThanOrEqualConstraints = NULL;
    seqNo = INT64_MAX;
    seqLengths = NULL;
    //end clean up

    return primeConstraints;
}

struct Constraints ***buildAllConstraints_StartFromMinusOne(struct hashtable **lessThanConstraints,
                                                            struct hashtable **lessThanOrEqualConstraints,
                                                            int64_t seqNo, int64_t *seqLengths) {
    //the basic constraINT_32 algorithm assumes an offset starting at 0, but we need a -1 offset, for the position
    //before the first position, this is what this dirty wrapper function achieves
    int64_t *temp;
    int64_t i;
    int64_t j;
    int64_t k;
    struct Constraints ***primeConstraintsMatrix;
    struct Constraints *primeConstraints;

    temp = st_malloc(sizeof(int64_t)*seqNo);
    for(i=0; i<seqNo; i++) {
        temp[i] = seqLengths[i]+2;
    }
    primeConstraintsMatrix = st_malloc(sizeof(void *)*seqNo);
    for(i=0; i<seqNo; i++) {
        primeConstraintsMatrix[i] = buildConstraints(lessThanConstraints, lessThanOrEqualConstraints, i, seqNo, temp);
    }
    //primeConstraints =  [ buildConstraints(lessThanConstraints, lessThanOrEqualConstraints,
    //                                       selectedSeq, seqNo, seqLengths) for selectedSeq in xrange(0, seqNo) ]
    //for seq1 in xrange(0, seqNo):
    for(i=0; i<seqNo; i++) {
        //for seq2 in xrange(0, seqNo):
        for(j=0; j<seqNo; j++) {
            primeConstraints = primeConstraintsMatrix[i][j];
            for(k=0; k<primeConstraints->length-1; k++) {
            //for i in xrange(1, len(constraints.xList)-1):
                primeConstraints->xList[k] -= 1;
                primeConstraints->yList[k] -= 1;
            }
        }
    }
    //memory clean up
    free(temp);
    //end clean up
    return primeConstraintsMatrix;
}

struct Constraints ***buildAllConstraints_FromAlignment(char **alignment, int64_t alignmentLength, int64_t seqNo, int64_t *seqLengths,
                                                        int64_t relaxValue, int64_t GAP, struct BinaryTree *binaryTree, int64_t totalConstraints) {
    int64_t i;
    struct hashtable **lessThanConstraints;
    struct hashtable **lessThanOrEqualConstraints;
    struct Constraints ***primeConstraintsMatrix;
    struct Chunks *intChunks;

    intChunks = constructChunks(MEDIUM_CHUNK_SIZE, sizeof(int64_t));
    if(alignment != NULL) {
        if(totalConstraints) {
            lessThanConstraints = convertAlignmentToInputConstraints(alignment, alignmentLength, seqNo, seqLengths, relaxValue, GAP, intChunks);
        }
        else {
            lessThanConstraints = convertAlignmentToPhylogeneticInputConstraints(alignment, alignmentLength, seqNo, seqLengths, relaxValue, GAP, intChunks, binaryTree);
        }
    }
    else {
        lessThanConstraints = getEmptyConstraints(seqNo, TRUE, intChunks);
    }
    lessThanOrEqualConstraints = getEmptyConstraints(seqNo, FALSE, intChunks);

    primeConstraintsMatrix = buildAllConstraints_StartFromMinusOne(lessThanConstraints, lessThanOrEqualConstraints, seqNo, seqLengths);

    //memory clean up
    for(i=0; i<seqNo*seqNo; i++) {
        hashtable_destroy(lessThanConstraints[i], FALSE, FALSE);
        hashtable_destroy(lessThanOrEqualConstraints[i], FALSE, FALSE);
    }
    free(lessThanConstraints);
    free(lessThanOrEqualConstraints);
    destructChunks(intChunks);
    //end memory clean up
    return primeConstraintsMatrix;
}
