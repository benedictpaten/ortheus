#ifndef CONSTRAINTSC_H_
#define CONSTRAINTSC_H_

#include "fastCMaths.h"
#include "commonC.h"

#define CONSTRAINT_LESS_THAN_OR_EQUAL 1
#define CONSTRAINT_LESS_THAN 2

struct Constraints {
    int32_t *xList;
    int32_t *yList;
    int32_t *constraintsList;
    int32_t length;
    int32_t maxLength;
};

void destructConstraints(struct Constraints *constraints);

struct Constraints ***buildAllConstraints_FromAlignment(char **alignment, int32_t alignmentLength, int32_t seqNo, int32_t *seqLengths, int32_t relaxValue, int32_t gap, struct BinaryTree *binaryTree, int32_t totalConstraints);

void getXConstraint(struct Constraints *constraints, int32_t x, int32_t *xConstraint, int32_t *yConstraint, int32_t *constraintType);

void getYConstraint(struct Constraints *constraints, int32_t y, int32_t *xConstraint, int32_t *yConstraint, int32_t *constraintType);

#endif /*CONSTRAINTSC_H_*/
