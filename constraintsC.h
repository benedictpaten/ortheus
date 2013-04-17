/*
 * Copyright (C) 2008-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */
/*
 * Copyright (C) 2008-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#ifndef CONSTRAINTSC_H_
#define CONSTRAINTSC_H_

#include "fastCMaths.h"
#include "commonC.h"

#define CONSTRAINT_LESS_THAN_OR_EQUAL 1
#define CONSTRAINT_LESS_THAN 2

struct Constraints {
    int64_t *xList;
    int64_t *yList;
    int64_t *constraintsList;
    int64_t length;
    int64_t maxLength;
};

void destructConstraints(struct Constraints *constraints);

struct Constraints ***buildAllConstraints_FromAlignment(char **alignment, int64_t alignmentLength, int64_t seqNo, int64_t *seqLengths, int64_t relaxValue, int64_t gap, struct BinaryTree *binaryTree, int64_t totalConstraints);

void getXConstraint(struct Constraints *constraints, int64_t x, int64_t *xConstraint, int64_t *yConstraint, int64_t *constraintType);

void getYConstraint(struct Constraints *constraints, int64_t y, int64_t *xConstraint, int64_t *yConstraint, int64_t *constraintType);

#endif /*CONSTRAINTSC_H_*/
