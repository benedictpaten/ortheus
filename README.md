Ortheus
=======

Evolutionary transducer based method for aligning and reconstructing indel
histories (Ensembl Compara modifications)

This repository holds the necessary changes of [Benedict Paten's version](https://github.com/benedictpaten/ortheus)
to run the latest Ensembl Compara pipeline.
Ensembl is **not** the official maintainer of this software.
Pull-requests can still be submitted, but we will only accept them if they
can provide a benefit to Ensembl.

You can find more documentation on the UCSC page: (http://hgwdev.cse.ucsc.edu/~benedict/code/Ortheus.html)

The main new features are:
* new `-B` option to refer to a file that contains the Newick tree, instead
  of passing the Newick tree itself as a command-line argument
* new `-A` option to refer to a file that contains the paths of the Fasta
  files to align, instead of passing these as command-line arguments.

Significant bugs that are now fixed are:
* wrong branch lengths when rerooting the tree
* off-by-one gaps in the alignments
* _AttributeError_ failure in `calculateProbableRootOfGeneTree`

## Branches and tags

There is a single branch (master) where all the development goes.
The version number stated in the source code (0.5.0) is not maintained.

`ensembl_production_XX` tags are used to refer to the version used for the
production of Ensembl version XX. Due to deployment constraints, these tags
may not include the latest changes of the master branch.
Instead, we provide `ensembl_release_candidate_Y` tags, Y starting from 1, for
the "next" version we will deploy in production.

## Installation

Installing Ortheus.

(1) Download and install Pecan

(2) Download and install sonLib

(3) Download and install Semphy, if tree estimation is required

(4) Place the directory containing Ortheus on your python path, i.e.
PYTHONPATH=${PYTHONPATH}:FOO
where FOO/sonLib is the path to the base directory of Ortheus. 

(5) Add the bin directory to your path:
PATH=${PATH}:foo/ortheus/bin
where foo/ortheus is the path to the base directory of Ortheus

(6) Compile the C code:
Modify the include.mk file to point at where you installed sonLib. You need not do this if 
you install sonLib and ortheus in the same parent directory.
In Ortheus type 'make all' 

(7) Test the installation
Type 'make test' in the base directory if ortheus is installed in the same parent directory as sonLib, else type 'python allTests.py'

The tests will not verify the installation of Pecan or Semphy currently, but rather just test the core Ortheus algorithms.

See https://github.com/benedictpaten/ for different projects.

For usage instructions see doc/README
