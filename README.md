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
