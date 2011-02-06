rootPath = ./
include ./include.mk

all : ${binPath}/ortheus_core ${binPath}/Ortheus.py

#These are parameters for composing/compiling the graph

#Use this line to point at a graphviz tool
makeGraph = dot

#This parameter can be adjusted (>=1) to control the silent-state number vs. the
#number of transitions.
collapseCoefficient = 1000 

#The location of the model files used
modelPath = ./models
model = ${modelPath}/affineModel.sxpr
paramModel = ${modelPath}/affineModel.param
 
clean :
	rm -f ${binPath}/ortheus_core model.otra xyzModelC.c xyzModelC.h ${binPath}/Ortheus.py

${binPath}/Ortheus.py : Ortheus.py old/Nester.py old/Stitcher.py old/EstimateTree.py old/bioio.py old/tree.py old/misc.py
	cp Ortheus.py ${binPath}/Ortheus.py
	chmod +x ${binPath}/Ortheus.py

${binPath}/ortheus_core : *.c *.h xyzModelC.c xyzModelC.h ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/ortheus_core *.c ${basicLibs}

# stdin from </dev/null works around stray stdin read on OS/X that hangs backgroud
# jobs
xyzModelC.c xyzModelC.h : ${model} ${paramModel} TransducerComposer.py TransducerCompiler.py
	rm -f model.otra xyzModelC.c xyzModelC.h ${modelPath}/model.dot ${modelPath}/ortheusmodel.pdf
	python TransducerComposer.py ${model} temp.otra ${modelPath}/model.dot ${collapseCoefficient} </dev/null
	cat  ${paramModel} temp.otra > model.otra
	rm temp.otra
	python  TransducerCompiler.py model.otra xyzModelC.c xyzModelC.h </dev/null
	#Use this line if you want the pretty picture
	#${makeGraph} ${modelPath}/model.dot -Tpdf > ${modelPath}/model.pdf
  
tests :
	#Running python allTests.py
	PYTHONPATH=.. PATH=../../bin:$$PATH python allTests.py --testLength=SHORT --logDebug
