binPath=${rootPath}bin
libPath=${rootPath}lib
#Modify this variable to set the location of sonLib
sonLibRootPath=${rootPath}../sonLib
sonLibPath=${sonLibRootPath}/lib

include  ${sonLibRootPath}/include.mk

cflags += -I ${sonLibPath}
basicLibs = ${sonLibPath}/sonLib.a ${sonLibPath}/cutest.a
basicLibsDependencies = ${sonLibPath}/sonLib.a ${sonLibPath}/cutest.a 

