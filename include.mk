#Location of sonLib
binPath=${rootPath}bin
libPath=${rootPath}lib
sonLibRootPath=${rootPath}../sonLib
sonLibPath=${sonLibRootPath}/lib

include  ${sonLibRootPath}/include.mk

cflags += -I ${sonLibPath}
basicLibs = ${sonLibPath}/sonLib.a ${sonLibPath}/cutest.a
basicLibsDependencies = ${sonLibPath}/sonLib.a ${sonLibPath}/cutest.a 

