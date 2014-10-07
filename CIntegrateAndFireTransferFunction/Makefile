#######################
#
#	Building CIntegrateAndFireTransferFunction for a Python module
#
#######################

CXX=gcc
ARCH=-arch x86_64
CFLAGS=  -I/System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7

FILE_0=../Tool/CTool
FILE_1=CIntegrateAndFireTransferFunction

all: ${FILE_1}_wrap.cpp ${FILE_1}.o ${FILE_1}_wrap.o ${FILE_1}.so 

${FILE_1}_wrap.cpp :${FILE_1}.swig ${FILE_1}.h ${FILE_0}.h
	@echo "=======prepare ${FILE_1} module for python==="
	swig -c++ -python -o ${FILE_1}_wrap.cpp ${FILE_1}.swig
	@echo "\n"

${FILE_1}.o :${FILE_1}.cpp ${FILE_1}.h ${FILE_0}.h
	@echo "=======compile ${FILE_1}.cpp=========="
	$(CXX) $(ARCH) -c -fPIC ${FILE_1}.cpp -o ${FILE_1}.o $(CFLAGS)
	@echo "\n"

${FILE_1}_wrap.o : ${FILE_1}_wrap.cpp ${FILE_1}.h ${FILE_0}.h
	@echo "=======compile ${FILE_1}_wrap.cpp=========="
	$(CXX) $(ARCH) -c -fPIC ${FILE_1}_wrap.cpp -o ${FILE_1}_wrap.o $(CFLAGS) 
	@echo "\n"

${FILE_1}.so : ${FILE_1}.o ${FILE_0}.o ${FILE_1}_wrap.o ${FILE_1}.h ${FILE_0}.h
	@echo "=======linking objects=========="
	$(CXX) $(ARCH) -bundle -flat_namespace -undefined suppress -o _${FILE_1}.so ${FILE_1}.o ${FILE_0}.o ${FILE_1}_wrap.o $(CFLAGS)
	@echo "\n"

clean : 
	$(RM) ${FILE_1}_wrap.cpp ${FILE_1}_wrap.o _${FILE_1}.so ${FILE_1}.o