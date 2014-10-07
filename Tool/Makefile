#######################
#
#	Building CTool Class for a Python module
#
#######################

CXX=gcc
ARCH=-arch x86_64
CFLAGS= -I/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7
FILE_0=CTool

all: ${FILE_0}_wrap.cpp ${FILE_0}.o ${FILE_0}_wrap.o ${FILE_0}.so
	
${FILE_0}_wrap.cpp : ${FILE_0}.swig ${FILE_0}.h
	@echo "=====prepare ${FILE_0}_wrap module for python====="
	swig -c++ -python -o ${FILE_0}_wrap.cpp ${FILE_0}.swig
	@echo "\n"

${FILE_0}.o : ${FILE_0}.cpp ${FILE_0}.h
	@echo "=====compile ${FILE_0}.cpp =================="
	$(CXX) $(ARCH) -c -fPIC ${FILE_0}.cpp -o ${FILE_0}.o $(CFLAGS)
	@echo "\n"

${FILE_0}_wrap.o : ${FILE_0}_wrap.cpp ${FILE_0}.h
	@echo "=====compile ${FILE_0}_wrap.cpp =================="
	$(CXX) $(ARCH) -c -fPIC ${FILE_0}_wrap.cpp -o ${FILE_0}_wrap.o $(CFLAGS)
	@echo "\n"

${FILE_0}.so : ${FILE_0}.o ${FILE_0}_wrap.o ${FILE_0}.h
	@echo "=====link ${FILE_0}_wrap.o ${FILE_0}.o ============"
	$(CXX) $(ARCH) -bundle -flat_namespace -undefined suppress -o _${FILE_0}.so ${FILE_0}.o ${FILE_0}_wrap.o $(CFLAGS)
	@echo "\n"

clean : 
	$(RM) ${FILE_0}_wrap.cpp ${FILE_0}.o _${FILE_0}.so