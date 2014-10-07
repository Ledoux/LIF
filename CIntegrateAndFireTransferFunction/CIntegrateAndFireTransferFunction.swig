%module CIntegrateAndFireTransferFunction
%include "../../SwigConfigs.swig"
%{ 
	#define SWIG_FILE_WITH_INIT
	#include "../../CTool.h"
	#include "CIntegrateAndFireTransferFunction.h"

%}
// Include the header file with above prototypes
%include "../../CTool.h"
%include "./CIntegrateAndFireTransferFunction.h"
