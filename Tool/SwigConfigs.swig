//introduce some C/C++ type to deal with
%include "std_vector.i"
%include "std_string.i"
%include "std_map.i"
%include "complex.i"

//Instantiate templates used by CGlobalDynamic
namespace std 
{
	%template(IntVector) vector<int>;
	%template(DoubleVector) vector<double>;
	%template(DoubleVectorVector) vector< vector<double> >;
	%template(StringVector) vector<string>;
	%template(MapStringInt) map<string,int>;
	%template(MapStringString) map<string,string>;
	%template(MapStringDouble) map<string,double>;
	%template(MapStringComplexDouble) map<string,std::complex<double> >;
}

