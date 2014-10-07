#ifndef CTOOL_H
#define CTOOL_H
#include <map>
#include <vector>
#include <string>
#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <iterator>
#include <typeinfo>
#include <complex>

class CToolClass;

class CToolClass
{
	public:
	
		/*** DICT ATTRIBUTES ***/
		std::map<std::string, double> DoubleDict;
		std::map<std::string, int> IntDict;
		std::map<std::string,std::string> StringDict;
		
		/*** SET METHODS ***/
		
		template<class KeyClass, class ValueClass>
		void setDict(std::map<KeyClass,ValueClass>* Dict,std::map<KeyClass,ValueClass> _Dict)
		{
			typedef typename std::map<KeyClass,ValueClass>::const_iterator const_iterator;
			for (const_iterator it=_Dict.begin(); it!=_Dict.end(); ++it)
			{
				(*Dict)[it->first]=it->second;
			}
		};
		
		void setDicts(
						std::map<std::string, double> _DoubleDict,
						std::map<std::string, int> _IntDict,
						std::map<std::string,std::string> _StringDict
					);
		
};


/***COMPLEX NUMBER METHODS***/

std::complex<double> getComplex(double re, double im);


//std::complex<double> getComplex(std::complex<double> X);

/*
std::complex<double> operator+(std::complex<double> a,double b){ return a+getComplex(b,0);};
std::complex<double> operator+(double b,std::complex<double> a){ return a+getComplex(b,0);};
std::complex<double> operator-(std::complex<double> a,double b){ return a-getComplex(b,0);};
std::complex<double> operator-(double b,std::complex<double> a){ return getComplex(b,0)-a;};
std::complex<double> operator*(std::complex<double> a,double b){ return a*getComplex(b,0);};
std::complex<double> operator*(double b,std::complex<double> a){ return a*getComplex(b,0);};
std::complex<double> operator/(std::complex<double> a,double b){ return a/getComplex(b,0);};
std::complex<double> operator/(double b,std::complex<double> a){ return getComplex(b,0)/a;};
*/

/***writing and reading text file functions***/
FILE * openFile(std::string nameFile);
void closeFile(FILE * pFile);
void getWord(FILE * pFile,char * word);
int readUntilLineEnd(FILE * pFile);
int readUntilFileEnd(FILE * pFile);
std::vector<fpos_t> getLinesStart(FILE * pFile);
int goToLine(FILE * pFile, int nLine);
int goToLastLine(FILE * pFile);
int getNLines(std::string nameFile);

/***generating new seed number***/
int updateRanSeed(std::string directoryPath);
double ran1(long *idum);
double gasdev(long *idum);
double expdev(long *idum);
double gamdev(int ia, long *idum);

/***combinatory analysis***/
int factorial(int n);
int arrangement(int k,int n);
int combinaison(int k,int n);


/***sum function***/
template<class T>
T sumVec(std::vector<T> * vec)
{
	T sum=(T)(0);
	for(int i=0;i<(int)((*vec).size());i++)
	{
		sum+=(*vec)[i];
	}
	return sum;
	//return accumulate((*vec).begin(),(*vec).end(),0.);
}

template<class T>
void printVec(std::vector<T> * vec)
{
	for(int i=0;i<(int)(*vec).size();i++)
	{
		std::cout<<(*vec)[i]<< " ";
	}
	std::cout<<std::endl;
}

void printVecFloat(std::vector<float> * vec);
void printVecInt(std::vector<int> * vec);

template<class T>
void printArray(std::vector< std::vector<T> > * array)
{
	for(int i=0;i<(int)(*array).size();i++)
	{
		for(int j=0;j<(int)(*array)[i].size();j++)
		{
			std::cout<<(*array)[i][j]<< " ";
		}
		std::cout<<std::endl;
	}
};

void printArrayFloat(std::vector< std::vector<float> >  * vec);
void printArrayInt(std::vector< std::vector<int> >  * vec);


class isEqualValue{

public:
	isEqualValue(int _value=0) : value(_value) {}
	bool operator() (int x) const
	{ 
		return (x==value); 
	}

private:
	int value;
};


class isEven{

public:
	bool operator() (int x) const
	{ 
		return ((x % 2) == 0); 
	}
};

class isOdd{

public:
	bool operator() (int x) const
	{ 
		return ((x % 2) == 1); 
	}
};

/***build a functor that computes the mean***/
class meanClass 
{
public:
	meanClass(int _vecSize=0): meanValue(0.), vecSize(_vecSize) {}
	int meanValue;
	int vecSize;
	void operator()(double v) { meanValue += v/vecSize; }
};


double meanValue(std::vector<double> * vec);
double stdValue(std::vector<double> * vec);
std::vector<double> meanVec(std::vector<std::vector<double> > * vec);
std::vector<double> stdVec(std::vector<std::vector<double> > * vec);
std::vector< int> findAll_whereValue(std::vector<int> vec, int value);
std::vector< std::vector< int> > findAll_even(std::vector<int> vec);
std::vector< std::vector<int> > getPermutations(std::vector<int> vec, int begin, int end);
std::vector<int> getRowFromVecVec(std::vector< std::vector<int> >  vecVec, int row);
std::vector<int> getColFromVecVec(std::vector< std::vector<int> >  vecVec, int col);
std::vector< std::vector<int> > getRowsFromVecVec(std::vector< std::vector<int> >  vecVec, std::vector<int> rows);
std::vector< std::vector<int> > getVecVecFromVecVec(std::vector< std::vector<int> >  vecVec, std::vector<int> rows, std::vector<int> cols);
std::vector< std::vector<int> > getCombinaisons(std::vector<std::vector<int> > permu, int size);

/***look for the indx of a listed name***/
int whichIndx(std::string name, std::vector<std::string> names);

/***MATHEMATIC TOOLS ON C++ VECTORS***/
//template <typename T> 
//void vecAdd(std::vector<T> * A,std::vector<T> * B)
//{
//	if((*A).size()!=(*B).size()){printf("ERROR in vecAdd : vectors don't have the same length\n");}
//	for(int i=0;i<(int)(*A).size();i++){(*A)[i]+=(*B)[i];}
//};
//template <typename T> 
//void vecAdd(std::vector<T> * A,T b)
//{
//	for(int i=0;i<(int)(*A).size();i++){(*A)[i]+=b;}
//};

void vecSub(std::vector<double> * A,std::vector<double> * B);
void vecSub(std::vector<double> * A,double b);

void vecSub(std::vector<double> * A,std::vector<double> * B);
void vecMul(std::vector<double> * A,std::vector<double> * B);

std::vector<double> vecFromSub(std::vector<double> * A,std::vector<double> * B);
std::vector<double> operator-(std::vector<double>  A,std::vector<double> B);

void vecSub(std::vector<double> * A,double b);
void vecMul(std::vector<double> * A,double b);

void vecVecAdd(std::vector< std::vector<double> > * A, std::vector< std::vector<double> > * B);

std::vector<int> range(int start,int end,int step);
std::vector<double> arange(double start,double end,double step);

std::vector<double> vecOnes(int sizeVec);
std::vector<double> vecZeros(int sizeVec);
std::vector< std::vector<double> > vecVecOnes(std::vector<int> sizeVec);
std::vector< std::vector<double> > vecVecZeros(std::vector<int> sizeVec);




/*
* Program to demostrate the use of template class and associative array.
* By Arnav Mukhopadhyay
*/

template <class T>
class AssocArray
{
  private:
	  typedef struct _Data
	  {
		  T data;
		  std::string name;
	  } Data ;
	  std::vector<Data> stack;
  public:
	  long Size()
	  {
		  return stack.size();
	  }

	  bool IsItem(std::string name)
	  {
		  for(int i=0; i<Size(); i++)
		  {
			  if(stack[i].name == name)
				  return true;
		  }
		  return false;
	  }

	  bool AddItem(std::string name, T data)
	  {
		  if(IsItem(name))
			  return false;
		  Data d;
		  d.name = name;
		  d.data = data;
		  stack.push_back(d);
		  return true;
	  }

	  T& operator [] (std::string name)
	  {
		  for(int i=0; i<Size(); i++)
		  {
			  if(stack[i].name == name)
				  return stack[i].data;
		  }
		  long idx = Size();
		  Data d;
		  d.name = name;
		  stack.push_back(d);
		  return stack[idx].data;
	  }

	  std::string GetItemName(long index)
	  {
		  if(index<0)
			  index = 0;
		  for(int i=0; i<Size(); i++)
			  if(i == index)
				  return stack[i].name;
		  return "";
	  }

	  T& operator [] (long index)
	  {
		  if(index < 0)
			  index = 0;
		  for(int i=0; i<Size(); i++)
		  {
			  if(i == index)
				  return stack[i].data;
		  }
		  return stack[0].data;
	  } 

};


#endif