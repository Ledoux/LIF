#ifndef CTOOL_CPP
#define  CTOOL_CPP
#include "CTool.h"
#include <Accelerate/Accelerate.h>
#include <vector>
#include <string>
#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>

using namespace std;


/***PYTHONIZE DICT ATTRIBUTES CONTAINER C++ CLASS LIKE***/

std::map<std::string,double> _DefaultStringDoubleMap;
std::map<std::string,int> _DefaultStringIntMap;
std::map<std::string,std::string> _DefaultStringStringMap;
void CToolClass::setDicts(
						  std::map<std::string,double> _DoubleDict=_DefaultStringDoubleMap,
						  std::map<std::string,int> _IntDict=_DefaultStringIntMap,
						  std::map<std::string,std::string> _StringDict=_DefaultStringStringMap
						  )
{
	setDict(&DoubleDict,_DoubleDict);
	setDict(&IntDict,_IntDict);
	setDict(&StringDict,_StringDict);
}

/***COMPLEX NUMBER METHODS***/


std::complex<double> getComplex(double re, double im)
{
  std::complex<double> out(re,im);
  return out;
}

/***TEXT FILE MANAGING FUNCTIONS***/

void getWord(FILE * pFile, char * word)
{
	int c;
	int i=0;
	while((c=getc(pFile))!=' ')
	{
		word[i]=c;
		i++;
	}
	word[i] = '\0';
}

/***open a File***/
FILE * openFile(string nameFile)
{
	FILE * pFile;
	char dataFileName[124];
	sprintf(dataFileName,"%s",nameFile.c_str());
	if((pFile=fopen(dataFileName,"r"))==NULL){printf("ERROR : No data called like %s.\n",dataFileName);exit(1);}
	return pFile;
}

/***close a File***/
void closeFile(FILE * pFile)
{
	fclose(pFile);
}

/***read until Line end and return number of data***/
int readUntilLineEnd(FILE * pFile)
{
	int c,n=0;

	/***read the end of line***/
	do
	{
		c=getc(pFile);
		if(c==32){n++;}
	}
	while(c != EOF&&c!=10);

	return n;
}

/***read until end and return the number of lines***/
int readUntilFileEnd(FILE * pFile)
{
	int c,n=0;

	/***read the end of line***/
	do
	{
		c=getc(pFile);
		/***count the lines***/
		if(c==10)
		{
			n++;
		}
		
	}
	while(c != EOF);

	
	return n;
}

//record the positions where a new line begins ***/
std::vector<fpos_t> getLinesStart(FILE * pFile)
{
	int c;
	std::vector<fpos_t> posLines;
	
	//be sure to be at the beginning of the file
	rewind(pFile);
	fpos_t posFile;
	fgetpos(pFile,&posFile);
	posLines.push_back(posFile);
	
	//scan until the end of the file
	while(c!=EOF)
	{
		c=getc(pFile);
		//if it's a begin of a line
		if(c==10)
		{
			fpos_t posFile;
			fgetpos(pFile,&posFile);
			posLines.push_back(posFile);
		}
	}
	
	//remove the last element that is pointer to the EOF
	posLines.pop_back();
	return posLines;
}

//read until the nLine and check that EOF is not reached before***/
int goToLine(FILE * pFile, int nLine)
{
	int c,n=0;
	while(c!=EOF&&n<nLine)
	{
		c=getc(pFile);
		if(c==10){n++;}	
	}
	if(nLine>0)
	{
		if(n<nLine){printf("WARNING : EOF reached at the %d Line before getting the %d Line\n",n,nLine);}
		return n;
	}
	else
	{
		return 0;
	}
}

/***go to the last line of a file***/
int goToLastLine(FILE * pFile)
{
	int nLines;
	nLines=readUntilFileEnd(pFile);
	rewind(pFile);
	goToLine(pFile,nLines-1);

	return nLines;
}

//give the number of lines of a file
int getNLines(string nameFile)
{
	FILE * pFile;
	int nLines;
	pFile=openFile(nameFile);
	nLines=readUntilFileEnd(pFile);
	fclose(pFile);
	return nLines;
}

/***RANDOM NUMBERS TOOLS***/

/***set an initial Seed number different for each simulation***/
int updateRanSeed(string directoryPath)
{
	FILE *pFileRandom;
	char bufferFileRandom[1000];
	int idumread=1200;
	sprintf(bufferFileRandom,"%sfileRandom.txt",directoryPath.c_str());
	if((pFileRandom=fopen(bufferFileRandom,"r"))==NULL)
	{
		printf("WARNING : %s not found\n",bufferFileRandom);
		pFileRandom=fopen(bufferFileRandom,"w");
		fprintf(pFileRandom,"%d",1);
		fclose(pFileRandom);
		pFileRandom=fopen(bufferFileRandom,"r");
	}
	fscanf(pFileRandom,"%d ",&idumread);
	fclose(pFileRandom);
	int idum=idumread+1;
	pFileRandom=fopen(bufferFileRandom,"w");
	fprintf(pFileRandom,"%d ",idum);
	fclose(pFileRandom);

	return idum;
}


/***ran1 numericl recipes, uniform pick point between 0 and 1***/

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

/***gasdev numericl recipes, gaussian pick point mean 0 and std 1***/

double gasdev(long *idum)
{
	double ran1(long *idum);
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;

	if  (iset == 0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}

/***expdev numericl recipes, exponential pick point lambda 1***/

double expdev(long *idum)
{
	double ran1(long *idum);
	double dum;

	do
		dum=ran1(idum);
	while (dum == 0.0);
	return -log(dum);
}

/***gammadev numericl recipes,gamma pick point  ia parameter***/

double gamdev(int ia, long *idum)
{
	double ran1(long *idum);
	void nrerror(char error_text[]);
	int j;
	double am,e,s,v1,v2,x,y;

	if (ia < 1) printf("Error in routine gamdev\n");
	if (ia < 6) {
		x=1.0;
		for (j=1;j<=ia;j++) x *= ran1(idum);
		x = -log(x);
	} else {
		do {
			do {
				do {
					v1=ran1(idum);
					v2=2.0*ran1(idum)-1.0;
				} while (v1*v1+v2*v2 > 1.0);
				y=v2/v1;
				am=ia-1;
				s=sqrt(2.0*am+1.0);
				x=s*y+am;
			} while (x <= 0.0);
			e=(1.0+y*y)*exp(am*log(x/am)-s*y);
		} while (ran1(idum) > e);
	}
	return x;
}

/***COMBINATORY ANALYSIS****/
int factorial(int n)
{
 	 return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}
int arrangement(int k,int n)
{
	return (k<=n)? factorial(n)/factorial(n-k):0;
}
int combinaison(int k,int n)
{
	return arrangement(k,n)/factorial(k);
}

/***MANIPULATING DATA VECTORS***/

vector<int> findAll_whereValue(vector<int> vec, int value)
{
	vector<int> vecSolu;
	for (std::vector<int>::iterator it = std::find_if(vec.begin(), vec.end(),isEqualValue(value));
	    it != vec.end();
	    it = std::find_if(++it, vec.end(), isEqualValue(value)))
	{
		size_t i = it - vec.begin();
		vecSolu.push_back((int) i);
	}
	return vecSolu;
}

vector< vector<int> > findAll_even(vector<int> vec)
{
	vector< vector<int> > vecSolu(2);
	for (std::vector<int>::iterator it = std::find_if(vec.begin(), vec.end(),isEven());
	    it != vec.end();
	    it = std::find_if(++it, vec.end(), isEven()))
	{
	 	vecSolu[0].push_back(*it);
		size_t i = it - vec.begin();
		vecSolu[1].push_back((int) i);
	}
	return vecSolu;
}


vector< vector<int> > getPermutations(vector<int> vec, int begin, int end)
{
	vector< vector<int> > vecVec;
	if(begin<end&&end<(int)vec.size())
	{
		int size=factorial(1+end-begin);
		vecVec.resize(size);
		int count=0;
		//using next_permutations from the std algorihtm library
		do
		{
			vecVec[count].resize(vec.size());
			vecVec[count]=vec;
			count++;
		}
		while (count<size&&next_permutation(vec.begin()+begin, vec.begin()+end+1)); 
	}
	else
	{
		printf("WARNING : begin >= end or end>= vec.size() \n");
	}
	return vecVec;
}

vector<int> getVecFromVec(vector<int> vec, vector<int> indx)
{
	vector<int> vecOut;
	for(int i=0;i<(int)indx.size();i++)
	{	
		if(indx[i]<(int)vec.size())
		{
			vecOut.push_back(vec[indx[i]]);
		}
		else
		{
			printf("WARNING : indx[i]>=vec.size\n");
		}
	}
	return vecOut;
}


vector<int> getRowFromVecVec(vector< vector<int> >  vecVec, int row)
{
	
	if(row<(int)vecVec.size())
	{	
		return vecVec[row];
	}
	else
	{
		vector<int> outNULL;
		printf("WARNING : row >=vecVec.size() \n");
		return outNULL;
	}
}

vector<int> getColFromVecVec(vector< vector<int> > vecVec, int col)
{
	vector<int> vec;
	for(int i=0;i<(int)vecVec.size();i++)
	{
		if(col<(int)vecVec[i].size())
		{
			vec.push_back(vecVec[i][col]);
		}
		else
		{
			printf("WARNING : col >=vecVec[%d].size() \n",i);
		}
	}
	return vec;
}


vector< vector<int> > getRowsFromVecVec(vector< vector<int> >  vecVec, vector<int> rows)
{
	vector< vector<int> > vecVecRows(rows.size());
	for(int i=0;i<(int)rows.size();i++)
	{
		vecVecRows[i].resize(vecVec[i].size());
		vecVecRows[i]=getRowFromVecVec(vecVec,rows[i]);
	}
	return vecVecRows;
}

vector< vector<int> > getVecVecFromVecVec(vector< vector<int> >  vecVec,vector<int> rows,vector<int> cols)
{

		vector< vector<int> > vecVecRows;
		vecVecRows=getRowsFromVecVec(vecVec,rows);

		vector< vector<int> > vecVecOut(rows.size());
		for(int i=0;i<(int)vecVecOut.size();i++)
		{
			vecVecOut[i].resize(cols.size());
			for(int j=0;j<(int)cols.size();j++)
			{
				if(cols[j]<(int)vecVecRows[i].size())
				{
					vecVecOut[i][j]=vecVecRows[i][cols[j]];
				}
				else
				{
					printf("WARNING : cols[%d]>=vecVecRows.size()\n",j);
				}
			}	
		}
		return vecVecOut;
}

vector< vector<int> > getCombinaisons(vector<vector<int> > permu, int size)
{
	vector< vector<int> > vecVec;
	int nCol=permu[0].size();
	if(size<=nCol)
	{
		int nCombi=arrangement(size,nCol);
		vector<int> rows(nCombi);
		rows[0]=0;
		for(int i=1;i<nCombi;i++)
		{
			rows[i]=rows[i-1]+factorial(nCol-size);
		}
		vector<int> cols(size);
		for(int i=0;i<size;i++)
		{
			cols[i]=i;
		}
		vecVec=getVecVecFromVecVec(permu,rows,cols);	
	}
	else
	{
		vecVec.resize(0);
		vecVec[0].resize(0);
		printf("WARNING size>permu[0].size()\n");
	}
	return vecVec;
}




/***look for the indx of a listed name***/
int whichIndx(string name, vector<string> names)
{
	int i=0;
	while(i<(int)names.size()&&names[i]!=name){i++;}
	if(i==(int)names.size()){printf("WARNING : name not found in the list\n");} 
	return i;
} 


/***MATHEMATIC TOOLS ON C++ VECTORS***/

void vecAdd(vector<double> * A,vector<double> * B)
{
	if((*A).size()!=(*B).size()){printf("ERROR in vecSub : vectors don't have the same length\n");}
	for(int i=0;i<(int)(*A).size();i++){(*A)[i]+=(*B)[i];}
}
void vecSub(vector<double> * A,vector<double> * B)
{
	if((*A).size()!=(*B).size()){printf("ERROR in vecSub : vectors don't have the same length\n");}
	for(int i=0;i<(int)(*A).size();i++){(*A)[i]-=(*B)[i];}
}
void vecMul(vector<double> * A,vector<double> * B)
{
	if((*A).size()!=(*B).size()){printf("ERROR in vecSub : vectors don't have the same length\n");}
	for(int i=0;i<(int)(*A).size();i++){(*A)[i]*=(*B)[i];}
}

void vecAdd(vector<double> * A,double b)
{
	for(int i=0;i<(int)(*A).size();i++){(*A)[i]+=b;}
}
void vecSub(vector<double> * A,double b)
{
	for(int i=0;i<(int)(*A).size();i++){(*A)[i]-=b;}
}
void vecMul(vector<double> * A,double b)
{
	for(int i=0;i<(int)(*A).size();i++){(*A)[i]*=b;}
}

vector<double> vecFromSub(vector<double> * A,vector<double> * B)
{
	if((*A).size()!=(*B).size()){printf("ERROR in vecSub : vectors don't have the same length\n");}
	vector<double> C((*A).size());
	for(int i=0;i<(int)(*A).size();i++){C[i]=(*A)[i]-(*B)[i];}
	return C;
}

vector<double> operator-(vector<double>  A,vector<double>  B)
{
	if(A.size()!=B.size()){printf("ERROR in vecSub : vectors don't have the same length\n");}
	vector<double> C(A.size());
	for(int i=0;i<(int)A.size();i++){C[i]=A[i]-B[i];}
	return C;
}

vector<int> range(int start,int end,int step)
{
	vector<int> output(0);
	if(start>end && step>0.)
	{
		printf("WARNING : start>end && step>0. in range\n");
		return output;
	}
	else if(start<end && step<0.)
	{
		printf("WARNING : start<end && step<0. in range\n");
		return output;
	}
	else
	{
		if(start<end)
		{
			output.push_back(start);
			while(output[output.size()-1]<end)
			{
				output.push_back(output[output.size()-1]+step);
			}
		}
		else
		{
			output.push_back(start);
			while(output[output.size()-1]>end)
			{
				output.push_back(output[output.size()-1]+step);
			}
		}
		return output;
	}
}

vector<double> arange(double start,double end,double step)
{
	vector<double> output(0);
	if(start>end && step>0.)
	{
		printf("WARNING : start>end && step>0. in arange\n");
		return output;
	}
	else if(start<end && step<0.)
	{
		printf("WARNING : start<end && step<0. in arange\n");
		return output;
	}
	else
	{
		if(start<end)
		{
			output.push_back(start);
			while(output[output.size()-1]<end)
			{
				output.push_back(output[output.size()-1]+step);
			}
		}
		else
		{
			output.push_back(start);
			while(output[output.size()-1]>end)
			{
				output.push_back(output[output.size()-1]+step);
			}
		}
		return output;
	}
}

void vecVecAdd(vector< vector<double> > * A, vector< vector<double> > * B)
{
	for(int i=0;i<(int)(*A).size();i++)
	{
		vecAdd(&(*A)[i],&(*B)[i]);
	}
}

vector<double> vecZeros(int sizeVec)
{
	vector<double> output(sizeVec);
	return output;
}

vector<double> vecOnes(int sizeVec)
{
	vector<double> output(sizeVec);
	for(int i=0;i<sizeVec;i++){output[i]=1.;}
	return output;
}

vector< vector<double> > vecVecZeros(vector<int> sizeVec)
{
	vector< vector<double> > output(sizeVec.size());
	for(int i=0;i<(int)sizeVec.size();i++)
	{
		output[i].resize(sizeVec[i]);
		output[i]=vecZeros(sizeVec[i]);
	}
	return output;
}

vector< vector<double> > vecVecOnes(vector<int> sizeVec)
{
	vector< vector<double> > output(sizeVec.size());
	for(int i=0;i<(int)sizeVec.size();i++)
	{
		output[i].resize(sizeVec[i]);
		output[i]=vecOnes(sizeVec[i]);
	}
	return output;
}


/***mean function***/
double meanValue(std::vector<double> * vec)
{
	return sumVec(vec)/(*vec).size();
}

/***std function***/
double stdValue(std::vector<double> * vec)
{
	double mean=meanValue(vec);
	double std = 0.;
        for(int i=0; i<(int)((*vec).size()); i++)
        {
             std += ((*vec)[i] - mean) * ((*vec)[i] - mean) ;
        }
	return std;
}
/***mean function for cols of a matrix***/
std::vector<double> meanVec(std::vector<std::vector<double> > * vec)
{
	std::vector<double> meanVec((*vec).size());
	for(int i=0;i<(int)(*vec).size();i++)
	{
		meanVec[i]=meanValue(&(*vec).at(i));
	}
	return meanVec;
}
/***std function for cols of a matrix***/
std::vector<double> stdVec(std::vector<std::vector<double> > * vec)
{
	std::vector<double> stdVec((*vec).size());
	for(int i=0;i<(int)(*vec).size();i++)
	{
		stdVec[i]=stdValue(&(*vec).at(i));
	}
	return stdVec;
}

void printVecInt(std::vector<int> * vec)
{
	for(int i=0;i<(int)(*vec).size();i++)
	{
		std::cout<<(*vec)[i]<< " ";
	}
	std::cout<<std::endl;
}

void printVecFloat(std::vector<float> * vec)
{
	for(int i=0;i<(int)(*vec).size();i++)
	{
		std::cout<<(*vec)[i]<< " ";
	}
	std::cout<<std::endl;
}


void printArrayFloat(std::vector< std::vector<float> > * array)
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

void printArrayInt(std::vector< std::vector<int> > * array)
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


#endif