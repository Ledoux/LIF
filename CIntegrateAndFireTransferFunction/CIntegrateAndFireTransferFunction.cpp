#ifndef CIntegrateAndFireTransferFunction_CPP
#define CIntegrateAndFireTransferFunction_CPP
#include "CIntegrateAndFireTransferFunction.h"
#include "../../CTool.h"
#include <map>
#include <array>
#include <complex>
#include <vector>
#include <string>
#include <stdio.h>
#include <iostream>

const double PSHRNK=-0.25;
const double PGROW=-0.2;
const double SAFETY=0.9;
const double ERRCON=1.89e-4;
const int MAXSTP=10000;
const double TINY=1e-30;
const double EPS=0.000001;
const double a1=  -1.26551223;
const double a2=  1.00002368;
const double a3=  .37409196;
const double a4=  .09678418;
const double a5=  -.18628806;
const double a6=  .27886087;
const double a7=  -1.13520398;
const double a8= 1.48851587;
const double a9= -.82215223;
const double a10= .17087277;
const double SQPI=sqrt(4.0e0*atan(1.0e0));
const double TWOPI=8.*atan(1.);
const std::string VariableNames[3]={"TotalStationaryCurrent","VoltageNoise","MembraneConstantTime"};

using namespace std;


/******************************************************/
/***************ZERO ORDER DYNAMIC*********************/


/*** compute upper and lower bounds of the integrals***/
void CIntegrateAndFireTransferFunctionClass::computeIntegralUpperBound()
{
	DoubleDict["IntegralUpperBound"]=(DoubleDict["VoltageThreshold"]-DoubleDict["TotalStationaryCurrent"])/DoubleDict["VoltageNoise"];
}
void CIntegrateAndFireTransferFunctionClass::computeIntegralLowerBound()
{
	DoubleDict["IntegralLowerBound"]=(DoubleDict["VoltageReset"]-DoubleDict["TotalStationaryCurrent"])/DoubleDict["VoltageNoise"];
}


/*** compute erf function***/
double CIntegrateAndFireTransferFunctionClass::getErrorFunction(double z)
{
	static double t,ef,at;
	static double w;
	w = fabs(z);
	t = 1.0e0/(1.0e0 + 0.5e0 * w);
	at=a1+t*(a2+t*(a3+t*(a4+t*(a5+t*(a6+t*(a7+t*(a8+t*(a9+t*a10))))))));
	ef=t*exp(at);
	if(z>0.0e0)
	ef = 2.0e0*exp(w*w)-ef;  
	return(ef);
}



/***LIF static transduction function***/
double CIntegrateAndFireTransferFunctionClass::getLIFStationaryRate()
{
	static double w,z,cont,ylow;
 	static int i,N;
	
	if(DoubleDict["VoltageNoise"]>0.)
	{
		N=10000;
		//algorithm integration
		w=0.;
		if(DoubleDict["IntegralUpperBound"]<-100.&&DoubleDict["IntegralLowerBound"]<-100.) 
		{
			w=log(DoubleDict["IntegralLowerBound"]/DoubleDict["IntegralUpperBound"])-0.25/pow(DoubleDict["IntegralUpperBound"],2.)+0.25/pow(DoubleDict["IntegralLowerBound"],2.);
			w=1./(DoubleDict["RefractoryPeriod"]+DoubleDict["MembraneConstantTime"]*w);
		}
		else if(DoubleDict["IntegralLowerBound"]<-100.) 
		{
			ylow=-100.; 
			N=(int)(100.*(DoubleDict["IntegralUpperBound"]-ylow));
			for(i=0;i<=N;i++) 
			{
				z=ylow+(DoubleDict["IntegralUpperBound"]-ylow)*(double)(i)/(double)(N);
				cont=getErrorFunction(z);
				if(i==0||i==N) {w+=0.5*cont;}
				else {w+=cont;}
			}
			w*=(DoubleDict["IntegralUpperBound"]-ylow)*SQPI/(double)(N);
			w+=log(-DoubleDict["IntegralLowerBound"]/100.)-0.000025+0.25/pow(DoubleDict["IntegralLowerBound"],2.);
			w=1./(DoubleDict["RefractoryPeriod"]+DoubleDict["MembraneConstantTime"]*w);
		}	
		else 
		{
			ylow=DoubleDict["IntegralLowerBound"];
			N=(int)(100.*(DoubleDict["IntegralUpperBound"]-ylow));
			for(i=0;i<=N;i++) 
			{
				z=ylow+(DoubleDict["IntegralUpperBound"]-ylow)*(double)(i)/(double)(N);
				cont=getErrorFunction(z);
				if(i==0||i==N){w+=0.5*cont;}
				else{w+=cont;}
			}
			w*=(DoubleDict["IntegralUpperBound"]-ylow)*SQPI/(double)(N);
			w=1./(DoubleDict["RefractoryPeriod"]+DoubleDict["MembraneConstantTime"]*w);
		}
	}
	else
	{
		//zero noise case
		if(DoubleDict["TotalStationaryCurrent"]>DoubleDict["VoltageThreshold"])
		{
			return 1./(DoubleDict["RefractoryPeriod"]+DoubleDict["MembraneConstantTime"]*log((DoubleDict["VoltageReset"]-DoubleDict["TotalStationaryCurrent"])/(DoubleDict["VoltageThreshold"]-DoubleDict["TotalStationaryCurrent"])));
		}
		else
		{
			return 0.;
		}
	}
	
	return w;
}



/******************************************************/
/***************FIRST ORDER DYNAMIC*********************/

/***LIF Linear Perturbative Transfer Function at zero frequency***/
std::map<std::string,double> CIntegrateAndFireTransferFunctionClass::getLIFPerturbativeRate0()
{
	std::map<std::string,double> OutputDict;
	static double DeltaVariable,Variable,DeltaFunction;
	static std::string VariableName;
	DeltaVariable=0.01;
	
	/**** partial derivative****/
	for(int VariableIdx=0;VariableIdx<sizeof(VariableNames)/sizeof(std::string);VariableIdx++)
	{
		/*** set the name***/
		VariableName=VariableNames[VariableIdx];
	
		/*** current***/
		Variable=DoubleDict[VariableName];
		DoubleDict[VariableName]+=DeltaVariable;
		computeIntegralUpperBound();
		computeIntegralLowerBound();
		DeltaFunction=getLIFStationaryRate();
		DoubleDict[VariableName]-=2.*DeltaVariable;
		computeIntegralUpperBound();
		computeIntegralLowerBound();
		DeltaFunction-=getLIFStationaryRate();
		DoubleDict[VariableName]=Variable;
		computeIntegralUpperBound();
		computeIntegralLowerBound();
		OutputDict[VariableName]=DeltaFunction/(2.*DeltaVariable);
	}
	/**** return OutputDict***/
	return OutputDict;
}

/***LIF Linear Perturbative Transfer Function***/

/***get the first order pass-filter***/
std::map<std::string,std::complex<double> > CIntegrateAndFireTransferFunctionClass::getRatePerturbativeRate(std::complex<double> lambda)
{
	/*******init output ********/
	std::map<std::string,std::complex<double> > OutputDict;
	OutputDict["TotalStationaryCurrent"]=1./(1.+getComplex(0.,1.)*lambda);
	OutputDict["VoltageNoise"]=1./(2.+getComplex(0.,1.)*lambda);
	OutputDict["MembraneConstantTime"]=1.;
	
	/**** return output****/
	return OutputDict;
}

/***Brunel methods****/

/***all complex function to compute RLIF***/
std::complex<double> CIntegrateAndFireTransferFunctionClass::gsurg(std::complex<double> xx,std::complex<double> yy) 
{	
	static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
	static int j;
  	static std::complex<double> xdemi,xxx,xgdemi,ydemi,ygdemi,yyy,expx,expy,ser,roro,nano,rano;
	if(xx.real()<=0.)
	{
		xdemi=xx+std::complex<double>(1.5,0.0);
		xgdemi=xx+std::complex<double>(6.5,0.0);
		xxx=xx+std::complex<double>(2.,0.0);
	}
	else if(xx.real()<=1.) 
	{
		xdemi=xx+std::complex<double>(0.5,0.0);
		xgdemi=xx+std::complex<double>(5.5,0.0);
		xxx=xx+std::complex<double>(1.,0.0);
	}
	else 
	{
		xdemi=xx+std::complex<double>(-0.5,0.0);
		xgdemi=xx+std::complex<double>(4.5,0.0);
		xxx=xx;
	}
	if(yy.real()<=0.)
	{
		ydemi=yy+std::complex<double>(1.5,0.0);
		ygdemi=yy+std::complex<double>(6.5,0.0);
		yyy=yy+std::complex<double>(2.,0.0);
	}
	else if(yy.real()<=1.)
	{
		ydemi=yy+std::complex<double>(0.5,0.0);
		ygdemi=yy+std::complex<double>(5.5,0.0);
		yyy=yy+std::complex<double>(1.0,0.0);
	}
	else 
	{
		ydemi=yy+std::complex<double>(-0.5,0.0);
		ygdemi=yy+std::complex<double>(4.5,0.0);
		yyy=yy;
	}
	expx=std::complex<double>(log(abs(xgdemi)),asin(xx.imag()/abs(xgdemi)));
	expy=std::complex<double>(log(abs(ygdemi)),asin(yy.imag()/abs(ygdemi)));
	roro=expx*xdemi-xgdemi-(expy*ydemi)+ygdemi;
	nano=std::complex<double>(exp(roro.real())*cos(roro.imag()),exp(roro.real())*sin(roro.imag()));
	ser=std::complex<double>(1.000000000190015,0.0);
	for (j=0;j<=5;j++) 
	{
		ser+=(cof[j]*(1./xxx));
		xxx+=1.;
	}
	rano=ser;
	ser=std::complex<double>(1.000000000190015,0.0);
	for (j=0;j<=5;j++) 
	{
		ser+=(cof[j]*(1./yyy));
		yyy+=1.;
	}
	rano/=ser;
	nano*=rano;
	if(xx.real()<=0.0) nano/=(xx*(1.+xx));
	else if(xx.real()<=1.0) nano/=xx;
	if(yy.real()<=0.0) nano*=(yy*(1.+yy));
	else if(yy.real()<=1.0) nano*=yy;
	return(nano);
}
	
void CIntegrateAndFireTransferFunctionClass::onefone(std::complex<double> a, std::complex<double> c, std::complex<double> z, std::complex<double> *series, std::complex<double> *deriv) 
{	
	static int n;
 	static std::complex<double> aa,cc,fac,temp;
	(*deriv)=std::complex<double>(0.0e0,0.0e0);
	fac=std::complex<double>(1.0e0,0.0e0);
	temp=fac;	
	aa=a;
	cc=c;
	for (n=1;n<=10000;n++) 
	{
		fac*=(aa/cc);
		(*deriv)=std::complex<double>(fac.real(),fac.imag());
		fac*=(1.0e0/n)*z;
		*series=temp+fac;
		if ((*series).real() == temp.real() && (*series).imag() == temp.imag()) return;
		temp= *series;
		aa+=1.;
		cc+=1.;
	}
}

void CIntegrateAndFireTransferFunctionClass::correctionlargey(std::complex<double> a, std::complex<double> c, std::complex<double> z, std::complex<double> *series) 
{
	static int n;
 	static std::complex<double> aa,cc,fac,temp;
	temp=1.;
	fac=1.;
	aa=a;
	cc=(a+1.)-c;
	for(n=1;n<=10;n++) 
	{
		fac*=aa*cc;
		fac/=((-(double)(n))*z);
		*series=temp+fac;
		if ((*series).real() == temp.real() && (*series).imag() == temp.imag()) return;
		temp=*series;
		aa+=1.;
		cc+=1.;
	}
}
	
void CIntegrateAndFireTransferFunctionClass::u01(std::complex<double> omc, double y, std::complex<double> *seriesu, int *indic)
{	
	static std::complex<double> ah,_sqom,sqom,series,deriv,asymptoty,asymptotom,ut;
  	static double alpha;
	ah=getComplex(0.,0.5)*omc;
	_sqom=sqrt(abs(omc))*exp(getComplex(0.,arg(omc)/2.));
	sqom=_sqom+getComplex(0.,1.)*_sqom;
	if(abs(omc)<10.5 && y>-3.5)
	{
	 /***** use normal series ******/
		*indic=1;
		onefone(ah,0.5,std::complex<double>(y*y,0.),&series,&deriv);
		ut=0.5*series*gsurg(ah,0.5+ah);
		onefone(0.5+ah,1.+0.5,std::complex<double>(y*y,0.),&series,&deriv);
		ut+=(y*series);
		ut*=2.*SQPI/gsurg(ah,1.);
		*seriesu=ut;
	}
	else {
		if(y<-3.5-0.15*(abs(omc)-10.5))
		{ 	
			/****** large y expansion *****/
			asymptoty=exp(-getComplex(0.,1.)*log(fabs(y))*omc);
			correctionlargey(ah,0.5,std::complex<double>(y*y,0.),&series);
			asymptoty*=series;
		}
		if(y>-3.5-0.25*(abs(omc)-10.5))
		{ 
			/******* large omega expansion *********/
			asymptotom=y*sqom;
			asymptotom+=(pow(y,3)/6.-y/2.)/sqom;
			asymptotom+=(-0.24*pow(y,2))/(sqom*sqom);
			asymptotom+=((-pow(y,5)/40.)+pow(y,3)/12.+y/8.)/(sqom*sqom*sqom);
			asymptotom+=((pow(y,4)/8.)-pow(y,2)/4.)/(sqom*sqom*sqom*sqom);
			asymptotom=SQPI*exp(0.5*y*y)*exp(asymptotom);
			asymptotom/=gsurg(0.5+ah,1.);
		}
		if(y<-3.5-0.25*(abs(omc)-10.5))
		{ 
			/****** large y expansion *****/
			*indic=2;
			*seriesu=asymptoty;
		}
		else if(y>-3.5-0.15*(abs(omc)-10.5))
		{
			 /******* large omega expansion *********/
			*indic=3;
			*seriesu=asymptotom;
		}
		else { 
			/****** interpolate between large y and large omega expansions ***/
			*indic=4;
			alpha=-((y+3.5)/(abs(omc)-10.5)+0.15)/0.1;
			*seriesu=alpha*asymptoty+(1.-alpha)*asymptotom;
		}
	}
}
	
void CIntegrateAndFireTransferFunctionClass::u23(std::complex<double> omc, double y, std::complex<double> *seriesu, int *indic)
{ 	
	static std::complex<double> ah,_sqom,sqom,series,deriv,asymptoty,asymptotom,ut,prefac;
  	static double alpha;
	ah=getComplex(0.,0.5)*omc;
	_sqom=sqrt(abs(omc))*exp(getComplex(0.,arg(omc)/2.));
	sqom=_sqom+getComplex(0.,1.)*_sqom;
	if(abs(omc)<10.5 && y>-3.5)
	{
	 	/***** use normal series ******/
		*indic=1;
		onefone(1.+ah,1.5,std::complex<double>(y*y,0.),&series,&deriv);
		ut=getComplex(0.,y)*omc*series*gsurg(ah,0.5+ah);
		onefone(0.5+ah,0.5,std::complex<double>(y*y,0.),&series,&deriv);
		ut+=series;
		ut*=(2.*SQPI)/gsurg(ah,1.);
		*seriesu=ut;
	}
	else 
	{
		if(y<-3.5-0.15*(abs(omc)-10.5))
 		{ 
			/****** large y expansion *****/
			asymptoty=exp(-getComplex(0.,1.)*log(fabs(y))*omc)/fabs(y);
			correctionlargey(ah+1.,1.5,std::complex<double>(y*y,0.),&series);
			asymptoty*=series;
			asymptoty*=std::complex<double>(0.,abs(omc));
		}
		if(y>-3.5-0.25*(abs(omc)-10.5))
		{
			 /******* large omega expansion *********/
			asymptotom=y*sqom;
			asymptotom+=(pow(y,3)/6.-y/2.)/sqom;
			asymptotom+=-0.24*pow(y,2)/(sqom*sqom);
			asymptotom+=(-pow(y,5)/40.+pow(y,3)/12.+y/8.)/(sqom*sqom*sqom);
			asymptotom+=(pow(y,4)/8.-pow(y,2)/4.)/(sqom*sqom*sqom*sqom);
			asymptotom=SQPI*exp(0.5*y*y)*exp(asymptotom);
			prefac=sqom+std::complex<double>(y,0.);
			prefac+=(0.5*y*y-0.5)/sqom;
			prefac+=-0.48*y/(sqom*sqom);
			prefac+=(-pow(y,4)/8.+pow(y,2)/4.+1./8.)/(sqom*sqom*sqom);
			prefac+=(pow(y,3)/2.-y/2.)/(sqom*sqom*sqom*sqom);
			asymptotom*=prefac/gsurg(0.5+ah,1.);
		}
		if(y<-3.5-0.25*(abs(omc)-10.5)) {
			 /****** large y expansion *****/
			*indic=2;
			*seriesu=asymptoty;
		}
		else if(y>-3.5-0.15*(abs(omc)-10.5)) {
			/******* large omega expansion *********/
			*indic=3;
			*seriesu=asymptotom;
		}
		else { 
			/****** interpolate between large y and large omega expansions ***/
			*indic=4;
			alpha=-((y+3.5)/(abs(omc)-10.5)+0.15)/0.1;
			*seriesu=(alpha*asymptoty)+(1.-alpha)*asymptotom;
		}
	}
}

/***compute the perturbative complex rate comp1.nt for LIF neuron***/ 
std::map<std::string,std::complex<double> > CIntegrateAndFireTransferFunctionClass::getBrunelLIFPerturbativeRate(std::complex<double> lambda)
{
	static std::complex<double> u01_yt,u01_yr,u23_yt,u23_yr;
	static std::complex<double> omt,series;
	static std::map<std::string,std::complex<double> > rateDict,OutputDict;
	static std::map<std::string,double > drate0Dict;
	static double rate0;
	static int indicator,indict,indich;
	omt=DoubleDict["MembraneConstantTime"]*lambda;
		
	/***** omega nul exception*****/
	if(abs(omt)==0.)
	{
		drate0Dict=getLIFPerturbativeRate0();
		OutputDict["TotalStationaryCurrent"]=drate0Dict["TotalStationaryCurrent"];
		OutputDict["VoltageNoise"]=drate0Dict["VoltageNoise"];
		DoubleDict["VoltageCapacitance"]=1.;
		OutputDict["MembraneConstantTime"]=-rate0*drate0Dict["MembraneConstantTime"]/(DoubleDict["VoltageCapacitance"]);
	}
	else
	{
		/****** U23 numerators*****/
		/***** numerator  numu *****/
		u23(omt,DoubleDict["IntegralUpperBound"],&series,&indicator);
		u23_yt=series;indict=indicator;
		u23(omt,DoubleDict["IntegralLowerBound"],&series,&indicator);
		u23_yr=series;indich=indicator;
		cout<<"BRUNEL u23 Complex "<<u23_yr.real()<<" "<<u23_yr.imag()<<" "<<indich<<endl;
		/***** numerator numu *****/
		OutputDict["TotalStationaryCurrent"]=u23_yt-u23_yr;
		OutputDict["VoltageNoise"]=DoubleDict["IntegralUpperBound"]*u23_yt-DoubleDict["IntegralLowerBound"]*u23_yr;
		OutputDict["MembraneConstantTime"]=1.;
		
		/****** U01 denominators*****/
		u01(omt,DoubleDict["IntegralUpperBound"],&series,&indicator);
		u01_yt=series;indict=indicator;
		u01(omt,DoubleDict["IntegralLowerBound"],&series,&indicator);
		u01_yr=series;indich=indicator;
		cout<<"BRUNEL u01 Complex "<<u01_yr.real()<<" "<<u01_yr.imag()<<" "<<indich<<endl;
		/***** denominator numu *****/
		OutputDict["TotalStationaryCurrent"]/=u01_yt-u01_yr;
		OutputDict["VoltageNoise"]/=u01_yt-u01_yr;
		OutputDict["VoltageNoise"]+=getComplex(0.,1.)*omt;
		OutputDict["MembraneConstantTime"]=1.;
			
		/******* rate low-pass filter*************/
		rate0=getLIFStationaryRate();
		rateDict=getRatePerturbativeRate(omt);
		cout<<"BRUNEL rate Complex "<<rateDict["TotalStationaryCurrent"].real()<<" "<<rateDict["TotalStationaryCurrent"].imag()<<" "<<indich<<endl;
		
		/******* build everything ***********/
		OutputDict["TotalStationaryCurrent"]*=rate0*rateDict["TotalStationaryCurrent"]/DoubleDict["VoltageNoise"];
		OutputDict["VoltageNoise"]*=rate0*rateDict["VoltageNoise"]/(DoubleDict["VoltageNoise"]*DoubleDict["VoltageNoise"]);
		DoubleDict["VoltageCapacitance"]=1.;
		OutputDict["MembraneConstantTime"]*=-rate0*rateDict["MembraneConstantTime"]/(DoubleDict["VoltageCapacitance"]);
		cout<<"BRUNEL rlif Complex "<<OutputDict["TotalStationaryCurrent"].real()<<" "<<OutputDict["TotalStationaryCurrent"].imag()<<" "<<indich<<endl;
		
	
	
	}
	
	//cout<<"BRUNEL rate Complex "<<mu_rate.real()<<" "<<mu_rate.imag()<<endl;
	/***** return output****/
	return OutputDict;
}

/****Hakim and Ostojic Methods***/
void CIntegrateAndFireTransferFunctionClass::rkck(std::vector<std::complex<double> > y, std::vector<std::complex<double> > dydx, double x,  double h, std::vector<std::complex<double> >& yout, std::vector<std::complex<double> >& yerr, std::complex<double> lamda, void (CIntegrateAndFireTransferFunctionClass::*derivs)(double, vector<std::complex<double> >&, std::vector<std::complex<double> >&, std::complex<double> )) {

  int i;
  double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
    b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
    b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
    b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
    b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
    c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
    dc5 = -277.00/14336.0;
  double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
    dc4=c4-13525.0/55296.0,dc6=c6-0.25;
  
  int n=y.size();
  std::vector<std::complex<double> >  ak2(n), ak3(n), ak4(n), ak5(n), ak6(n);
  std::vector<std::complex<double> >  ytemp(n);
  
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+b21*h*dydx[i];
  (*this.*derivs)(x+a2*h,ytemp,ak2,lamda);
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
  (*this.*derivs)(x+a3*h,ytemp,ak3,lamda);
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
  (*this.*derivs)(x+a4*h,ytemp,ak4,lamda);
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
  (*this.*derivs)(x+a5*h,ytemp,ak5,lamda);
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
  (*this.*derivs)(x+a6*h,ytemp,ak6,lamda);
  for (i=0;i<n;i++){
    yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
  }
  for (i=0;i<n;i++)
    yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
}

void CIntegrateAndFireTransferFunctionClass::rkqs(std::vector<std::complex<double> >& y, std::vector<std::complex<double> > dydx,  double& x, const double htry, const double eps, const std::vector<double>& yscal, double& hdid, double& hnext, std::complex<double> lamda, void (CIntegrateAndFireTransferFunctionClass::*derivs)(double, std::vector<std::complex<double> >&, std::vector<std::complex<double> >&, std::complex<double> )) {
  int n=y.size();
  int i;
  double errmax,h,htemp,xnew;
  std::vector<std::complex<double> > yerr(n), ytemp(n);
  h=htry;
  for (;;) {
    rkck(y,dydx,x,h,ytemp,yerr,lamda,derivs);
    errmax=0.0;
    for (i=0;i<n;i++) errmax = max(errmax, abs(yerr[i]/yscal[i]));
    errmax /= eps;
    if (errmax <= 1.0) break;
    htemp=SAFETY*h*pow(errmax,PSHRNK);
    h=(h >= 0.0 ? max(htemp,0.1*h) : min(htemp,0.1*h));
    xnew=(x)+h;
    if (xnew == x){
      cout<<"stepsize underflow in rkqs; if running crossoverfield quintessence, J may be too low\n";
      exit(0);
    }
  }  
  if (errmax > ERRCON) hnext=SAFETY*h*pow(errmax,PGROW);
  else hnext=5.0*h;
  x += (hdid=h);
  for (i=0;i<n;i++) y[i]=ytemp[i];  
}

double CIntegrateAndFireTransferFunctionClass::inside_integ(double x)
{
  float w;
  float y,ymin=-20.;
  int i,N=10000;
  w=0.0e0;
  for(i=0;i<=N;i++) {       
    y=ymin+(x-ymin)*(double)(i)/(double)(N);
    if(i==0||i==N) w+=0.5*pow(getErrorFunction(y)*exp(0.5*(x*x-y*y)),2.);
    else w+=pow(getErrorFunction(y)*exp(0.5*(x*x-y*y)),2.);
  }
  w*=(x-ymin)/(double)(N);
  return(w);
}


void CIntegrateAndFireTransferFunctionClass::phi_derivs(double y, std::vector<std::complex<double> >& phi, std::vector<std::complex<double> >& dphi, std::complex<double> lamda){

  dphi[0]=phi[1];
  dphi[1]=2.*((lamda-1.)*phi[0]-y*phi[1]);
}

void CIntegrateAndFireTransferFunctionClass::phi_t_derivs(double y, std::vector<std::complex<double> >& phi, std::vector<std::complex<double> >& dphi, std::complex<double> lamda){

  dphi[0]=phi[1];
  dphi[1]=2.*(lamda*phi[0]+y*phi[1]);
}

double CIntegrateAndFireTransferFunctionClass::odeint(vector<std::complex<double> >& ystart, const double x1, const double x2, const double eps, const double h1,  const double hmin, std::complex<double> lamda, void (CIntegrateAndFireTransferFunctionClass::*derivs)(double, std::vector<std::complex<double> >&, std::vector<std::complex<double> >&, std::complex<double>)){
  int nstp,i;
  double x,hdid,h,hnext;
  double sum=0;  
  int nvar=ystart.size();
  std::vector<double> yscal(nvar);
  std::vector<std::complex<double> > y(nvar), dydx(nvar);
  x=x1;
  (h1*(x2-x1)<0)?h=-1:h=1;
  for (i=0;i<nvar;i++) y[i]=ystart[i];
  for (nstp=1;nstp<=MAXSTP;nstp++) {
 
    (*this.*derivs)(x,y,dydx,lamda);
    for (i=0;i<nvar;i++) yscal[i]=abs(y[i])+abs(dydx[i]*h)+TINY;
    
    if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
    
    rkqs(y,dydx,x,h,eps, yscal ,hdid,hnext,lamda,derivs);
    sum=sum+hdid*real(y[0]);

    if ((x-x2)*(x2-x1) >= 0.0) {
      for (i=0;i<nvar;i++) ystart[i]=y[i];
      return sum;
    }
    if (fabs(hnext) <= hmin){
      cout<<"Step size too small in odeint\n";
    }
    h=hnext;
  }
  return 0.;
}


std::complex<double> CIntegrateAndFireTransferFunctionClass::phi_t_rka(double y, std::complex<double> lambda, std::complex<double>& phi_pr){
  int nb_steps=10000; 
  double eps=0.000001;
  double h1=0.2;
  double hmin=0.000000001;
 
  std::complex<double> out;
  double y_start=DoubleDict["IntegralLowerBound"]-15.;
  if(y<y_start){
    cout<<"phi2: y too small\n";
  }
  std::vector< std::complex<double> > phi_start(2);
  phi_start[0]=std::complex<double>(1.,0.);
  phi_start[1]=std::complex<double>(2.*y_start,0.);
  deriv=&CIntegrateAndFireTransferFunctionClass::phi_t_derivs;
  double nstp=odeint(phi_start,y_start,y,eps,h1,hmin,lambda,deriv);
  out=phi_start[0];
  phi_pr=phi_start[1];

  return out;
}

std::map<std::string,std::complex<double> > CIntegrateAndFireTransferFunctionClass::getHakimLIFPerturbativeRate(std::complex<double> lambda){
	
	static std::map<std::string,std::complex<double> > OutputDict,rateDict;
	static std::complex<double> omt;
	static double rate0;
	omt=getComplex(0.,1.)*DoubleDict["MembraneConstantTime"]*lambda;
	
	/****** hypergeo compute*****/
	std::complex<double> phi_yth_pr;
	std::complex<double> phi_yth=phi_t_rka(DoubleDict["IntegralUpperBound"],omt,phi_yth_pr);
	std::complex<double> phi_yr_pr;
	std::complex<double> phi_yr=phi_t_rka(DoubleDict["IntegralLowerBound"],omt,phi_yr_pr);
	cout<<"Class "<<DoubleDict["IntegralUpperBound"]<<" "<<DoubleDict["IntegralLowerBound"]<<" "<<phi_yth.real()<<" "<<phi_yth.imag()<<endl;
	
	OutputDict["TotalStationaryCurrent"]=(phi_yr_pr-phi_yth_pr)/((1.+omt)*(phi_yr-phi_yth));
	OutputDict["VoltageNoise"]=((DoubleDict["IntegralLowerBound"]*phi_yr_pr-DoubleDict["IntegralUpperBound"]*phi_yth_pr)/(phi_yr-phi_yth)-2.)/(2.+omt)+1.;
	OutputDict["MembraneConstantTime"]=1.;
	
	/******* rate low-pass filter*************/
	rate0=getLIFStationaryRate();
	OutputDict["TotalStationaryCurrent"]*=rate0/(DoubleDict["VoltageNoise"]);
	OutputDict["VoltageNoise"]*=rate0/(DoubleDict["VoltageNoise"]*DoubleDict["VoltageNoise"]);
	DoubleDict["VoltageCapacitance"]=1.;
	OutputDict["MembraneConstantTime"]*=-rate0*rateDict["MembraneConstantTime"]/(DoubleDict["VoltageCapacitance"]);

	/******* return output************/
	return OutputDict;

}




#endif




