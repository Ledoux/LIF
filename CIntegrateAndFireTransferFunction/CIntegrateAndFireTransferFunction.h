#ifndef CIntegrateAndFireTransferFunction_H
#define  CIntegrateAndFireTransferFunction_H
#include "../../CTool.h"
#include <map>
#include <iostream>
#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <math.h>
#include <string>
#include <vector>

class CIntegrateAndFireTransferFunctionClass : public CToolClass
{
	public:
	
		/*** COMPUTE METHODS ***/
		
		/*** INTEGRALS BOUNDARIES***/
		void computeIntegralUpperBound();
		void computeIntegralLowerBound();
		
		/*** GET METHODS ***/
		
		/*** STATIONARY RATES ****/
		double getErrorFunction(double z);
		double getLIFStationaryRate();
		
		/*** LINEAR PERTURBATIONS ***/
		
		/*** RatePerturbation***/
		std::map<std::string,std::complex<double> > getRatePerturbativeRate(std::complex<double> lambda);
		
		/*** LIF perturbation at zero frequency***/
		std::map<std::string,double> getLIFPerturbativeRate0();
		
		/*** Brunel Methods ***/
		std::complex<double> gsurg(std::complex<double> xx,std::complex<double> yy);
		void onefone(std::complex<double> a, std::complex<double> c, std::complex<double> z, std::complex<double> *series, std::complex<double> *deriv);
		void correctionlargey(std::complex<double> a, std::complex<double> c, std::complex<double> z, std::complex<double> *series); 
		void u01(std::complex<double> omc, double y, std::complex<double> *seriesu, int *indic);
		void u23(std::complex<double> omc, double y, std::complex<double> *seriesu, int *indic);
		std::map<std::string,std::complex<double> > getBrunelLIFPerturbativeRate(std::complex<double> lambda);
		
		/*** Hakim Ostojic Methods ***/
		double inside_integ(double x);
		void (CIntegrateAndFireTransferFunctionClass::*deriv)(double, std::vector<std::complex<double> >&, std::vector<std::complex<double> >&, std::complex<double>);//pointer to functions like derivs
		void phi_derivs(double y, std::vector<std::complex<double> >& phi, std::vector<std::complex<double> >& dphi, std::complex<double> lamda);
		void phi_t_derivs(double y, std::vector<std::complex<double> >& phi, std::vector<std::complex<double> >& dphi, std::complex<double> lamda);
		void rkck(std::vector<std::complex<double> > y, std::vector<std::complex<double> > dydx, double x,  double h, std::vector<std::complex<double> >& yout, std::vector<std::complex<double> >& yerr, std::complex<double> lamda, void (CIntegrateAndFireTransferFunctionClass::*derivs)(double, std::vector<std::complex<double> >&, std::vector<std::complex<double> >&, std::complex<double> ));
		void rkqs(std::vector<std::complex<double> >& y, std::vector<std::complex<double> > dydx,  double& x, const double htry, const double eps, const std::vector<double>& yscal, double& hdid, double& hnext, std::complex<double> lamda, void (CIntegrateAndFireTransferFunctionClass::* derivs)(double, std::vector<std::complex<double> >&, std::vector<std::complex<double> >&, std::complex<double> ));
		double odeint(std::vector<std::complex<double> >& ystart, const double x1, const double x2, const double eps, const double h1,  const double hmin, std::complex<double> lamda,
						void (CIntegrateAndFireTransferFunctionClass::*derivs)(double, std::vector<std::complex<double> >&, std::vector<std::complex<double> >&, std::complex<double>));
		std::complex<double> phi_t_rka(double y, std::complex<double> lamda, std::complex<double>& phi_pr);
		std::map<std::string,std::complex<double> > getHakimLIFPerturbativeRate(std::complex<double> lamda);
		

};

//class CToolClass;
//class CIntegrateAndFireTransferFunctionClass;

//class CIntegrateAndFireTransferFunctionClass: public CToolClass {
//	public:
//
//}


#endif



