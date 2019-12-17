
#include "TMath.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"

#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/AllIntegrationTypes.h"
#include "Math/Functor.h"
#include "Math/GaussIntegrator.h"

#include <iostream>
using namespace std;

double Z = 79.;
double pi = 3.1415926;
double hbarc = 0.197327053;
double A = 197.0;


// why is there no way to pass parameters?
double _q2 = 0;
double _r = 0;
ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,0,0.001,100);

ROOT::Math::Integrator igrho(ROOT::Math::IntegrationOneDim::kADAPTIVE,0,0.001,100);

double rho( double r, double R=6.38, double a=0.535, double w = 0.0 ) {
	// R is the radius
	// a is the skin depth
	// w is the shape parameter

	// cout << "r = " << r << endl;
	// cout << "R = " << R << ", a = " << a << ", w = " << w << endl;

	double b = exp( (r - R) / a );
	double c = 1 + (w*w * r*r / R*R);
	double rho0 = (3 * Z) / ( 4 * pi * pow(R, 3) );
	return (rho0 * c) / ( 1 + b ) *0.9362;
	// this 0.936 is to normalize the FF to be = 1 at q^2 = 0
}

double f1rho( double *x, double*par ){
	return rho( x[0], par[0], par[1], par[2] );
}


// perform 1D integral over r at a given q
double FormFactorIntegral(double x ){
	double r = x;
	double q = sqrt(_q2);

	double rv = r * rho(r) * sin( q * r );
	return rv;
}

double FormFactor( double *x, double *par ){
	double q2 = x[0];
	_q2 = q2 / ( hbarc * hbarc );
	double q = sqrt(_q2);


	double a0 = 4 * pi / (Z * q);
	
	double v = ig.Integral( 0, 20 );
	// cout << "v = " << v << endl;
	return a0 * v;
}

double FormFactorQ( double *x, double *par ){
	double xp = x[0] * x[0];
	double parp = 0;
	return pow(FormFactor( &xp, &parp ), 2);
}

double ZhaFormFactor( double *x, double *par ) {

	double A23 = pow(A, -2. / 3.);
	double R = pow(A, (1./3.))*1.16 * (1. - 1.16 * A23);
	
	double q = sqrt(x[0]);
	// q = q / 5.;
	double q2 = q*q;
	double a1 = q * R / hbarc;
	double a2 = hbarc / (q * 1.16 * (1. - 1.16 * A23 ));
	double sph  = (sin(a1) - a1  * cos(a1 )) * 3. * a2 * a2 * a2 / A;
	double a0 = 0.70;

	return sph / (1. + (a0 * a0 * q2) / (hbarc * hbarc));
}

double ZhaFormFactorQ( double *x, double *par ) {

	double A23 = pow(A, -2. / 3.);
	double R = pow(A, (1./3.))*1.16 * (1. - 1.16 * A23);
	
	double q = x[0];
	// q = q / 5.;
	double q2 = q*q;
	double a1 = q * R / hbarc;
	double a2 = hbarc / (q * 1.16 * (1. - 1.16 * A23 ));
	double sph  = (sin(a1) - a1  * cos(a1 )) * 3. * a2 * a2 * a2 / A;
	double a0 = 0.70;

	return sph / (1. + (a0 * a0 * q2) / (hbarc * hbarc));
}


double invFFIntegral( double x ){
	double q = x;
	// double px[] = { q*q / (5.07*5.07) };
	// double vFF = ZhaFormFactor( px, 0 );
	double px[] = { q*q/ (5.07*5.07) };
	double vFF = FormFactor( px, 0 );
	double a0 = sin( q * _r ) * (q / _r);
	return vFF * a0;

}

double invFF( double *x, double *par ) {
	double r = x[0];
	_r = r;

	double v = igrho.Integral( 0, 10 );
	return v / ( 2 * pi * pi ) *1e2 * 0.8;
}


double photonkt( double *x, double *par ){
	double alpha = 1.0 / 137.0;
	
	double kt = x[0];
	double k = 0.775 / 2.0; // y = ln( 2k / M_v ) -> mid rapidity, k = M_v / 2;
	double gamma = 108;
	double q2 = kt*kt + k*k / (gamma*gamma);

	// prefactors are to normalize maximum to unity
	return 1e-3 * 1.81 * pow( kt * alpha * Z * FormFactor( &q2, par ), 2) / pow( pi * q2, 2 );
}


int main() {
	TCanvas * can = new TCanvas( "can", "", 900, 600 );
	can->Print( "rp.pdf[" );

	ROOT::Math::Functor1D FFInt(&FormFactorIntegral);
	// ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kVEGAS,0,0.001,50000);
	ig.SetRelTolerance(0.001);
	ig.SetAbsTolerance(0.0);
	ig.SetFunction(FFInt);

	ROOT::Math::Functor1D invFFInt(&invFFIntegral);
	igrho.SetRelTolerance(0.001);
	igrho.SetAbsTolerance(0.0);
	igrho.SetFunction(invFFInt);

	TH1 * hFrame = new TH1F( "hFrame", "hFrame", 500, 0, 25 );
	hFrame->GetYaxis()->SetRangeUser( 0, 0.21 );
	hFrame->Draw();

	TF1 * fRho = new TF1( "fRho", &f1rho, 0, 25, 3 );
	fRho->SetParameters( 6.38, 0.535, 0.0 );
	fRho->Draw("same");
	cout << "Irho = " << fRho->Integral(0, 100) << endl;
	can->Print( "rp.pdf" );
	cout << fRho->Eval( 1 ) << endl;

	TF1 * fPhotonkt = new TF1( "fPhotonkt", &photonkt, 1e-8, 1.0, 1 );
	fPhotonkt->SetRange( 1e-8, 0.25 );
	fPhotonkt->SetLineColor(kOrange);
	fPhotonkt->SetNpx( 1000 );

	TF1 * fFFQ = new TF1( "fFFQ", &FormFactorQ, 1e-8, 1.0, 1 );
	fFFQ->SetRange( 1e-8, 0.25 );
	fFFQ->SetNpx( 2000 );

	TF1 * fFF = new TF1( "fFF", &FormFactor, 1e-8, 1.0, 1 );
	fFF->SetRange( 1e-8, 5 );
	fFF->SetNpx( 1000 );
	
	fFF->SetLineStyle( 7 );
	

	TF1 * fZhaFFQ = new TF1( "fZhaFFQ", &ZhaFormFactorQ, 1e-8, 1.0, 1 );
	fZhaFFQ->SetRange( 1e-8, 5 );
	fZhaFFQ->SetLineColor(kOrange);
	fZhaFFQ->SetNpx( 5000 );


	TF1 * fZhaFF = new TF1( "fZhaFF", &ZhaFormFactor, 1e-8, 1.0, 1 );
	fZhaFF->SetRange( 1e-8, 5 );
	fZhaFF->SetLineColor(kOrange);
	fZhaFF->SetNpx( 5000 );
	
	// hFrame->GetYaxis()->SetRangeUser( -0.1, 1.2 );
	// hFrame->GetXaxis()->SetRangeUser( 1e-8, 1.0 );
	// hFrame->Draw();
	
	fZhaFF->Draw("");
	fFF->Draw("same");
	
	gPad->SetLogx(1);


	can->Print( "rp.pdf" );
	
	cout << fFF->Eval( 1e-5 ) << endl;


	TF1 * fIFF = new TF1( "fIFF", &invFF, 0, 25, 1 );
	fIFF->SetNpx(500);
	fIFF->Draw();
	fIFF->SetLineColor(kBlue);
	fIFF->SetLineStyle(7);
	gPad->SetLogx(0);
	fRho->Draw("same");
	can->Print( "rp.pdf" );


	TFile * fOut = new TFile( "SnapshotFormFactor.root", "RECREATE" );
	fPhotonkt->Write();
	fFF->Write();
	fFFQ->Write();
	fZhaFF->Write();
	fZhaFFQ->Write();
	fRho->Write();
	fIFF->Write();
	fOut->Write();

	can->Print( "rp.pdf]" );

}