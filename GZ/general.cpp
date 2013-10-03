/*
 *  general.cpp
 *  GZ
 *
 *  Created by Sven Lafèbre on 20-12-07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include </user/bdegier/hisparc/lafebre2/GZ/GZ/general.h>

using namespace std;

//int sign(double x) { 
//   if (x > 0) return  1;
//   return -1;
//}

double ZenithFunc(double x) {
	return asin(2*x-1);
}

double GalacticFlux(int Z, double E) {
	E /= 1.e12; 
	return J0[Z]/1.e12 * pow(E,-gammaZ[Z]) * pow(1. + pow(E/(4.49e3*(double)Z),1.9),-1.1);
}

double TotalFlux(double E) {
	double TF = 0.; 
        for (int i = 2; i < 93; i ++) {
	TF += GalacticFlux (i, E);
        }

	return TF;
}

/*double ExtragalacticFlux(int Z, double E) {
	return 0.;
}*/

int ChargeToMass(int Z) {

	if (Z >= 0 && Z <= 92)
		return AtomicMasses[Z];

	// Return an estimate for larger values of Z
	cerr << "Warning: using estimated value A = " << int(Z*2.56) << " for Z = " << Z << endl;
	return int(Z*2.56);
}

vector<double> ProjectedDistance(const vector<double> dist, const double theta, const double phi,
								 double xi, double zeta, double latitude, double lon) {
	// Return empty vector if dimensions are wrong
	if (dist.size() != 2) {
		vector<double> A;
		return A;
	}
	
	vector<double> delta(2);
	
	double st = sin(theta);
	double ct = cos(theta);
	double sp = sin(phi);
	double cp = cos(phi);
	double sx = sin(xi);
	double cx = cos(xi);
	double sz = sin(zeta);
	double cz = cos(zeta);
	double sb = sin(latitude);
	double cb = cos(latitude);
	double sl = sin(lon);
	double cl = cos(lon);
	double se = sin(precession);
	double ce = cos(precession);
	double sxpp = sin(xi+phi);
	double cxpp = cos(xi+phi);
	double clpz = cos(lon+zeta);
	double slpz = sin(lon+zeta);

	double denominator = 4.*(sb*(ce*st + ct*se*sxpp) + 
			 cb*(-(sl*(-(cz*se*st) + ct*(cx*cp*sz + 
						ce*cz*cp*sx + ce*cz*cx*sp - sz*sx*sp))) + 
				cl*(cz*ct*cxpp + sz*(se*st - ce*ct*sxpp)
				   )));
	double numerator = 4*dist[1]*cb*cz*ct*sl*se + 
			 4*dist[1]*cb*cl*ct*se*sz - 
			 4*cxpp*(dist[0]*sb*se + 
				dist[1]*cb*clpz*st) - 
			 2*dist[0]*cb*sin(lon+zeta-xi-phi) + 
			 dist[0]*cb*sin(lon-precession+zeta-xi-phi) + 
			 dist[0]*cb*sin(lon+precession+zeta-xi-phi) - 
			 4*dist[1]*sb*se*st*sxpp + 
			 4*dist[1]*ce*(ct*sb + 
				cb*slpz*st*sxpp) + 
			 2*dist[0]*cb*sin(lon+zeta+xi+phi) + 
			 dist[0]*cb*sin(lon-precession+zeta+xi+phi) + 
			 dist[0]*cb*sin(lon+precession+zeta+xi+phi);

	delta[0] = dist[0]*(sl*(ce*cxpp*sz + cz*sxpp) + 
		   cl*(-(ce*cz*cxpp) + sz*sxpp))
		  + dist[1]*(sl*(ct*se*sz + 
			  st*(-(cz*cxpp) + ce*sz*sxpp)) - 
		   cl*(cxpp*sz*st + cz*(ct*se + ce*st*sxpp))) - 
		((cz*(-(cl*se*st) + ct*(cxpp*sl + cl*ce*sxpp)) + 
			 sz*(cl*ct*cxpp + sl*(se*st - ce*ct*sxpp)))*
		   numerator)/denominator;

	delta[1] = dist[0]*(cb*cxpp*se + 
		   sb*(ce*cxpp*slpz + clpz*sxpp)) + 
		dist[1]*(cb*(-(ce*ct) + se*st*sxpp) + 
		   sb*(sz*(cl*ct*se + cxpp*sl*st + cl*ce*st*sxpp) + 
			  cz*(ct*sl*se - cl*cxpp*st + ce*sl*st*sxpp))) - 
		((-(cb*(ce*st + ct*se*sxpp)) + 
			 sb*sl*(cz*se*st - ct*(cxpp*sz + ce*cz*sxpp)) + 
			 cl*sb*(cz*ct*cxpp + sz*(se*st - ce*ct*sxpp)))*
		   numerator)/denominator;
			
	return delta;
}

double ZenithAngle(double theta, double phi, double xi, double zeta, double b, double l) {
	double cl = cos(l+zeta);
	double sl = sin(l+zeta);
	double ce = cos(precession);
	double se = sin(precession);
	double ct = cos(theta);
	double st = sin(theta);
	double cf = cos(phi+xi);
	double sf = sin(phi+xi);
	double cb = cos(b);
	double sb = sin(b);

	return st*(sb*ce + cb*sl*se) + ct*(cf*cb*cl + sf*(sb*se - cb*sl*ce));
}

vector<double> operator+ (const vector<double> &a, const vector<double> &b) {
	int length = a.size();
	if (a.size() > b.size()) length = b.size();
	
	vector<double> r(length);
	for (int i = 0; i < length; i ++) {
		r[i] = a[i] + b[i];
	}
	return r;
}

vector<double> operator- (const vector<double> &a, const vector<double> &b) {
	int length = a.size();
	if (a.size() > b.size()) length = b.size();
	
	vector<double> r(length);
	for (int i = 0; i < length; i ++) {
		r[i] = a[i] - b[i];
	}
	return r;
}

double operator* (const vector<double> &a, const vector<double> &b) {
	int length = a.size();
	if (a.size() > b.size()) length = b.size();

	double p;
	for (int i = 0; i < length; i ++) {
		p += a[i]*b[i];
	}
	return p;
}
