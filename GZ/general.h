/*
 *  general.h
 *  GZ
 *
 *  Created by Sven Lafèbre on 26-11-07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __GENERAL_H__
#define __GENERAL_H__ 1

#include <iostream>
#include <vector>
#include <math.h>

# define pi         3.141592653589793
# define Mdip       1.007050342e-11//T*AU^3
# define Mspot      5.0352522e-12//T*AU^3
# define rho0       1.495979e11//m
# define rho02      2.237953e22//m^2
# define Brho0    (-3.5e-9)//T
# define Bphi0      3.5e-9//T
# define Rsun       6.96e8//m
# define Rearth     6.367e6//m
# define c          2.997924e8//m/s
# define eV         1.602176e-19//C
# define kBT        0.5//eV
# define MaxDistance (4*rho0)
# define MinEnergy  1.e16//eV
# define MaxEnergy  1.e20//eV

# define precession pi/2.    // precession angle of Earth
# define cose       0.0      // cosine of precession
# define sine       1.0      // sine of precession

using namespace std;

const int AtomicMasses[] = {
	  1,												//  0 (neutron)
	  1,   4,   7,   9,  11,  12,  14,  16,  19,  20,   //  1 - 10 
	 23,  24,  27,  28,  31,  32,  35,  40,  39,  40,	// 11 - 20
	 45,  48,  51,  52,  55,  56,  59,  59,  64,  65,	// 21 - 30
	 70,  73,  75,  79,  80,  84,  85,  88,  89,  91,	// 31 - 40
	 93,  96,  98, 101, 103, 106, 108, 112, 115, 119,	// 41 - 50	
	122, 128, 127, 131, 133, 137, 139, 140, 141, 144,	// 51 - 60
	145, 150, 152, 157, 159, 163, 165, 167, 169, 173,	// 61 - 70
	175, 178, 181, 184, 186, 190, 192, 195, 197, 201,	// 71 - 80
	204, 207, 209, 209, 210, 222, 223, 226, 227, 232,	// 81 - 90
	231, 238											// 91 - 92
};

/// Absolute flux J_0(Z) per m^2 sr s TeV
const double J0[] = { 0.,			// Padding for J0[0]
	8.73e-2, 5.71e-2, 2.08e-3, 4.74e-4, 8.95e-4, 1.06e-2, 2.35e-3, 1.57e-2, 3.28e-4, 4.60e-3,
	7.54e-4, 8.01e-3, 1.15e-3, 7.96e-3, 2.70e-4, 2.29e-3, 2.94e-4, 8.36e-4, 5.36e-4, 1.47e-3,
	3.04e-4, 1.14e-3, 6.31e-4, 1.36e-3, 1.35e-3, 2.04e-2, 7.51e-5, 9.96e-4, 2.18e-5, 1.66e-5,
	2.75e-6, 4.02e-6, 9.99e-7, 2.11e-6, 1.34e-6, 1.30e-6, 6.93e-7, 2.11e-6, 7.82e-7, 8.42e-7,
	5.05e-7, 7.79e-7, 6.98e-8, 3.01e-7, 3.77e-7, 5.10e-7, 4.54e-7, 6.30e-7, 1.61e-7, 7.15e-7,
	2.03e-7, 9.10e-7, 1.34e-7, 5.74e-7, 2.79e-7, 1.23e-6, 1.23e-7, 5.10e-7, 9.52e-8, 4.05e-7,
	8.30e-8, 3.68e-7, 1.58e-7, 6.99e-7, 1.48e-7, 6.27e-7, 8.36e-8, 3.52e-7, 1.02e-7, 4.15e-7,
	1.72e-7, 3.57e-7, 2.16e-7, 4.16e-7, 3.35e-7, 6.42e-7, 6.63e-7, 1.03e-6, 7.70e-7, 7.43e-7,
	4.28e-7, 8.06e-7, 3.25e-7, 3.99e-7, 4.08e-8, 1.74e-7, 1.78e-8, 7.54e-8, 1.97e-8, 8.87e-8,
	1.71e-8, 3.54e-7,
};

const double gammaZ[] = { 0.,		// Padding for gammaZ[0]
	2.71, 2.64, 2.54, 2.75, 2.95, 2.66, 2.72, 2.68, 2.69, 2.64,
	2.66, 2.64, 2.66, 2.75, 2.69, 2.55, 2.68, 2.64, 2.65, 2.70,
	2.64, 2.61, 2.63, 2.67, 2.46, 2.59, 2.72, 2.51, 2.57, 2.56,
	2.55, 2.54, 2.54, 2.53, 2.52, 2.51, 2.51, 2.50, 2.49, 2.48,
	2.47, 2.46, 2.46, 2.45, 2.44, 2.43, 2.42, 2.41, 2.40, 2.39,
	2.38, 2.37, 2.37, 2.36, 2.35, 2.34, 2.33, 2.32, 2.31, 2.30,
	2.29, 2.28, 2.27, 2.25, 2.24, 2.23, 2.22, 2.21, 2.20, 2.19,
	2.18, 2.17, 2.16, 2.15, 2.13, 2.12, 2.11, 2.10, 2.09, 2.08,
	2.06, 2.05, 2.04, 2.03, 2.02, 2.00, 1.99, 1.98, 1.97, 1.96,
	1.94, 1.93,
};

//int sign(double x);

double ZenithFunc(double x);

/// \brief Calculate the flux density of a cosmic particle
///
/// Calculate the galactic flux of cosmic ray particles given a certain particle species and energy, given
/// the parametrisation (J.Hörandel, Aph 19 (2003) 193):
/// \f[
///    \Phi_0(Z,E) = J_0(Z)
///                    \left(\frac{E}{1\,\hbox{TeV}}\right)^-\gamma(Z)
///                    \left[1 + \left(\frac{E}{4490Z}\right)^{1.9}\right]^{-1.1}
/// \f]
/// \param	Z	Particle charge number Z
/// \param	E	Particle energy (in eV)
/// \return	Flux in particles per square meter per second per electronvolt
double GalacticFlux(int Z, double E);

double TotalFlux(double E);


/*double ExtragalacticFlux(int Z, double E);*/

int ChargeToMass(int Z);

vector<double> ProjectedDistance(const vector<double> dist, const double theta, const double phi,
								 double xi, double zeta, double latitude, double lon);

double ZenithAngle(double theta, double phi, double xi, double zeta, double b, double l);

vector<double> operator+(const vector<double> &a, const vector<double> &b);

vector<double> operator-(const vector<double> &a, const vector<double> &b);

double         operator*(const vector<double> &a, const vector<double> &b);

#endif // __GENERAL_H__
