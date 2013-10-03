/*
 *  DetLofar.cpp
 *  GZ
 *
 *  Created by Sven Lafèbre on 21-12-07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "DetLofar.h"

DetLofar::DetLofar() : Detector(0) {
	_latitude  = 0.9234;
	_longitude = 0.1199;

	_remoteSize = 200.;
	_coreSize   = 5.e3;

	_minDistance[0] = -1.e4;
	_minDistance[1] = -1.e4;
	_maxDistance[0] =  1.e4;
	_maxDistance[1] =  1.e4;
	
	_effectiveAreaIterations = 10;
	_detectionProbabilityIterations = 10;
}

DetLofar::DetLofar(int s) : Detector(s) {
	_latitude  = 0.9234;
	_longitude = 0.1199;

	_remoteSize = 200.;
	_coreSize   = 5.e3;

	_minDistance[0] = -1.e4;
	_minDistance[1] = -1.e4;
	_maxDistance[0] =  1.e4;
	_maxDistance[1] =  1.e4;
	
	_effectiveAreaIterations = 10;
	_detectionProbabilityIterations = 10;
}

bool DetLofar::IsDetected(const std::vector<double> &r) const {
	// No detection for close events closer together than 1 km
	if (r*r < 1e6) return false;
	
	vector<double> pos(2, 0.);

	if (InCircle(r, pos, _coreSize))
		return true;

	double rho = sqrt(r*r);
	for (double i = 0.75; i < 3.5; i = 0.1+i*1.12) {
		double dist = 1.e6*exp(-1.5*pi+pi*i/3.);

		// These two lines make the algortihm much more efficient
		if (dist > rho + _remoteSize) return false;
		if (fabs(dist-rho) > _remoteSize) continue;

/*		for (int j = 0; j < 5; j ++) {
			pos[0] = dist*cos(pi*(i/3.+j*2/5.));
			pos[1] = dist*sin(pi*(i/3.+j*2/5.));
			
			if (InCircle(r, pos, _remoteSize))
				return true;
		}*/
	}

	return false;
}

double DetLofar::TotalArea() const {
	// ~ 8.41947e7 m^2
	return 45*pi*_remoteSize*_remoteSize + pi*_coreSize*_coreSize;
}

double DetLofar::AngularSensitivity(double zenith) const {
	// Zenith angle must be smaller than 80 degrees
	return zenith < 4.*pi/9.;
}

double DetLofar::EnergyEfficiency(double energy) const {
	return erf((log10(energy)-17.)*1.7)/2.+.5;
//	return energy > 1.2e16;
}

bool DetLofar::InCircle(const vector<double> &pos,
						const vector<double> &center,
						double radius) const {
	return (pos-center)*(pos-center) <= radius*radius;
}
