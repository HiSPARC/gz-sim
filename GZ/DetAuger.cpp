/*
 *  DetAuger.cpp
 *  GZ
 *
 *  Created by Sven Lafèbre on 21-12-07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "DetAuger.h"

DetAuger::DetAuger() : Detector(0) {
	_latitude  = -0.6150;
	_longitude = -1.2092;

	_minDistance[0] = -3.3e4;
	_minDistance[1] = -3.3e4;
	_maxDistance[0] =  3.3e4;
	_maxDistance[1] =  3.3e4;
	
	_effectiveAreaIterations = 10;
	_detectionProbabilityIterations = 10;
}

DetAuger::DetAuger(int s) : Detector(s) {
	_latitude  = -0.6150;
	_longitude = -1.2092;

	_minDistance[0] = -3.3e4;
	_minDistance[1] = -3.3e4;
	_maxDistance[0] =  3.3e4;
	_maxDistance[1] =  3.3e4;
	
	_effectiveAreaIterations = 10;
	_detectionProbabilityIterations = 10;
}

bool DetAuger::IsDetected(const std::vector<double> &r) const {
	// No detection for close events closer together than 5 km
	if (r*r < 2.5e7) return false;
	
	if (r[1] < 3.2e4 &&
		(r[0] - 2.6e4)*(r[0] - 2.6e4) + (r[1] - 3.3e4)*(r[1] - 3.3e4) > 2e4 &&
		sqrt(3.)*r[0] + r[1] < 7.7e4 &&
		r[0] < 3.28e4 &&
		sqrt(3.)*r[0] - r[1] < 4.6e4 &&
		r[0] - 1.8*r[1] < 5.9e4 &&
		r[1] > -3.24e4 &&
		1.8*r[1] + r[0] > -6.6e4 &&
		sqrt(3.)*r[0] + r[1] > -5.6e4 &&
		(r[0] + 2.5e4)*(r[0] + 2.5e4) + (r[1] + 1.4e4)*(r[1] + 1.4e4) > 8e3 &&
		r[0] > -3.28e4 &&
		sqrt(3.)*r[0] - r[1] > -6e4 &&
		r[0] - 2.8*r[1] > -8.7e4) return true;

	return false;
}

double DetAuger::TotalArea() const {
	return 3e9;
}

double DetAuger::AngularSensitivity(double zenith) const {
	// Zenith angle must be smaller than 60 degrees
	return zenith < pi/3.;
}

double DetAuger::EnergyEfficiency(double energy) const {
	return erf((log10(energy)-18.)*1.7)/2.+.5;
//	return energy >= 1e18;
}
