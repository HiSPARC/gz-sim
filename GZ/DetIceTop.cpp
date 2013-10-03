/*
 *  DetIceTop.cpp
 *  GZ
 *
 *  Created by Sven Lafebre on 5-05-09.
 *  Copyright 2009 Pennsylvania State University. All rights reserved.
 *
 */

#include "DetIceTop.h"

DetIceTop::DetIceTop() : Detector(0) {
	// Location of the detector on Earth
	_latitude  = -0.6150;
	_longitude = -1.2092;

	// Maximum extent of the detector (m)
	_minDistance[0] = -6.e2;
	_minDistance[1] = -6.e2;
	_maxDistance[0] =  6.e2;
	_maxDistance[1] =  6.e2;
	
	_coreSize = 6.e2;
	
	// Number of iterations to run
	_effectiveAreaIterations = 10;
	_detectionProbabilityIterations = 10;
}

DetIceTop::DetIceTop(int s) : Detector(s) {
	// Location of the detector on Earth
	_latitude  = -0.6150;
	_longitude = -1.2092;
	
	// Maximum extent of the detector (m)
	_minDistance[0] = -6.e2;
	_minDistance[1] = -6.e2;
	_maxDistance[0] =  6.e2;
	_maxDistance[1] =  6.e2;
	
	_coreSize = 6.e2;
	
	// Number of iterations to run
	_effectiveAreaIterations = 10;
	_detectionProbabilityIterations = 10;
}

bool DetIceTop::IsDetected(const std::vector<double> &r) const {
	// No detection for events closer together than 125 m
	if (r*r < 15625) return false;

	vector<double> pos(2, 0.);
	
	if (InCircle(r, pos, _coreSize))
		return true;
	
	return false;
}

double DetIceTop::TotalArea() const {
	return 1.13e6;
}

double DetIceTop::AngularSensitivity(double zenith) const {
	// Zenith angle must be smaller than 45 degrees
	return erf((zenith-pi/4.)*10.)/2.+.5;
}

double DetIceTop::EnergyEfficiency(double energy) const {
	return erf((log10(energy)-18.)*1.7)/2.+.5;
//	return energy >= 1e18;
}

bool DetIceTop::InCircle(const vector<double> &pos,
						const vector<double> &center,
						double radius) const {
	return (pos-center)*(pos-center) <= radius*radius;
}
