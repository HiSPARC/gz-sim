/*
 *  DetDense.cpp
 *  GZ
 *
 *  Created by Sven Lafèbre on 21-12-07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "DetDense.h"

DetDense::DetDense() : Detector(0) {
	_latitude  =  pi/4.;
	_longitude =  0.;
	_radius = 1e-3;

	_minDistance[0] = -_radius;
	_minDistance[1] = -_radius;
	_maxDistance[0] =  _radius;
	_maxDistance[1] =  _radius;
	
	_effectiveAreaIterations = 1;
	_detectionProbabilityIterations = 100;
}

DetDense::DetDense(int s) : Detector(s) {
	_latitude  =  pi/4.;
	_longitude =  0.;
	_radius = 110e3;

	_minDistance[0] = -_radius;
	_minDistance[1] = -_radius;
	_maxDistance[0] =  _radius;
	_maxDistance[1] =  _radius;
	
	_effectiveAreaIterations = 10;
	_detectionProbabilityIterations = 10;
}

bool DetDense::IsDetected(const std::vector<double> &r) const {
	// No detection for close events closer together than 1 km
	if (r*r < 1e6) return false;
	
	if (r*r < _radius*_radius) return true;

	return false;
}

double DetDense::TotalArea() const {
	return pi*_radius*_radius;
}

double DetDense::AngularSensitivity(double zenith) const {
	// Zenith angle must be smaller than 60 degrees
	return zenith < pi/3.;
}

double DetDense::EnergyEfficiency(double energy) const {
	return erf((log10(energy)-17.5)*1.7)/2.+.5;

//	return energy > 1e17;
}
