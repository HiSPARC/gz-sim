/*
 *  Detector.cpp
 *  GZ
 *
 *  Created by Sven Lafèbre on 17-12-07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "Detector.h"

using namespace std;

Detector::Detector() :
		_randomGenerator(0),
		_maxDistance(2),
		_minDistance(2) {
	
}

Detector::Detector(int s) :
		_randomGenerator(s),
		_maxDistance(2),
		_minDistance(2) {
	
}

double Detector::EffectiveArea(double theta, double phi, std::vector<double> delta) {
	double flux = 0;
	double xi, zeta;

	if (fabs(delta[0]) > _maxDistance[0]-_minDistance[0] ||
		fabs(delta[1]) > _maxDistance[0]-_minDistance[0]) {
		return 0;
	}

	for (int i = 0; i < _effectiveAreaIterations; i ++) {
		// Choose random value for rotation angle
		xi   = _randomGenerator.Random()*2.*pi;
		// Choose random value for projected precession angle
		zeta = _randomGenerator.Random()*2.*pi;

		double       zenith = ZenithAngle(theta, phi, xi, zeta, _latitude, _longitude);
		
		if (zenith < 0) {
			i --;
			continue;
		}

		// Calculate projected distance onto detector
		vector<double> dist = ProjectedDistance(delta, theta,phi, xi,zeta, _latitude,_longitude);
		flux +=  SecondDetectionProbability(dist)
				*AngularSensitivity(zenith);
	}

	flux /= (double)_effectiveAreaIterations;
	
	// Return the efective area
	return TotalArea() * flux;
}

double Detector::SecondDetectionProbability(const std::vector<double> &dist) {
	int sum = 0;
	int tot = 0;
	vector<double> r(2);

	if (fabs(dist[0]) > _maxDistance[0]-_minDistance[0] ||
		fabs(dist[1]) > _maxDistance[0]-_minDistance[0]) {
		return 0;
	}

	for (int i = 0; i < _detectionProbabilityIterations; i ++) {
		// Choose random value for position vector
		r = _minDistance;
		r[0] += _randomGenerator.Random()*(_maxDistance[0]-_minDistance[0]-dist[0]);
		r[1] += _randomGenerator.Random()*(_maxDistance[1]-_minDistance[1]-dist[1]);
		
		if (IsDetected(r)) {
			tot ++;
			sum += IsDetected(r+dist);
		}
	}

	if (tot == 0) return 0;
	
//	cout << (double)sum << "\t" << (double)tot << endl;
	return (double)sum/(double)tot;
}
