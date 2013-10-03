/*
 *  DetLofar.h
 *  GZ
 *
 *  Created by Sven Lafèbre on 21-12-07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __DETLOFAR_H__
#define __DETLOFAR_H__ 1

#include "Detector.h"

class DetLofar : public Detector {
	double _coreSize;
	double _remoteSize;

	bool InCircle(const vector<double> &pos, const vector<double> &center, double radius) const;

public:
	DetLofar();
	DetLofar(int s);

	bool   IsDetected(const std::vector<double> &r) const;
	double TotalArea() const;
	double AngularSensitivity(double zenith) const;
	double EnergyEfficiency(double energy) const;
	
};

#endif // __DETLOFAR_H__
