/*
 *  DetDense.h
 *  GZ
 *
 *  Created by Sven Lafèbre on 21-12-07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __DETDENSE_H__
#define __DETDENSE_H__ 1

#include "Detector.h"

class DetDense : public Detector {
	double _radius;

public:
	DetDense();
	DetDense(int s);

	bool   IsDetected(const std::vector<double> &r) const;
	double TotalArea() const;
	double AngularSensitivity(double zenith) const;
	double EnergyEfficiency(double energy) const;
	
};

#endif // __DETDENSE_H__
