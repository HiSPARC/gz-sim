/*
 *  DetIceTop.h
 *  GZ
 *
 *  Created by Sven Lafebre on 5-05-09.
 *  Copyright 2009 Pennsylvania State University. All rights reserved.
 *
 */

#ifndef __DETICETOP_H__
#define __DETICETOP_H__ 1

#include "Detector.h"

class DetIceTop : public Detector {
	double     _coreSize;
	
	bool InCircle(const vector<double> &pos, const vector<double> &center, double radius) const;
	
public:
	DetIceTop();
	DetIceTop(int s);
	
	bool   IsDetected(const std::vector<double> &r) const;
	double TotalArea() const;
	double AngularSensitivity(double zenith) const;
	double EnergyEfficiency(double energy) const;
	
};

#endif // __DETICETOP_H__
