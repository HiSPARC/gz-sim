/*
 *  DetAuger.h
 *  GZ
 *
 *  Created by Sven Lafèbre on 21-12-07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __DETAUGER_H__
#define __DETAUGER_H__ 1

#include "Detector.h"

class DetAuger : public Detector {

public:
	DetAuger();
	DetAuger(int s);

	bool   IsDetected(const std::vector<double> &r) const;
	double TotalArea() const;
	double AngularSensitivity(double zenith) const;
	double EnergyEfficiency(double energy) const;
	
};

#endif // __DETAUGER_H__
