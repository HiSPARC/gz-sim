/*
 *  DetHisparc.h
 *  GZ
 *
 *  Created by Jeffrey Wouda and Bas de Gier on 22-05-2013.
 *
 */

#ifndef __DETHISPARC_H__
#define __DETHISPARC_H__ 1

#include "Detector.h"

class DetHisparc : public Detector {
	double _coreSize;

	bool InCircle(const vector<double> &pos, const vector<double> &center, double radius) const;


public:
	DetHisparc();
	DetHisparc(int s);

	bool   IsDetected(const std::vector<double> &r) const;
	double TotalArea() const;
	double AngularSensitivity(double zenith) const;
	double EnergyEfficiency(double energy) const;
	
};

#endif // __DETHISPARC_H__
