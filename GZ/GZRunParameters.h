/*
 *  GZRunParameters.h
 *  GZ
 *
 *  Created by Sven Lafèbre on 08-11-07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __GZRUNPARAMETERS_H__
#define __GZRUNPARAMETERS_H__

class GZRunParameters {
	int    _mass;          // In atomic mass units
	int    _charge;        // In elementary charge units
	bool   _proton;
	double _energy;        // In eV
	
	double _distance;
	double _azimuth;
	double _zenith;
	double _prob;
	
public:
	int Mass() const
		{return _mass; }
	int SetMass(const int m)
		{ _mass = m; return _mass; }
	int Charge() const
		{ return _charge; }
	int SetCharge(const int q)
		{ _charge = q; return _charge; }
	int Proton() const
		{ return _proton; }
	int SetProton(const int p)
		{ _proton = p; return _proton; }
	double Energy() const
		{ return _energy; }
	double SetEnergy(const double e)
		{ _energy = e; return _energy; }
	double Probability() const
		{ return _prob; }
	double SetProbability(const double p)
		{ _prob = p; return _prob; }

	double Distance() const
		{ return _distance; }
	double SetDistance(const double d)
		{ _distance = d; return _distance; }
	double Azimuth() const
		{ return _azimuth; }
	double SetAzimuth(const double a)
		{ _azimuth = a; return _azimuth; }
	double Zenith() const
		{ return _zenith; }
	double SetZenith(const double z)
		{ _zenith = z; return _zenith; }

	GZRunParameters();
	GZRunParameters(const GZRunParameters&);
	GZRunParameters& operator= (const GZRunParameters&);
	bool        operator==(const GZRunParameters&) const;
	bool        operator!=(const GZRunParameters&) const;
	bool        operator==(const bool&) const;
	bool        operator!=(const bool&) const;
};

#endif
