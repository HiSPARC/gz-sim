/*
 *  GZRunParameters.cpp
 *  GZ
 *
 *  Created by Sven Lafèbre on 08-11-07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "/user/bdegier/hisparc/lafebre2/GZ/GZ/GZRunParameters.h"

GZRunParameters::GZRunParameters() {

}

GZRunParameters::GZRunParameters(const GZRunParameters& that) {
	SetMass       (that.Mass());
	SetCharge     (that.Charge());
	SetEnergy     (that.Energy());
	SetDistance   (that.Distance());
	SetAzimuth    (that.Azimuth());
	SetZenith     (that.Zenith());
	SetProbability(that.Probability());
}

GZRunParameters& GZRunParameters::operator=(const GZRunParameters& that) {
	if (this != &that) {
		SetMass       (that.Mass());
		SetCharge     (that.Charge());
		SetEnergy     (that.Energy());
		SetDistance   (that.Distance());
		SetAzimuth    (that.Azimuth());
		SetZenith     (that.Zenith());
		SetProbability(that.Probability());
	}
	
	return *this;
}

bool GZRunParameters::operator==(const GZRunParameters& that) const {
	return (
		Mass()        == that.Mass() &&
		Charge()      == that.Charge() &&
		Energy()      == that.Energy() &&
		Distance()    == that.Distance() &&
		Azimuth()     == that.Azimuth() &&
		Zenith()      == that.Zenith() &&
		Probability() == that.Probability());
}

bool GZRunParameters::operator!=(const GZRunParameters& that) const {
	return !(*this == that);
}

bool GZRunParameters::operator==(const bool& that) const {
	return ((bool)Charge() == that);
}

bool GZRunParameters::operator!=(const bool& that) const {
	return (Charge() != that);
}
