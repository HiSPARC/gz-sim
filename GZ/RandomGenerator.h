/*
 *  RandomGenerator.h
 *  GZ
 *
 *  Created by Sven Lafèbre on 20-12-07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include </user/bdegier/hisparc/lafebre2/GZ/GZ/randomc.h>

double RandomUniform();

int RandomInteger(int min, int max);

double RandomDistribution(double (*func)(double));

double RandomUniform() {
	// make random floating point numbers in interval from 0 to 1:
	return TheGenerator.Random();
}

int RandomInteger(int min, int max) {
	// make random integer numbers in interval min <= x <= max:
	return TheGenerator.IRandom(min, max);
}

double RandomDistribution(double (*func)(double)) {
	return func(TheGenerator.Random());
}

