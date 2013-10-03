#include </user/bdegier/hisparc/lafebre/GZ_magnetic/general.h>
#include </user/bdegier/hisparc/lafebre/GZ_magnetic/randomc.h>
#include </user/bdegier/hisparc/lafebre/GZ_magnetic/GZRunParameters.h>

CRandomMersenne TheGenerator(0);	// make instance of random number generator
#include "/user/bdegier/hisparc/lafebre/GZ_magnetic/RandomGenerator.h"

using namespace std;

int sign(const double &f) {
	if (f > 0) return 1;
	return -1;
}

int main () {

double z=4.;
double r=0.;

for (int j=0;j<4000; j++){
	double Br=0., Bf=0., Bz=0.;
	r = j/1000.;	

	double rr   = r*r;
	double zz   = z*z;
	double zr52 = pow(zz+rr,-5./2.);
	double zr32 = pow(zz+rr,-3./2.);
	// Dipole component
	Bz += (2.*pow(0.0046491,3.)/2.)*(rr-2.*zz)*zr52;
	Br += -3.*(2.*pow(0.0046491,3.)/2.)*z*r*zr52;

	// Dynamo component
	Bf += sign(z)*(3.5*pow(10.,-5.))/r;

	// Ring current component
	Bz += (-3.5*pow(10.,-5.))*fabs(z)*zr32;
	Br += sign(z)*(-3.5*pow(10.,-5.))*r*zr32;

	double R  = .8*0.0046491;
	double RR = R*R;
	// Sunspot component
	for (int i = 0; i < 180; i ++) {
		double r0  = sqrt(rr+RR-2*r*R*cos(i*pi/90.));
		double rr0 = r0*r0;
		double zr052 = pow(zz+rr0,-5./2.);
		Bz += (1000.*pow(0.1*0.0046491,3.)/2.)*(rr0-2.*zz)*zr052;
		Br += -3.*(1000.*pow(0.1*0.0046491,3.)/2.)*z*r0*zr052;
	}

//In microGauss
	Bz = Bz * 1000000.;
	Br = Br * 1000000.;
	Bf = Bf * 1000000.;

cout << r << "\t" << Br << "\t" << Bz << "\t" << Bf << endl;
	}
return 0;
}
