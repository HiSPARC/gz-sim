#include <general.h>
#include <randomc.h>
#include <GZRunParameters.h>

CRandomMersenne TheGenerator(0);    // make instance of random number generator
#include "RandomGenerator.h"

using namespace std;

vector<double> set(double x, double y, double z) {
    vector<double> vec(3);
    vec[0] = x;
    vec[1] = y;
    vec[2] = z;
    return vec;
}

vector<double> cross(vector<double> A, vector<double> B) {
    vector<double> C(3);
    C[0] = A[1]*B[2] - A[2]*B[1];
    C[1] = A[2]*B[0] - A[0]*B[2];
    C[2] = A[0]*B[1] - A[1]*B[0];
    return C;
}

double abs(vector<double> A) {
    return sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]);
}

int sign(const double &f) {
    if (f > 0) return 1;
    return -1;
}

vector<double> difference(vector<double> A, vector<double> B) {
    A[0] -= B[0];
    A[1] -= B[1];
    A[2] -= B[2];
    return A;
}

double angle(vector<double> A, vector<double> B) {
    return acos((A[0]*B[0] + A[1]*B[1] + A[2]*B[2])/abs(A)/abs(B));
}

vector<double> B(double r, double f, double z) {
    double Br=0., Bf=0., Bz=0.;
    r /= rho0;
    z /= rho0;

    double rr = r*r;
    double zz = z*z;
    double zr52 = pow(zz+rr,-5./2.);
    double zr32 = pow(zz+rr,-3./2.);
    // Dipole component
    Bz += Mdip*(rr-2.*zz)*zr52;
    Br += -3.*Mdip*z*r*zr52;

    // Dynamo component
    Bf += sign(z)*Bphi0/r;

    // Ring current component
    Bz += Brho0*fabs(z)*zr32;
    Br += sign(z)*Brho0*r*zr32;

    double R = .8*Rsun/rho0;
    double RR = R*R;
    // Sunspot component
    for (int i = 0; i < 180; i ++) {
        double r0 = sqrt(rr+RR-2*r*R*cos(i*pi/90.-f));
        double rr0 = r0*r0;
        double zr052 = pow(zz+rr0,-5./2.);
        Bz += Mspot*(rr0-2.*zz)*zr052;
        Br += -3.*Mspot*z*r0*zr052;
    }

    //cout << Br << "\t" << Bf << "\t" << Bz << "\t" << endl;

    vector<double> field(3);
    field[0] = sqrt(Br*Br+Bf*Bf)*cos(f+atan2(Bf,Br));
    field[1] = sqrt(Br*Br+Bf*Bf)*sin(f+atan2(Bf,Br));
    field[2] = Bz;
    return field;
}

double crosssection (double e, int A) {
    //if (e > 3e7)   return  7;
    //if (e > 1.3e7) return 80;
    //return 0.1;

    double e2 = e*e;
    double e2T2 = 6.4e13*e2;
    double e02;
    if (A > 4) {
        e02 = 1.819e15*pow(A,-0.42);
    } else {
        e02 = 8.556e11*pow(A,4.866);
    }
    double sigma = 1.45*A*e2T2/((e2-e02)*(e2-e02)+e2T2);

    if (e > 3e7)
        return max(double(A)/8., sigma);
    else
        return sigma;
}

double integrate (int A, double alpha, double gamma) {
    double total = 0;
    double de = 1e-2;
    for (double epsilon = 0.01; epsilon < 10; epsilon += de) {
        total += crosssection(epsilon*gamma*(1.+cos(alpha)),A)*epsilon*epsilon/(exp(epsilon/kBT)-1.);
        //cout << epsilon*gamma*(1.+cos(alpha)) << "\t" << crosssection(epsilon*gamma*(1.+cos(alpha)),A) << endl;
    }
    return total*1e-31*de;
}

double fragment (int A, vector<double> x, vector<double> v, double e_in) {
    double alpha = pi-angle(x,v);
    double n0 = 7.2e13;
    double r = abs(x);
    double gamma = e_in*eV/(A*1.672622e-27*c*c);

    if (r < Rsun) return 0;
    return (1.+cos(alpha))*n0*(rho0*rho0)/(r*r)*integrate(A,alpha,gamma);
}

vector<double> deflect (GZRunParameters par, vector<double> x, vector<double> v) {
    double q = par.Charge()*eV;//
    double m = par.Mass()*1.672622e-27;//kg
    double E = par.Energy()*eV;//eV
    vector<double> a(3);

    double gamma = E/(m*c*c);
    double t = 0.;
    double dt = 1.0;
    bool within_region = true;
//    bool   entered_region = false;
//    bool   good = true;

    a = set(0, 0, 0);

    long int i = 0;
    while (within_region) {
        t = i*dt;
        i ++;
        // Condition for stopping. Add a small delay to account for deflection.
        if (par.Distance() < c*t/(1+1e-5-1.15907e7/par.Distance())) within_region = false;
        /*if (t > distance/c*1.1) {
            within_region = false;
            //cout << "######";
            cerr << "Aborting deflection run (too long)" << endl;;
            good = false;
        }
        if (abs(x) < Rsun) {
            within_region = false;
            //cout << "######";
            cerr << "Aborting deflection run (solar injection)" << endl;;
            good = false;
        }
        if ((x[0]-rho0)*(x[0]-rho0)+x[1]*x[1]+x[2]*x[2] < Rearth*Rearth)
            within_region = false;*/
        a = cross(v, B(sqrt(x[0]*x[0]+x[1]*x[1]), atan2(x[1],x[0]), x[2]));

        v[0] += dt*q/(gamma*m)*a[0];
        v[1] += dt*q/(gamma*m)*a[1];
        v[2] += dt*q/(gamma*m)*a[2];

        x[0] += dt*v[0];
        x[1] += dt*v[1];
        x[2] += dt*v[2];
        
        // debug: print particle coordinates
        //cout << t << "\t" << x[0] << "\t" << x[1] << "\t" << x[2] << "\t" << "\t" << x2[0] << "\t" << x2[1] << "\t" << x2[2] << endl;
    }
//    double dr = sqrt((x[0]-rho0)*(x[0]-rho0)+x[1]*x[1]+x[2]*x[2]);
//    if (dr > Rearth) cout << "# Earth miss: " << dr << " m" << endl;
    return x;
}

double EnergyFunc(double x) {
    double Elo = 1.e16;
    double Ehi = 1.e20;
    double index = 3.;
    double norm = pow(Ehi,1.-index) - pow(Elo,1.-index);
    
    return pow(pow(Elo,1.-index) + x/norm,1./(1.-index));
}

GZRunParameters SelectRandomParameters() {
    GZRunParameters par;

    // Actual disintegration location will be determined later
    par.SetDistance(RandomUniform()*MaxDistance);
    // Azimuth angle is assumed homogeneous
    par.SetAzimuth (RandomUniform()*2*pi);
    // Zenith angle is assumed homogeneous (corrected for spherical deformation)
    par.SetZenith  (RandomDistribution(&ZenithFunc));
    
    // Choose primary particle properties
    par.SetEnergy (pow(10.,RandomUniform()*4.+16));
    par.SetCharge (RandomInteger(2,92));
    par.SetMass (ChargeToMass(par.Charge()));
    return par;
}

void Run(GZRunParameters par) {
    vector<double> x(3);
    vector<double> v(3);

    x = set(rho0/50., 0, 0);
    v = set(0, c, 0);

    vector<double> f1, f2, distance;
//    double dr = .1*rho0;

//    for (double r = dr; r < 4.*rho0+dr/2.; r += dr) {
        //Convert Earth-centered spherical cooridinates to Heliocenrtic cartesian
        x[0] = par.Distance()*cos(par.Azimuth())*cos(par.Zenith()) + rho0; 
        x[1] = par.Distance()*sin(par.Azimuth())*cos(par.Zenith());
        x[2] = par.Distance()*sin(par.Zenith())+100.;
        
        v[0] = -c*cos(par.Azimuth())*cos(par.Zenith());
        v[1] = -c*sin(par.Azimuth())*cos(par.Zenith());
        v[2] = -c*sin(par.Zenith());
        
        for (int remn = 0; remn <= 1; remn ++) {
            GZRunParameters remnant1 = par;
            GZRunParameters remnant2 = par;

            remnant1.SetMass(1);
            remnant1.SetCharge(remn);
            remnant1.SetEnergy(par.Energy()*(double)remnant1.Mass()/(double)par.Mass());

            remnant2.SetMass(par.Mass()-1);
            remnant2.SetCharge(par.Charge()-remn);
            remnant2.SetEnergy(par.Energy()*(double)remnant2.Mass()/(double)par.Mass());
            
            // First, deflect the smaller remnant
            f1 = deflect(remnant1, x, v);
            if (abs(f1) <= 0.01) continue;
            
            // Then the bigger remnant
            f2 = deflect(remnant2, x, v);
            if (abs(f2) <= 0.01) continue;
            
            // If both arrive at Earth, we can find their distance and angle
            double angle = pi/2.;
            if (f1[2] != f2[2])
                angle = atan2(sqrt((f2[0]-f1[0])*(f2[0]-f1[0])+(f2[1]-f1[1])*(f2[1]-f1[1])),
                              (f1[2]-f2[2]));
            distance = f1 - f2;
            
/*printf("%d        %d        %lg        %lg    %lg        %lg    %lg        %lg        %lg        %lg        %lg        %lg        %lg        %lg        %lg\n", par.Charge(), remn, log10(par.Energy()), par.Zenith()/pi, par.Azimuth()/pi+(par.Azimuth()/pi<0?1:-1), par.Distance()/rho0, abs(distance), angle/pi, fragment((int)par.Mass(),x,v,par.Energy()), f1[0], f1[1], f1[2], f2[0], f2[1], f2[2]);*/
            cout << par.Charge() << "\t"
                 << remn << "\t"
                 << log10(par.Energy()) << "\t"
                 << par.Zenith()/pi << "\t"
                 << par.Azimuth()/pi+(par.Azimuth()/pi<0?1:-1) << "\t"
                 << par.Distance()/rho0 << "\t"
                 << abs(distance) << "\t"
                 << angle/pi << "\t"
                 << fragment((int)par.Mass(), x, v, par.Energy()) << "\t"
                 << f1[0] << "\t" << f1[1] << "\t" << f1[2] << "\t"
                 << f2[0] << "\t" << f2[1] << "\t" << f2[2] << endl;
        }
//    }
//    cout << endl << endl;
}

int main () {
    int32 seed = (int32)time(0);  // random seed
    TheGenerator.RandomInit(seed);  // initialize random number generator

    cout << "# "
         << "Z" << "\t"
         << "Zproduct" << "\t"
         << "log10(E)" << "\t"
         << "theta/pi" << "\t"
         << "phi/pi" << "\t"
         << "r[AU]" << "\t"
         << "distance[m]" << "\t"
         << "angle/pi" << "\t"
         << "prob" << endl;
    for (int i = 0; i < 200000; i ++) {
        GZRunParameters par = SelectRandomParameters();
        Run(par);
    }

    return 0;
}

/*
int main2 () {
    GZRunParameters par;

    // Actual disintegration location will be determined later
    par.SetDistance(RandomUniform()*MaxDistance);
    // Azimuth angle is assumed homogeneous
    par.SetAzimuth (RandomUniform()*2*pi);
    // Zenith angle is assumed homogeneous (corrected for spherical deformation)
    par.SetZenith  (RandomDistribution(&ZenithFunc));
    
    par.SetEnergy  (pow(10.,RandomUniform()*5.+14));
    par.SetCharge  (RandomInteger(2,92));
    par.SetMass    (ChargeToMass(par.Charge()));

    Run(par);

    par.SetZenith  (-par.Zenith());

    Run(par);
    
    return 0;
}
*/
