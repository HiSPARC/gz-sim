/*
 *  id.cpp
 *  GZ
 *
 *  Created by Sven Lafèbre on 22-11-07.
 *
 */

#define AngularResolution  64
#define Ebins               4

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

#include "general.h"
#include "siginthandler.h"
#include "signalhandler.h"
#include "randomc.h"
#include "GZRunParameters.h"
#include "Hist2D.h"
#include "DetHisparc.h"

SigintHandler   TheSigintHandler;	// global sigint handler

CRandomMersenne TheGenerator(0);	// make instance of random number generator
#include "RandomGenerator.h"

using namespace std;

/*!
 *  \brief DetermineAxes
 */
vector<double> DetermineAxes(double xi, double zeta, double latitude, double lon) {
	vector<double> T(3);
	double theta = 0.;

	T[0] = cos(lon+zeta)*cos(theta)*cos(xi)*sin(latitude)
			+ sin(lon+zeta)*sin(latitude)*(sine*sin(theta) - cose*cos(theta)*sin(xi))
			+cos(latitude)*(cose*sin(theta)+cos(theta)*sine*sin(xi));
	T[1] = cos(lon-zeta)*sine*sin(theta)
			+ cos(theta)*(cos(xi)*sin(lon-zeta)-cose*cos(lon-zeta)*sin(xi));
	T[2] = sin(latitude)*(cose*sin(theta) + cos(theta)*sine*sin(xi))
			+ cos(latitude)*(cos(lon+zeta)*cos(theta)*cos(xi)
				+ sin(lon+zeta)*(sine*sin(theta)-cose*cos(theta)*sin(xi)));

	// Return empty vector if we're below ground
	if (T[2] < 0) {
		vector<double> A;
		return A;
	}

	double norm = sqrt(T[0]*T[0] + T[1]*T[1] + T[2]*T[2]);
	cout << "\t" << T[0]/norm << "\t" << T[1]/norm << "\t" << T[2]/norm << endl;
		
	return T;
}

/// \brief Determine zenith and azimuth angle in the detector frame
///
/// \return A two-element vector of the arccosines of the zenith angle and azimuth angle
vector<double> DetermineAzimuthZenith(double theta, double phi, double xi, double zeta, double b, double l) {
	double cl = cos(l+zeta);
	double sl = sin(l+zeta);
	double ce = cos(precession);
	double se = sin(precession);
	double ct = cos(theta);
	double st = sin(theta);
	double cf = cos(phi+xi);
	double sf = sin(phi+xi);
	double cb = cos(b);
	double sb = sin(b);

	vector<double> a(2);
	a[0] = st*(sb*ce + cb*sl*se) + ct*(cf*cb*cl + sf*(sb*se - cb*sl*ce));
	a[1] = st*cl*se - ct*(cf*sl + sf*cl*se);
	
	return a;
}

bool GetParametersFromFile(ifstream* file, vector<double> *d0, GZRunParameters *par) {
	char            line[1024];

	// Read a line with a work package
	do {
		file->getline(line, 1024, '\n');
	} while (file->good() && line[0] == '#');
	if (!file->good()) {
		return 0;
	}

	int    Z, p;
	double E, th, ph, r0, ds, an, pr, x1, y1, z1, x2, y2, z2;
	sscanf(line, "%d %d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
	              &Z,&p, &E,&th,&ph,&r0,&ds,&an,&pr,&x1,&y1,&z1,&x2,&y2,&z2);
				  
	par->SetCharge     (Z				);
	par->SetMass       (AtomicMasses[Z]	);
	par->SetProton     (p				);
	par->SetEnergy     (pow(10.,E)		);
	par->SetDistance   (r0				);
	par->SetAzimuth    (ph*pi			);
	par->SetZenith     (th*pi			);
	par->SetProbability(pr				);

	(*d0)[0] = -sin(ph*pi)*(x1-x2)
	          + cos(ph*pi)*(y1-y2);
	(*d0)[1] =  sin(th*pi)*cos(ph*pi)*(x1-x2)
	          + sin(th*pi)*sin(ph*pi)*(y1-y2)
	          + cos(th*pi)           *(z1-z2);

//	(*d0)[0] = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
//	(*d0)[1] = 0;

	return 1;
}

int energy2bin (double E) {
	double Emin = 1.e16, Emax = 1.e20;
	return int(log10(E/Emin)/log10(Emax/Emin)*Ebins)+1;
}

int main (int argc, char** argv) {
	// install signal handler
	SignalHandler::get_instance().register_handler(SIGINT, &TheSigintHandler);
	
	ifstream file_input;
	string   working_dir = "";
	
	string   file_input_name = working_dir+"GZ_allruns.txt";
//	if (argc > 1) {
//		file_input_name = argv[1];
//	}

	cerr << "Reading file " << file_input_name.c_str() << "..." << endl;
	
	file_input.open(file_input_name.c_str(), ios::in);
	if (!file_input.is_open()) return 1;
 
	vector<double> dist(2);
	GZRunParameters par;
	
	vector<Hist2D> skymap_prob(92);
	vector<Hist2D> skymap_dist(92);
	Hist2D earth_dist         (184,  0.5, 184.5,
							    20, 20.0, 24.0);
	Hist2D earth_map          (100, -1e5, 1e5,
							   100, -1e5, 1e5);
	Hist2D event_rate_limit    (92,  0.5, 92.5,		// Particle species
								80, 16.0, 20.0);	// log(Energy/eV)
	Hist2D event_rate_perfect  (92,  0.5, 92.5,		// Particle species
								80, 16.0, 20.0);	// log(Energy/eV)
	Hist2D event_rate_2ndcheap (92,  0.5, 92.5,		// Particle species
								80, 16.0, 20.0);	// log(Energy/eV)
	Hist2D event_rate_2ndstrict(92,  0.5, 92.5,		// Particle species
								80, 16.0, 20.0);	// log(Energy/eV)

	for (int i = 0; i < 92; i ++) {
		skymap_prob[i].SetBins(AngularResolution+1,   -pi,    pi*65./64.,
							   AngularResolution/2+1, -pi/2., pi/2.*33./32.);
		skymap_dist[i].SetBins(AngularResolution+1,   -pi,    pi*65./64.,
							   AngularResolution/2+1, -pi/2., pi/2.*33./32.);
	}

	DetHisparc detector(0);
	int line = 0;
	while (GetParametersFromFile(&file_input, &dist, &par)) {
		line ++;
		if (!(line % (int)1e5)) cerr << "[" << line << "]" << endl;
		
//		if (par.Energy() < pow(10.,17.9) || par.Energy() > pow(10.,18.1)) continue;

		double d = sqrt(dist*dist);
		double p = par.Probability() * MaxDistance
				*GalacticFlux(par.Charge(), par.Energy()) /* / TotalFlux(par.Energy()) + ExtragalacticFlux(par.Charge(), par.Energy())*/;

		double r = 0.00952171*pow(0.700158,(15.9964-(2.57536*log10( par.Energy()*(par.Mass()-1)/par.Mass() ))));

		double e_proneu = par.Energy()/par.Mass();

		double p_limit     = p;
		double p_perfect   = p;
		double p_2ndcheap  = p;
		double p_2ndstrict = p;

#define DO_DETECTOR

//		---------------------- Cut from here for detector-independent case
#ifdef  DO_DETECTOR
		double effarea     = detector.EffectiveArea(par.Zenith(), par.Azimuth(), dist, r , e_proneu);

		p_limit     = p*(95.*pi*r*r);
		p_perfect   = p*effarea*(95.*pi*r*r);
cout<<effarea<<endl;

		p_2ndcheap  = p
					 *effarea*(95.*pi*r*r)
					 *detector.EnergyEfficiency(10.*par.Energy()/double(par.Mass()))
					 *detector.EnergyEfficiency(par.Energy()*double(par.Mass()-1.)/double(par.Mass()))
				   ;
		p_2ndstrict = p
					 *effarea*(95.*pi*r*r)
					 *detector.EnergyEfficiency(par.Energy()/double(par.Mass()))
				   ;
//		---------------------- Cut up to here
#endif


		// Add event to probability skymap
		skymap_prob[0]
			.Add(par.Azimuth(), par.Zenith(), p_2ndstrict);
		skymap_prob[par.Charge()-1]
			.Add(par.Azimuth(), par.Zenith(), p_2ndstrict);
		
		// Add event to distance skymap
		skymap_dist[0]
			.Add(par.Azimuth(), par.Zenith(), d*p_2ndstrict, p_2ndstrict);
		skymap_dist[1]
			.Add(par.Azimuth(), par.Zenith(), d*p_2ndstrict, p_2ndstrict);

		// Add event to Earth map
		earth_dist
			.Add(2-par.Proton(),              log10(d*par.Energy()/par.Mass()), p);
		earth_dist
			.Add(2*par.Charge()-par.Proton(), log10(d*par.Energy()/par.Mass()), p);

		// Add event to position map
		earth_map
			.Add(dist[0],        dist[1], p);
		
		// Add event rate to event rate arrays
		event_rate_limit
			.Add(1,            log10(par.Energy()), p_limit);
		event_rate_limit
			.Add(par.Charge(), log10(par.Energy()), p_limit);

		// Add event rate to event rate arrays
		event_rate_perfect
			.Add(1,            log10(par.Energy()), p_perfect);
		event_rate_perfect
			.Add(par.Charge(), log10(par.Energy()), p_perfect);

		// Add event rate to event rate arrays
		event_rate_2ndcheap
			.Add(1,            log10(par.Energy()), p_2ndcheap);
		event_rate_2ndcheap
			.Add(par.Charge(), log10(par.Energy()), p_2ndcheap);

		// Add event rate to event rate arrays
		event_rate_2ndstrict
			.Add(1,            log10(par.Energy()), p_2ndstrict);
		event_rate_2ndstrict
			.Add(par.Charge(), log10(par.Energy()), p_2ndstrict);

		// Handle signals
		if (TheSigintHandler.SigintReceived()) {
			cerr << endl << "Analysis aborted by user after line " << line << endl;
			break;	// leave the loop
		}
	}
	file_input.close();

	if (!TheSigintHandler.SigintReceived()) {
		cerr << "Analysis finished after line " << line << endl;
	}
	
#ifdef DO_DETECTOR
	for (int i = 79; i >= 0; i --) {
		for (int j =  0; j < 92; j ++) {
			event_rate_limit    .SetBinValue(j, i, event_rate_limit    .GetBinValue(j, i) * pow(10.,(double)i/20.+16)*(pow(10.,.05)-1.));
			event_rate_perfect  .SetBinValue(j, i, event_rate_perfect  .GetBinValue(j, i) * pow(10.,(double)i/20.+16)*(pow(10.,.05)-1.));
			event_rate_2ndcheap .SetBinValue(j, i, event_rate_2ndcheap .GetBinValue(j, i) * pow(10.,(double)i/20.+16)*(pow(10.,.05)-1.));
			event_rate_2ndstrict.SetBinValue(j, i, event_rate_2ndstrict.GetBinValue(j, i) * pow(10.,(double)i/20.+16)*(pow(10.,.05)-1.));
			if (i < 79) {
				event_rate_limit    .AddBin(j, i, event_rate_limit    .GetBinValue(j, i+1));
				event_rate_perfect  .AddBin(j, i, event_rate_perfect  .GetBinValue(j, i+1));
				event_rate_2ndcheap .AddBin(j, i, event_rate_2ndcheap .GetBinValue(j, i+1));
				event_rate_2ndstrict.AddBin(j, i, event_rate_2ndstrict.GetBinValue(j, i+1));
		} else {
			//	event_rate_limit    .AddBin(j, i, 1.785*event_rate_limit    .GetBinValue(j, i));
			//	event_rate_perfect  .AddBin(j, i, 1.785*event_rate_perfect  .GetBinValue(j, i));
			//	event_rate_2ndcheap .AddBin(j, i, 1.785*event_rate_2ndcheap .GetBinValue(j, i));
			//	event_rate_2ndstrict.AddBin(j, i, 1.785*event_rate_2ndstrict.GetBinValue(j, i));
			}
		}
	}
	//cout << "Final event rate: " << 18043684627.1364*event_rate_2ndstrict.GetBinAverage(0,0) << endl;
	cout << "Final event rate (per year): " << 31556926.*event_rate_2ndstrict.GetBinAverage(0,0) << endl;

#endif
	
	// Write disintegration probability skymap to file
	ofstream file_prob;
	string   file_prob_name = working_dir+"skymap-probability-dense.dat";
	file_prob.open(file_prob_name.c_str(), ios::trunc | ios::out);
	if (!file_prob.is_open()) return 1;
	cerr << "Writing probability skymap... " << file_prob_name << endl;
	for (int i = 0; i < 92; i ++) {
		skymap_prob[i].PrintAverages(&file_prob);
	}
	file_prob.close();

	// Write distance skymap to file
	ofstream file_dist;
	string   file_dist_name = working_dir+"skymap-distance-dense.dat";
	file_dist.open(file_dist_name.c_str(), ios::trunc | ios::out);
	if (!file_dist.is_open()) return 1;
	cerr << "Writing distance skymap...    " << file_dist_name << endl;
	for (int i = 0; i < 92; i ++) {
		skymap_dist[i].PrintAverages(&file_dist);
	}
	file_dist.close();

	// Write Earth distance map to file
	ofstream file_edist;
	string   file_edist_name = working_dir+"earth-distance-dense.dat";
	file_edist.open(file_edist_name.c_str(), ios::trunc | ios::out);
	if (!file_edist.is_open()) return 1;
	cerr << "Writing distance map...       " << file_edist_name << endl;
	earth_dist.PrintValues(&file_edist);
	file_edist.close();

	// Write Earth distance map to file
	ofstream file_map;
	string   file_map_name = working_dir+"earth-map-dense.dat";
	file_map.open(file_map_name.c_str(), ios::trunc | ios::out);
	if (!file_map.is_open()) return 1;
	cerr << "Writing position map...       " << file_map_name << endl;
	earth_map.PrintValues(&file_map);
	file_map.close();

	// Write event rate to file
	ofstream file_rate_limit;
	string   file_rate_limit_name = working_dir+"event-rate-dense-limit.dat";
	file_rate_limit.open(file_rate_limit_name.c_str(), ios::trunc | ios::out);
	if (!file_rate_limit.is_open()) return 1;
	cerr << "Writing event rate...         " << file_rate_limit_name << endl;
	event_rate_limit.PrintAverages(&file_rate_limit);
	file_rate_limit.close();

	// Write event rate to file
	ofstream file_rate_perfect;
	string   file_rate_perfect_name = working_dir+"event-rate-dense-perfect.dat";
	file_rate_perfect.open(file_rate_perfect_name.c_str(), ios::trunc | ios::out);
	if (!file_rate_perfect.is_open()) return 1;
	cerr << "Writing event rate...         " << file_rate_perfect_name << endl;
	event_rate_perfect.PrintAverages(&file_rate_perfect);
	file_rate_perfect.close();

	// Write event rate to file
	ofstream file_rate_2ndcheap;
	string   file_rate_2ndcheap_name = working_dir+"event-rate-dense-2ndcheap.dat";
	file_rate_2ndcheap.open(file_rate_2ndcheap_name.c_str(), ios::trunc | ios::out);
	if (!file_rate_2ndcheap.is_open()) return 1;
	cerr << "Writing event rate...         " << file_rate_2ndcheap_name << endl;
	event_rate_2ndcheap.PrintAverages(&file_rate_2ndcheap);
	file_rate_2ndcheap.close();

	// Write event rate to file
	ofstream file_rate_2ndstrict;
	string   file_rate_2ndstrict_name = working_dir+"event-rate-dense-2ndstrict.dat";
	file_rate_2ndstrict.open(file_rate_2ndstrict_name.c_str(), ios::trunc | ios::out);
	if (!file_rate_2ndstrict.is_open()) return 1;
	cerr << "Writing event rate...         " << file_rate_2ndstrict_name << endl;
	event_rate_2ndstrict.PrintAverages(&file_rate_2ndstrict);
	file_rate_2ndstrict.close();

	return 0;
}
