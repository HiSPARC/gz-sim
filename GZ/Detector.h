/*
 *  Detector.h
 *  GZ
 *
 *  Created by Sven Lafèbre on 17-12-07.
 *
 */

#ifndef __DETECTOR_H__
#define __DETECTOR_H__ 1

#include <vector>
#include "randomc.h"
#include "general.h"

//#include "general.h"

///
/// \brief    Generic Detector geometry class
///
class Detector {
protected:
    double _latitude;
    double _longitude;
    
    std::vector<double> _minDistance;
    std::vector<double> _maxDistance;
    
    CRandomMersenne _randomGenerator;

    int _effectiveAreaIterations;
    int _detectionProbabilityIterations;
    
public:
    Detector();
    Detector(int s);

    /// \brief    Returns the latitude of the detector
    double Latitude()  const { return _latitude;  }
    /// \brief    Returns the longitude of the detector
    double Longitude() const { return _longitude; }

    /// \brief    Returns the total area covered by the detector
    virtual double TotalArea() const { return 0.; }
    
    /// \brief    Returns whether detection is possible at a certain position
    virtual bool IsDetected(const std::vector<double> &r) const { return false; }
    
    /// \brief    Calculate the effective detector area given a GZ event
    ///
    /// \param    theta    Zenith arrival direction in Solar reference frame
    /// \param    phi        Azimuthal arrival direction in Solar reference frame
    /// \param    delta    Two-element vector containing the EW and NS separation
    ///                    of the GZ elements.
    /// \return    The effective area for the configuration in square kilometers
    double EffectiveArea(double theta, double phi, std::vector<double> delta, double r, double e_proneu);
    
    /// \brief    Returns the probability of detection of the second GZ event
    double SecondDetectionProbability(const std::vector<double> &dist, double r , double e_proneu);
    
    /// \brief    Calculate the angular sensitivity given the zenith angle of a GZ event
    ///
    /// \param    zenith    Zenith angle of the incoming GZ particles
    /// \return    The relative sensitivity of the detector
    virtual double AngularSensitivity(double zenith) const { return 0; }

    /// \brief    Calculate the detector efficiency given a certain energy
    ///
    /// \param    energy    Energy of an incoming GZ remnant
    /// \return    The relative efficiency of the detector
    virtual double EnergyEfficiency(double energy) const { return 0; }

    bool InCircle(const vector<double> &pos, const vector<double> &center, double radius) const;

};

#endif // __DETECTOR_H__
