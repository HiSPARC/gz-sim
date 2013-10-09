/*
 *  DetHisparc.cpp
 *  GZ
 *
 *  Created by Jeffrey Wouda and Bas de Gier on 13-05-2013.
 *
 */

#include "DetHisparc.h"
#include <iostream>
#include <fstream>

DetHisparc::DetHisparc() : Detector(0) {
    //geographic position in radians of the center of the rectangle
    _latitude = 0.93736207546;
    _longitude = 0.06625163375;

    _minDistance[0] = -4.268e5;
    _minDistance[1] = -2.78e5;
    _maxDistance[0] = 4.268e5;
    _maxDistance[1] = 2.78e5;

    _effectiveAreaIterations = 10;
    _detectionProbabilityIterations = 2500;
}

DetHisparc::DetHisparc(int s) : Detector(s) {
    _latitude = 0.93736207546;
    _longitude = 0.06625163375;

    _minDistance[0] = -4.268e5;
    _minDistance[1] = -2.78e5;
    _maxDistance[0] = 4.268e5;
    _maxDistance[1] = 2.78e5;

    _effectiveAreaIterations = 10;
    _detectionProbabilityIterations = 2500;
}

bool DetHisparc::IsDetected(const std::vector<double> &r) const {

    // EW distances for detector from centre of square
    const double EW[] = {70333.,70527.,74426.,71146.,70486.,69724.,67716.,81273.,70945.,67783.,66936.,66022.,67017.,66953.,53783.,58294.,55052.,58648.,56957.,57404.,60994.,76033.,75963.,76151.,76246.,75851.,76088.,76055.,76332.,83797.,90168.,87760.,89049.,89117.,87998.,85638.,136410.,136839.,135549.,133060.,134478.,136704.,119661.,119776.,43615.,45386.,56559.,55831.,55753.,48456.,39405.,23892.,23985.,24041.,58133.,56801.,56445.,2817.,179730.,183996.,183133.,201552.,201528.,201488.,197571.,188495.,154569.,193792.,111364.,108662.,111418.,111714.,94489.,111343.,117172.,111174.,111386.,84777.,83325.,69988.,83161.,143620.,126485.,125943.,124312.,-318728.,-420953.,-421691.,-421116.,-422111.,-419029.,-405019.,-405336.,-365510.,421307.,421516.,421553.,
    };

    // NS distances for detector from centre of square
    const double NS[] = {-150766.,-153877.,-145861.,-150651.,-152783.,-157824.,-149489.,-155370.,-150943.,-137035.,-139826.,-137005.,-139862.,-139911.,-149500.,-155305.,-151162.,-128434.,-134830.,-129572.,-118737.,-150220.,-150291.,-150183.,-150080.,-150073.,-150080.,-150212.,-150373.,-117809.,-180099.,-180307.,-180666.,-180462.,-180669.,-181435.,-209453.,-209224.,-206771.,-207914.,-211301.,-210226.,-202359.,-152736.,-171001.,-170544.,-209547.,-208196.,-210124.,-207173.,-210009.,-217166.,-216943.,-217027.,-173932.,-174751.,-173512.,-265029.,-50805.,-58115.,-52188.,-163255.,-163307.,-163304.,-159371.,-148783.,-159266.,-172169.,-251022.,-247848.,-255286.,-251303.,-259664.,-247331.,-253444.,-251367.,-250785.,-240755.,-238553.,-230390.,-241855.,-265336.,-272519.,-271713.,-271677.,92426.,-250034.,-246558.,-250072.,-250009.,-243841.,-259606.,-256884.,-238959.,273546.,273600.,273397.,
    };

    vector<double> pos(2, 0.);

    /*
    //This part is not relevant since seconddetectionprob is calculated in detector.cpp
    for (int i = 0; i < 97; i ++) {
        pos[0] = EW[i];
        pos[1] = NS[i];
        if (InCircle(r, pos, 400.)) {
    return true;}
    }
        
    return false;
    */
}

double DetHisparc::TotalArea() const {
    //total area of the HiSPARC array is 556 km x 829.2 km
    return 2.38327718545e11;
}

double DetHisparc::AngularSensitivity(double zenith) const {
    //maximum value for the zenith angle
    return zenith < pi/4.;
}

double DetHisparc::EnergyEfficiency(double energy) const {
    //minimum energy which HiSPARC is able to detect (10^15 eV)
    return erf((log10(energy)-15.)*1.7)/2.+.5;
}

bool DetHisparc::InCircle(const vector<double> &pos,
                          const vector<double> &center,
                          double radius) const {
    return (pos-center)*(pos-center) <= radius*radius;    
}
