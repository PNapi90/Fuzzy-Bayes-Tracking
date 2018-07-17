#include "Detector_Geometry.h"
#include <cmath>

using namespace std;

Detector_Geometry::Detector_Geometry(string det_type){
    this->det_type = det_type;
    load_detector_type();
}

Detector_Geometry::~Detector_Geometry(){}


void Detector_Geometry::load_detector_type(){
    if(det_type == "LNL"){
        theta_max = 37.5*M_PI/180.;
        phi_max = 2.*M_PI;
    }
    else if(det_type == "GANIL"){
        return;
    }
    else if(det_type == "GSI"){
        return;
    }
}


bool Detector_Geometry::check_inside(double r, double th, double phi){
    if(r > r_min && r < r_max){
        if(th < theta_max){
            if(phi < phi_max) return true;
        }
    }
    return false;
}
