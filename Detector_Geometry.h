#ifndef DETECTOR_GEOM_H
#define DETECTOR_GEOM_H 1

#include <string>

class Detector_Geometry{

private:

    const double r_min = 0.235;
    const double r_max = 0.325;

    double theta_max;
    double phi_max;

    std::string det_type;

    void load_detector_type();

public:
    Detector_Geometry(std::string);
    ~Detector_Geometry();
    
    bool check_inside(double,double,double);

};


#endif
